# Choice module: conditional logit with latent positions.
# par = (zeta[1..S], vec(U) [n*D], vec(V) [n*D]).

choice_unpack <- function(par, n, D, S) {
  list(zeta = if (S > 0) par[seq_len(S)] else numeric(0),
       U    = matrix(par[S + seq_len(n * D)], n, D),
       V    = matrix(par[S + n * D + seq_len(n * D)], n, D))
}

# Linear predictors eta[m, j] = zeta's(i_m, j; t_m) + u_{i_m}'v_j (eta = -Inf
# for the sender), choice probabilities pi and the log normalising constants.
choice_eta_pi <- function(zeta, U, V, sl, senders, cnames) {
  M <- length(senders); n <- nrow(U)
  eta <- matrix(0, M, n)
  for (s in seq_along(cnames)) eta <- eta + zeta[s] * sl[[cnames[s]]]
  if (ncol(U) > 0) eta <- eta + U[senders, , drop = FALSE] %*% t(V)
  eta[cbind(seq_len(M), senders)] <- -Inf
  # row maximum without an R-level loop over the M events; ties.method =
  # "first" keeps the random number stream untouched
  row_max <- eta[cbind(seq_len(M), max.col(eta, ties.method = "first"))]
  exp_eta <- exp(eta - row_max)
  row_sum <- rowSums(exp_eta)
  list(eta = eta, pi = exp_eta / row_sum, log_psi = row_max + log(row_sum))
}

# Negative log-posterior of the choice module (Normal priors).
choice_negloglik <- function(par, n, D, S, sl, senders, receivers, cnames,
                              sigma2_zeta, sigma2_u, sigma2_v) {
  pk <- choice_unpack(par, n, D, S)
  ep <- choice_eta_pi(pk$zeta, pk$U, pk$V, sl, senders, cnames)
  M  <- length(senders)
  loglik  <- sum(ep$eta[cbind(seq_len(M), receivers)] - ep$log_psi)
  penalty <- 0.5 * sum(pk$zeta^2) / sigma2_zeta +
             0.5 * sum(pk$U^2)    / sigma2_u    +
             0.5 * sum(pk$V^2)    / sigma2_v
  -(loglik - penalty)
}

# Its gradient.
choice_neggrad <- function(par, n, D, S, sl, senders, receivers, cnames,
                            sigma2_zeta, sigma2_u, sigma2_v) {
  pk <- choice_unpack(par, n, D, S)
  ep <- choice_eta_pi(pk$zeta, pk$U, pk$V, sl, senders, cnames)
  M  <- length(senders)

  g_zeta <- numeric(S)
  for (s in seq_along(cnames)) {
    g_zeta[s] <- sum(sl[[cnames[s]]][cbind(seq_len(M), receivers)]) -
                 sum(ep$pi * sl[[cnames[s]]])
  }
  g_zeta <- g_zeta - pk$zeta / sigma2_zeta
  if (D == 0) return(-g_zeta)   # no latent positions (reference model M0)

  # vectorized replacement for a per-actor loop with which(senders==i) /
  # which(senders!=i) inside it -- that was O(n*M) redundant scanning per
  # gradient call (found to roughly double per-call cost at n=167,
  # M~25000). Mathematically exact rewrite (verified identical to machine
  # precision against the loop version), not an approximation.
  full_rowsum <- function(x, group, n) {
    out <- matrix(0, n, ncol(x))
    ug  <- sort(unique(group))
    out[ug, ] <- rowsum(x, group, reorder = TRUE)
    out
  }

  # g_U: positive = sum over i's own outgoing events of V[receiver,];
  #      negative = sum over i's own outgoing events of (pi[m,] %*% V)
  g_U_pos <- full_rowsum(pk$V[receivers, , drop = FALSE], senders, n)
  g_U_neg <- full_rowsum(ep$pi %*% pk$V, senders, n)
  g_U <- g_U_pos - g_U_neg - pk$U / sigma2_u

  # g_V: positive = sum over events where i IS receiver, of U[sender,];
  #      negative = sum over ALL events of pi[m,i]*U[sender,] -- pi[m,sender_m]
  #      is exactly 0 by construction (self-choice excluded), so no
  #      explicit sender!=i filtering is needed
  g_V_pos <- full_rowsum(pk$U[senders, , drop = FALSE], receivers, n)
  g_V_neg <- t(ep$pi) %*% pk$U[senders, , drop = FALSE]
  g_V <- g_V_pos - g_V_neg - pk$V / sigma2_v

  -c(g_zeta, as.vector(g_U), as.vector(g_V))
}

# Log-posterior of the choice module and its gradient from one evaluation of
# the linear predictors (used by the sampler).
choice_val_grad <- function(par, n, D, S, sl, senders, receivers, cnames,
                            sigma2_zeta, sigma2_u, sigma2_v) {
  pk <- choice_unpack(par, n, D, S)
  ep <- choice_eta_pi(pk$zeta, pk$U, pk$V, sl, senders, cnames)
  M  <- length(senders)
  idx <- cbind(seq_len(M), receivers)
  value <- sum(ep$eta[idx] - ep$log_psi) -
    0.5 * sum(pk$zeta^2) / sigma2_zeta - 0.5 * sum(pk$U^2) / sigma2_u - 0.5 * sum(pk$V^2) / sigma2_v
  g_zeta <- numeric(S)
  for (s in seq_along(cnames)) {
    g_zeta[s] <- sum(sl[[cnames[s]]][idx]) - sum(ep$pi * sl[[cnames[s]]])
  }
  g_zeta <- g_zeta - pk$zeta / sigma2_zeta
  if (D == 0) return(list(value = value, gradient = g_zeta))
  Y <- ep$pi; Y[idx] <- Y[idx] - 1                       # pi - y
  g_U <- -rowsum(Y %*% pk$V, senders, reorder = TRUE)    # sum over events sent by i of (y - pi) V
  g_U_full <- matrix(0, n, D); g_U_full[sort(unique(senders)), ] <- g_U
  g_V <- -crossprod(Y, pk$U[senders, , drop = FALSE])    # sum over events of (y - pi)_j u_sender
  list(value = value,
       gradient = c(g_zeta, as.vector(g_U_full - pk$U / sigma2_u), as.vector(g_V - pk$V / sigma2_v)))
}

# Starting configuration from the singular value decomposition of the
# row-centred matrix of log(1 + number of events from i to j). Row-centring
# removes the sender effect, which cancels in the choice model.
spectral_start_choice <- function(events_df, n, D) {
  tab <- table(factor(events_df$sender, levels = seq_len(n)),
               factor(events_df$receiver, levels = seq_len(n)))
  X <- log1p(matrix(as.numeric(tab), n, n))
  diag(X) <- NA
  X <- X - rowMeans(X, na.rm = TRUE)
  diag(X) <- 0
  sv <- svd(X, nu = D, nv = D)
  sc <- sqrt(sv$d[seq_len(D)])
  list(U = sweep(sv$u, 2, sc, "*"), V = sweep(sv$v, 2, sc, "*"))
}

# Posterior mode of the choice module by BFGS from several starts (n_starts
# random starts and, with spectral_start = TRUE, one spectral start). The
# objective is not concave in (U, V) and has local modes; the best start is
# kept and all of them are returned ($starts, $start_par), with their
# convergence codes (1 = iteration limit reached). The covariance of the
# Laplace approximation floors the eigenvalues of the Hessian at half the
# prior precision of the positions before inversion, so that the flat
# directions (rotations, for D >= 2) get a large but finite variance.
fit_choice_model <- function(events_df, n, D, sl,
                              sigma2_zeta = 4, sigma2_u = 1, sigma2_v = 1,
                              seed = 1, hessian = TRUE, n_starts = 1,
                              spectral_start = FALSE, maxit = 500) {
  cnames <- choice_stat_names(sl)
  S <- length(cnames)
  senders   <- events_df$sender
  receivers <- events_df$receiver

  if (S + 2 * n * D == 0) {
    # no choice statistics and no positions: receivers are chosen uniformly
    return(list(zeta = numeric(0), U = matrix(0, n, 0), V = matrix(0, n, 0), se = NULL,
                vcov = matrix(0, 0, 0), hessian = matrix(0, 0, 0),
                value = length(senders) * log(n - 1), convergence = 0, counts = c(0, 0),
                n_starts = 0, start_values = numeric(0), n_floored = 0))
  }

  labels <- as.character(seq_len(n_starts))
  if (spectral_start && D > 0) labels <- c(labels, "spectral")
  fit <- NULL
  start_par <- matrix(NA_real_, length(labels), S + 2 * n * D, dimnames = list(labels, NULL))
  starts <- data.frame(start = labels, value = NA_real_, convergence = NA_integer_,
                       gradient_calls = NA_integer_, stringsAsFactors = FALSE)
  for (k in seq_along(labels)) {
    if (labels[k] == "spectral") {
      sp <- spectral_start_choice(events_df, n, D)
      par0 <- c(rep(0, S), as.vector(sp$U), as.vector(sp$V))
    } else {
      set.seed(seed + k - 1)
      par0 <- c(rep(0, S), rnorm(n * D, 0, 0.1), rnorm(n * D, 0, 0.1))
    }
    fit_k <- optim(par0, choice_negloglik, gr = choice_neggrad,
                   n = n, D = D, S = S, sl = sl, senders = senders, receivers = receivers,
                   cnames = cnames, sigma2_zeta = sigma2_zeta, sigma2_u = sigma2_u, sigma2_v = sigma2_v,
                   method = "BFGS", hessian = FALSE,
                   control = list(maxit = maxit, reltol = 1e-10))
    starts$value[k] <- fit_k$value; starts$convergence[k] <- fit_k$convergence
    starts$gradient_calls[k] <- fit_k$counts[2]
    start_par[k, ] <- fit_k$par
    if (is.null(fit) || fit_k$value < fit$value) fit <- fit_k
  }
  start_values <- setNames(starts$value, labels)
  if (any(starts$convergence != 0)) {
    warning(sprintf("fit_choice_model: BFGS stopped at the iteration limit (maxit = %d) for start(s) %s; increase maxit.",
                    maxit, paste(labels[starts$convergence != 0], collapse = ", ")))
  }

  pk <- choice_unpack(fit$par, n, D, S)
  names(pk$zeta) <- cnames

  se <- NULL
  vcov <- NULL
  H <- NULL
  n_floored <- NA
  if (hessian) {
    H <- optimHess(fit$par, choice_negloglik, gr = choice_neggrad,
                   n = n, D = D, S = S, sl = sl, senders = senders, receivers = receivers,
                   cnames = cnames, sigma2_zeta = sigma2_zeta, sigma2_u = sigma2_u, sigma2_v = sigma2_v)
    # floor the Hessian's eigenvalues, then invert: for D >= 2 the Hessian is
    # singular along the rotations of the latent space
    hc <- hessian_to_cov(H, abs_floor = if (D > 0) 0.5 * min(1 / sigma2_u, 1 / sigma2_v) else 0)
    vcov <- hc$cov
    n_floored <- hc$n_floored
    se_all <- sqrt(pmax(diag(vcov), 0))
    se <- list(zeta = se_all[seq_len(S)],
               U    = matrix(se_all[S + seq_len(n * D)], n, D),
               V    = matrix(se_all[S + n * D + seq_len(n * D)], n, D))
  }

  list(zeta = pk$zeta, U = pk$U, V = pk$V, se = se, vcov = vcov, hessian = H,
       value = fit$value, convergence = fit$convergence, counts = fit$counts,
       n_starts = n_starts, start_values = start_values, starts = starts, start_par = start_par,
       n_floored = n_floored)
}

# Orthogonal Procrustes alignment of every draw to a reference configuration
# (U and V jointly, since they are subject to the same transformation).
align_choice_draws <- function(U_draws, V_draws, U_ref, V_ref) {
  n <- dim(U_draws)[2]; D <- dim(U_draws)[3]
  if (D == 0) return(list(U = U_draws, V = V_draws))
  for (s in seq_len(dim(U_draws)[1])) {
    U_s <- matrix(U_draws[s, , ], n, D)
    V_s <- matrix(V_draws[s, , ], n, D)
    M_proc <- t(U_s) %*% U_ref + t(V_s) %*% V_ref
    if (D > 1) {
      sv <- svd(M_proc)
      R_proc <- sv$u %*% t(sv$v)
    } else {
      R_proc <- matrix(if (M_proc[1, 1] < 0) -1 else 1, 1, 1)
    }
    U_draws[s, , ] <- U_s %*% R_proc
    V_draws[s, , ] <- V_s %*% R_proc
  }
  list(U = U_draws, V = V_draws)
}

# Rotation of a choice-module parameter vector to the orientation of a
# reference configuration, for the aligned MALA sampler.
make_choice_aligner <- function(n, D, S, U_ref, V_ref) {
  iu <- S + seq_len(n * D); iv <- S + n * D + seq_len(n * D)
  function(par) {
    Mp <- crossprod(matrix(par[iu], n, D), U_ref) + crossprod(matrix(par[iv], n, D), V_ref)
    W <- if (D > 1) { sv <- svd(Mp); sv$u %*% t(sv$v) } else matrix(if (Mp[1, 1] < 0) -1 else 1, 1, 1)
    rot <- function(p) {
      p[iu] <- as.vector(matrix(p[iu], n, D) %*% W)
      p[iv] <- as.vector(matrix(p[iv], n, D) %*% W)
      p
    }
    list(par = rot(par), rot = rot)
  }
}

#' Stage 1: the choice module
#'
#' Estimates the choice sub-process, with the coefficients \eqn{\zeta} of the
#' choice statistics and the latent positions, from the choice of receivers
#' alone. The posterior mode is found by BFGS from several starting
#' configurations; the posterior is then simulated by a Metropolis-adjusted
#' Langevin algorithm (MALA) preconditioned with the covariance of the Laplace
#' approximation at the mode (`inference = "mcmc"`), or approximated by
#' draws from the Laplace approximation (`inference = "map"`). The draws are
#' aligned to the mode by orthogonal Procrustes rotations, and with
#' `align_state = TRUE` the state of the chain is kept in the orientation of
#' the mode, which keeps the preconditioner appropriate for \eqn{D \ge 2}.
#'
#' Usually called through [lsdynam_fit()]; call it directly to fit several
#' rate modules on the same positions (argument `choice` of
#' [lsdynam_fit()] and [lsdynam_rate_stage()]).
#'
#' @inheritParams lsdynam_fit
#' @param sl statistics computed by [compute_stats()].
#' @param n_draws_laplace number of draws from the Laplace approximation when
#'   `inference = "map"`.
#' @param align_state keep the state of the MALA chain in the orientation of
#'   the mode (recommended).
#'
#' @return A list with the draws of `zeta` (matrix), `U` and `V` (arrays of
#'   dimension (draws, \eqn{n}, \eqn{D}), aligned to the mode), the posterior
#'   means `U_mean` and `V_mean` (the mode for `inference = "map"`), the
#'   posterior mode and Laplace approximation `map` (with the solution and
#'   convergence code of every start in `map$starts` and `map$start_par`),
#'   the acceptance rate `accept_rate`, the final step size `eps`, the trace
#'   of the log-posterior `log_post`, `D`, `n` and the computing time
#'   `seconds`.
#' @seealso [lsdynam_fit()], [mode_stability_check()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 300, seed = 2)
#' sl <- compute_stats(ev, n = 15, stats = c("rate_indeg", "choice_recip"))
#' ch <- lsdynam_choice_stage(ev, n = 15, D = 1, sl = sl, inference = "map")
#' ch$map$zeta
#' @export
lsdynam_choice_stage <- function(events_df, n, D, sl, priors = lsdynam_priors(),
                                 inference = c("mcmc", "map"), sampler = c("mala", "rwm"),
                                 n_iter = 6000, burn_in = 1000, thin = 5, n_starts = 1,
                                 n_draws_laplace = 500, align_state = TRUE,
                                 spectral_start = TRUE, maxit = 5000, seed = 1, verbose = FALSE) {
  inference <- match.arg(inference); sampler <- match.arg(sampler)
  cnames <- choice_stat_names(sl); S <- length(cnames)
  senders <- events_df$sender; receivers <- events_df$receiver
  t0 <- Sys.time()
  map_fit <- fit_choice_model(events_df, n, D, sl, sigma2_zeta = priors$sigma2_zeta,
                              sigma2_u = priors$sigma2_u, sigma2_v = priors$sigma2_v,
                              seed = seed, n_starts = n_starts,
                              spectral_start = spectral_start, maxit = maxit)
  if (verbose) {
    cat(sprintf("  choice MAP: %.1f s, convergence %d, %s Hessian eigenvalue(s) floored\n",
                as.numeric(Sys.time() - t0, units = "secs"), map_fit$convergence, format(map_fit$n_floored)))
    if (!is.null(map_fit$starts) && nrow(map_fit$starts) > 1) {
      st <- map_fit$starts
      cat(sprintf("  starts: objective minus best = %s | convergence codes %s\n",
                  paste(sprintf("%s: %.2f", st$start, st$value - min(st$value)), collapse = ", "),
                  paste(st$convergence, collapse = " ")))
    }
  }
  d <- S + 2 * n * D
  par_map <- c(map_fit$zeta, as.vector(map_fit$U), as.vector(map_fit$V))
  accept <- NA; eps <- NA; lp_trace <- NULL
  if (d == 0) {
    samples <- matrix(0, n_draws_laplace, 0)
  } else if (inference == "map") {
    set.seed(seed + 500)
    Lc <- t(chol(make_pd(map_fit$vcov)))
    samples <- matrix(par_map, n_draws_laplace, d, byrow = TRUE) +
      t(Lc %*% matrix(rnorm(d * n_draws_laplace), d, n_draws_laplace))
  } else if (sampler == "mala") {
    fn <- function(p) choice_val_grad(p, n, D, S, sl, senders, receivers, cnames,
                                      priors$sigma2_zeta, priors$sigma2_u, priors$sigma2_v)
    mc <- mala_sample(fn, par_map, map_fit$vcov, n_iter = n_iter, burn_in = burn_in, thin = thin,
                      align = if (D > 0 && align_state) make_choice_aligner(n, D, S, map_fit$U, map_fit$V) else NULL,
                      seed = seed + 1000, verbose = verbose)
    samples <- mc$samples; accept <- mc$accept_rate; eps <- mc$eps; lp_trace <- mc$log_post
  } else {
    log_post <- function(p) -choice_negloglik(p, n, D, S, sl, senders, receivers, cnames,
                                              priors$sigma2_zeta, priors$sigma2_u, priors$sigma2_v)
    mc <- rwm_sample(log_post, par_map, map_fit$vcov, n_iter = n_iter, burn_in = burn_in, thin = thin,
                     seed = seed + 1000)
    samples <- mc$samples; accept <- mc$accept_rate
  }
  ns <- nrow(samples)
  zeta_draws <- samples[, seq_len(S), drop = FALSE]; colnames(zeta_draws) <- cnames
  U_draws <- array(samples[, S + seq_len(n * D)], dim = c(ns, n, D))
  V_draws <- array(samples[, S + n * D + seq_len(n * D)], dim = c(ns, n, D))
  al <- align_choice_draws(U_draws, V_draws, map_fit$U, map_fit$V)
  U_mean <- if (D > 0) apply(al$U, c(2, 3), mean) else matrix(0, n, 0)
  V_mean <- if (D > 0) apply(al$V, c(2, 3), mean) else matrix(0, n, 0)
  if (inference == "map") { U_mean <- map_fit$U; V_mean <- map_fit$V }
  list(zeta = zeta_draws, U = al$U, V = al$V, map = map_fit, U_mean = U_mean, V_mean = V_mean,
       accept_rate = accept, eps = eps, log_post = lp_trace, inference = inference,
       sampler = if (inference == "mcmc") sampler else "laplace", D = D, n = n,
       seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
