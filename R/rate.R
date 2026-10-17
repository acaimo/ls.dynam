# Rate module. The rate statistics of actor i change only when i sends or
# receives, so the exposure sum_m Delta_m exp{xi's(i;t_m) + a_i} over all M
# events collapses to one term per spell during which the statistics of i are
# constant, at a cost of order M + n per evaluation instead of M * n.
#   rows       : event indices to use (default all); used for held-out scores
#   rate_event : logical, TRUE for events that count as a sending decision
#                (default all)
compress_rate_stats <- function(sl, senders, n, rows = NULL, rate_event = NULL) {
  rnames <- rate_stat_names(sl); R <- length(rnames)
  M <- length(sl$delta_t)
  if (is.null(rows)) rows <- seq_len(M)
  if (is.null(rate_event)) rate_event <- rep(TRUE, M)
  dt <- sl$delta_t[rows]
  snd <- senders[rows]; ev <- rate_event[rows]
  N <- tabulate(snd[ev], nbins = n)
  if (R == 0) {
    return(list(actor = seq_len(n), S = matrix(0, n, 0), expo = rep(sum(dt), n), N = N,
                S_event = numeric(0), R = 0, rnames = rnames, n = n,
                expo_actor = rep(sum(dt), n), n_events = sum(ev)))
  }
  S_event <- vapply(rnames, function(rn) sum(sl[[rn]][cbind(rows, snd)][ev]), numeric(1))
  actor_l <- vector("list", n); S_l <- vector("list", n); expo_l <- vector("list", n)
  Mr <- length(rows)
  for (i in seq_len(n)) {
    X <- vapply(rnames, function(rn) sl[[rn]][rows, i], numeric(Mr))
    X <- matrix(X, nrow = Mr)
    chg <- if (Mr > 1) c(TRUE, rowSums(abs(X[-1, , drop = FALSE] - X[-Mr, , drop = FALSE])) > 0) else TRUE
    grp <- cumsum(chg)
    expo_l[[i]]  <- as.numeric(rowsum(dt, grp, reorder = TRUE))
    S_l[[i]]     <- X[chg, , drop = FALSE]
    actor_l[[i]] <- rep(i, sum(chg))
  }
  S <- do.call(rbind, S_l); colnames(S) <- rnames
  list(actor = unlist(actor_l), S = S, expo = unlist(expo_l), N = N, S_event = S_event,
       R = R, rnames = rnames, n = n, expo_actor = rep(sum(dt), n), n_events = sum(ev))
}

# Log-posterior of the rate module, with all normalising constants (needed
# for the marginal likelihood), its gradient and, with deriv = 2, the Hessian
# of the NEGATIVE log-posterior. par = (xi[1..R], a[1..n], mu_a, kappa[1..K]);
# B is the n x K link basis (K = 0: no link); nu < Inf gives Student-t
# activity effects.
rate_mod_eval <- function(par, cr, B, hp, sigma2_a, nu = Inf, deriv = 1) {
  n <- cr$n; R <- cr$R; K <- ncol(B)
  xi <- par[seq_len(R)]; a <- par[R + seq_len(n)]; mu <- par[R + n + 1]
  kappa <- par[R + n + 1 + seq_len(K)]

  lin <- a[cr$actor]
  if (R > 0) lin <- lin + as.numeric(cr$S %*% xi)
  Lam <- cr$expo * exp(lin)
  loglik <- sum(cr$N * a) - sum(Lam)
  if (R > 0) loglik <- loglik + sum(cr$S_event * xi)

  fitted <- rep(mu, n)
  if (K > 0) fitted <- fitted + as.numeric(B %*% kappa)
  resid <- a - fitted
  if (is.infinite(nu)) {
    lp_a <- -0.5 * n * log(2 * pi * sigma2_a) - 0.5 * sum(resid^2) / sigma2_a
    psi  <- resid / sigma2_a
    w    <- rep(1 / sigma2_a, n)
  } else {
    lp_a <- n * (lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi * sigma2_a)) -
            0.5 * (nu + 1) * sum(log1p(resid^2 / (nu * sigma2_a)))
    psi  <- (nu + 1) * resid / (nu * sigma2_a + resid^2)
    w    <- (nu + 1) * (nu * sigma2_a - resid^2) / (nu * sigma2_a + resid^2)^2
  }
  lp <- lp_a - 0.5 * log(2 * pi * hp$sigma2_mu) - 0.5 * mu^2 / hp$sigma2_mu
  if (R > 0) lp <- lp - 0.5 * R * log(2 * pi * hp$sigma2_xi) - 0.5 * sum(xi^2) / hp$sigma2_xi
  if (K > 0) lp <- lp - 0.5 * K * log(2 * pi * hp$sigma2_kappa) - 0.5 * sum(kappa^2) / hp$sigma2_kappa
  value <- loglik + lp
  if (!is.finite(value)) return(list(value = -Inf, gradient = rep(0, length(par)), loglik = loglik))
  if (deriv == 0) return(list(value = value, loglik = loglik))

  Lam_actor <- as.numeric(rowsum(Lam, cr$actor, reorder = TRUE))
  g_xi <- if (R > 0) cr$S_event - as.numeric(crossprod(cr$S, Lam)) - xi / hp$sigma2_xi else numeric(0)
  g_a  <- cr$N - Lam_actor - psi
  g_mu <- sum(psi) - mu / hp$sigma2_mu
  g_k  <- if (K > 0) as.numeric(crossprod(B, psi)) - kappa / hp$sigma2_kappa else numeric(0)
  out <- list(value = value, gradient = c(g_xi, g_a, g_mu, g_k), loglik = loglik)

  if (deriv >= 2) {
    d <- R + n + 1 + K
    H <- matrix(0, d, d)
    ia <- R + seq_len(n); im <- R + n + 1
    if (R > 0) {
      ix <- seq_len(R)
      LS <- Lam * cr$S
      H[ix, ix] <- crossprod(cr$S, LS) + diag(1 / hp$sigma2_xi, R)
      SL <- rowsum(LS, cr$actor, reorder = TRUE)           # n x R
      H[ia, ix] <- SL; H[ix, ia] <- t(SL)
    }
    H[cbind(ia, ia)] <- Lam_actor + w
    H[ia, im] <- -w; H[im, ia] <- -w
    H[im, im] <- sum(w) + 1 / hp$sigma2_mu
    if (K > 0) {
      ik <- R + n + 1 + seq_len(K)
      wB <- w * B
      H[ia, ik] <- -wB; H[ik, ia] <- t(-wB)
      H[im, ik] <- colSums(wB); H[ik, im] <- colSums(wB)
      H[ik, ik] <- crossprod(B, wB) + diag(1 / hp$sigma2_kappa, K)
    }
    out$hessian <- H
  }
  out
}

.rate_unpack_par <- function(par, cr, B) {
  n <- cr$n; R <- cr$R; K <- ncol(B)
  xi <- par[seq_len(R)]; names(xi) <- cr$rnames
  kap <- par[R + n + 1 + seq_len(K)]; names(kap) <- colnames(B)
  list(xi = xi, a = par[R + seq_len(n)], mu_a = par[R + n + 1], kappa_vec = kap)
}

# Posterior mode of the rate module by Newton's method with step halving, and
# the Laplace approximation to the log marginal likelihood. For the Gaussian
# link and fixed sigma_a the log-posterior is strictly concave.
rate_mod_fit <- function(cr, B, hp, sigma2_a, nu = Inf, start = NULL, maxit = 200, tol = 1e-8) {
  n <- cr$n; R <- cr$R; K <- ncol(B)
  d <- R + n + 1 + K
  if (is.null(start)) {
    a0 <- log((cr$N + 0.5) / pmax(cr$expo_actor, 1e-300))
    start <- c(rep(0, R), a0, mean(a0), rep(0, K))
  }
  par <- start
  ev <- rate_mod_eval(par, cr, B, hp, sigma2_a, nu, deriv = 2)
  if (!is.finite(ev$value)) stop("rate_mod_fit: non-finite log-posterior at the starting point")
  converged <- FALSE
  for (it in seq_len(maxit)) {
    step <- if (is.infinite(nu)) {
      tryCatch(solve(ev$hessian, ev$gradient), error = function(e) NULL)
    } else NULL
    if (is.null(step)) {
      # Student-t link (Hessian not necessarily positive definite) or a
      # numerically singular system: use the regularised inverse
      step <- as.numeric(hessian_to_cov(ev$hessian, rel_floor = 1e-6)$cov %*% ev$gradient)
    }
    slope <- sum(ev$gradient * step)
    if (!is.finite(slope) || slope <= 0) { step <- ev$gradient; slope <- sum(step^2) }
    tt <- 1; ok <- FALSE
    while (tt > 1e-12) {
      ev_new <- rate_mod_eval(par + tt * step, cr, B, hp, sigma2_a, nu, deriv = 2)
      if (is.finite(ev_new$value) && ev_new$value >= ev$value + 1e-4 * tt * slope) { ok <- TRUE; break }
      tt <- tt / 2
    }
    if (!ok) break
    par <- par + tt * step
    ev <- ev_new
    if (slope < tol || max(abs(ev$gradient)) < tol * (1 + abs(ev$value))) { converged <- TRUE; break }
  }
  hc <- hessian_to_cov(ev$hessian, rel_floor = 1e-10)
  logdet <- sum(log(pmax(hc$eigenvalues, hc$floor)))
  up <- .rate_unpack_par(par, cr, B)
  c(up, list(par = par, value = ev$value, loglik = ev$loglik, hessian = ev$hessian, vcov = hc$cov,
             log_evidence = ev$value + 0.5 * d * log(2 * pi) - 0.5 * logdet,
             sigma2_a = sigma2_a, converged = converged, iterations = it))
}

# sigma_a, and optionally sigma_kappa, by maximum marginal likelihood (Laplace
# approximation). The joint mode over (a, sigma_a) does not exist, so the
# standard deviations are not optimised together with the other parameters.
# sigma_a: NULL = estimate; a number = keep fixed.
rate_mod_fit_eb <- function(cr, B, hp, nu = Inf, interval = c(0.02, 10), start = NULL,
                            sigma_a = NULL, eb_kappa = FALSE, interval_kappa = c(0.01, 20)) {
  cache <- new.env()
  cache$start <- start
  K <- ncol(B)
  eb_kappa <- eb_kappa && K > 0
  fit_at <- function(ls_a, ls_k = NULL) {
    hp2 <- hp
    if (!is.null(ls_k)) hp2$sigma2_kappa <- exp(2 * ls_k)
    f <- rate_mod_fit(cr, B, hp2, sigma2_a = exp(2 * ls_a), nu = nu, start = cache$start)
    cache$start <- f$par
    f$sigma2_kappa <- hp2$sigma2_kappa
    f
  }
  la <- log(interval); lk <- log(interval_kappa)
  if (!eb_kappa) {
    if (!is.null(sigma_a)) return(fit_at(log(sigma_a)))
    op <- optimize(function(x) fit_at(x)$log_evidence, interval = la, maximum = TRUE, tol = 1e-4)
    fit <- fit_at(op$maximum)
    fit$sigma_a_at_boundary <- min(abs(op$maximum - la)) < 1e-2
  } else if (!is.null(sigma_a)) {
    op <- optimize(function(x) fit_at(log(sigma_a), x)$log_evidence, interval = lk, maximum = TRUE, tol = 1e-4)
    fit <- fit_at(log(sigma_a), op$maximum)
    fit$sigma_kappa_at_boundary <- min(abs(op$maximum - lk)) < 1e-2
  } else {
    # the surface can have a second local maximum at the lower bound of
    # sigma_kappa (no link), so the search is started from two values
    best <- NULL
    for (k0 in c(log(0.5), log(0.05))) {
      op <- optim(c(log(0.5), k0), function(th) -fit_at(th[1], th[2])$log_evidence, method = "L-BFGS-B",
                  lower = c(la[1], lk[1]), upper = c(la[2], lk[2]), control = list(factr = 1e9))
      if (is.null(best) || op$value < best$value) best <- op
    }
    fit <- fit_at(best$par[1], best$par[2])
    fit$sigma_a_at_boundary <- min(abs(best$par[1] - la)) < 1e-2
    fit$sigma_kappa_at_boundary <- min(abs(best$par[2] - lk)) < 1e-2
  }
  fit
}

#' Stage 2: the rate module
#'
#' Estimates the rate sub-process given the latent positions estimated in
#' stage 1: the coefficients \eqn{\xi} of the rate statistics, the activity
#' effects \eqn{a_i} and, with a link, the coefficients of the position-rate
#' link
#' \deqn{a_i \sim N(\mu_a + \textstyle\sum_l B_l(x_{il})^\top \kappa_l, \sigma_a^2),}{a_i ~ N(mu_a + sum_l B_l(x_il)' kappa_l, sigma_a^2),}
#' where \eqn{x_{il}}{x_il} are summaries of the positions of actor \eqn{i}
#' ([position_summaries()]) and \eqn{B_l} a natural cubic spline basis, with
#' columns centred, so that \eqn{\mu_a} is the average level, and scaled to
#' unit standard deviation across actors.
#'
#' The summaries are plugged in at the posterior means of the aligned
#' positions (plug-in posterior). With `n_cut > 0`, the cut posterior is
#' approximated as well: for each of `n_cut` draws of the positions the
#' summaries are recomputed, the basis is evaluated with the knots of the
#' plug-in fit, the mode of the rate module is found and one draw is taken from
#' its Laplace approximation.
#'
#' @inheritParams lsdynam_fit
#' @param sl statistics computed by [compute_stats()] (not needed when `cr`
#'   is given).
#' @param choice result of [lsdynam_choice_stage()] (or `fit$choice` of a
#'   fitted model).
#' @param n_draws_laplace number of draws from the Laplace approximation when
#'   `inference = "map"`.
#' @param exclude_from_link actors whose basis row is set to zero: they keep
#'   the exchangeable prior \eqn{N(\mu_a, \sigma_a^2)}{N(mu_a, sigma_a^2)} and
#'   do not inform the link (used by [link_leverage_check()]).
#' @param cr compressed rate statistics from a previous fit
#'   (`fit$rate$cr`), to avoid recomputing them.
#'
#' @return A list with the draws of `xi`, `a` (draws \eqn{\times n}{x n}),
#'   `mu_a` and `kappa_vec`, the posterior mode and Laplace approximation
#'   `map`, `sigma_a` and `sigma_kappa` with the method used to set them,
#'   the log marginal likelihood of the module (Laplace approximation)
#'   `log_evidence`, the acceptance rate, the link specification `link` (with
#'   the plug-in summaries `X`, the basis `B` and its specification), the
#'   cut-posterior draws `cut` (or `NULL`), the compressed statistics `cr`,
#'   the priors and the computing time.
#' @seealso [lsdynam_fit()], [link_comparison_table()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 300, seed = 2)
#' sl <- compute_stats(ev, n = 15, stats = c("rate_indeg", "choice_recip"))
#' ch <- lsdynam_choice_stage(ev, n = 15, D = 1, sl = sl, inference = "map")
#' rt <- lsdynam_rate_stage(ev, n = 15, sl = sl, choice = ch, link = "pop",
#'                          inference = "map")
#' rt$sigma_a
#' @export
lsdynam_rate_stage <- function(events_df, n, sl, choice,
                               link = if (choice$D > 0) c("pop", "rho") else "none",
                               link_type = c("spline", "linear"), df = 4, winsorize_q = 1,
                               priors = lsdynam_priors(), sigma_a = "eb", nu = Inf,
                               inference = c("mcmc", "map"), sampler = c("mala", "rwm"),
                               n_iter = 6000, burn_in = 1000, thin = 5, n_draws_laplace = 1000,
                               n_cut = 0, cut_refit_sigma = FALSE, rate_event = NULL,
                               exclude_from_link = NULL, plugin = c("positions", "summary"),
                               scale_basis = TRUE, sigma_kappa = NULL,
                               seed = 1, verbose = FALSE, cr = NULL) {
  inference <- match.arg(inference); sampler <- match.arg(sampler); link_type <- match.arg(link_type)
  plugin <- match.arg(plugin)
  t0 <- Sys.time()
  if (is.null(cr)) cr <- compress_rate_stats(sl, events_df$sender, n, rate_event = rate_event)
  hp <- priors
  has_link <- !identical(link, "none") && length(link) > 0
  if (has_link && choice$D == 0) stop("a link needs latent positions (D >= 1)")
  basis <- NULL; X <- NULL
  if (has_link) {
    ps <- position_summaries(choice$U_mean, choice$V_mean)
    if (!all(link %in% names(ps))) stop("unknown link variable(s): ", paste(setdiff(link, names(ps)), collapse = ", "))
        X <- as.matrix(ps[, link, drop = FALSE])
    if (plugin == "summary") {
      # posterior mean of the summary itself, E(x_i | events), instead of the
      # summary of the posterior-mean positions
      X <- apply(summaries_draws(choice$U, choice$V, link), c(2, 3), mean)
      X <- matrix(X, nrow = n, dimnames = list(NULL, link))
    }
    basis <- build_link_basis(X, type = link_type, df = df, winsorize_q = winsorize_q, scale = scale_basis)
    B <- basis$B
    if (length(exclude_from_link) > 0) B[exclude_from_link, ] <- 0
  } else {
    B <- matrix(0, n, 0)
  }
  K <- ncol(B)

    # sigma_kappa: NULL = prior variance priors$sigma2_kappa; a number = that
  # prior standard deviation; "eb" = maximum marginal likelihood
  eb_kappa <- identical(sigma_kappa, "eb") && K > 0
  if (is.numeric(sigma_kappa)) hp$sigma2_kappa <- sigma_kappa^2
  sigma_method <- if (identical(sigma_a, "eb")) "eb" else "fixed"
  if (sigma_method == "eb" || eb_kappa) {
    map <- rate_mod_fit_eb(cr, B, hp, nu = nu, sigma_a = if (sigma_method == "eb") NULL else sigma_a,
                           eb_kappa = eb_kappa)
    if (eb_kappa) hp$sigma2_kappa <- map$sigma2_kappa
  } else {
    map <- rate_mod_fit(cr, B, hp, sigma2_a = sigma_a^2, nu = nu)
  }
  sigma2_a <- map$sigma2_a
    if (verbose) cat(sprintf("  rate mode: sigma_a = %.3f (%s)%s, Newton iterations %d, converged %s\n",
                           sqrt(sigma2_a), sigma_method,
                           if (eb_kappa) sprintf(", sigma_kappa = %.3f (eb)", sqrt(hp$sigma2_kappa)) else "",
                           map$iterations, map$converged))
  d <- length(map$par)
  accept <- NA; eps <- NA
  if (inference == "map") {
    set.seed(seed + 2500)
    Lr <- t(chol(make_pd(map$vcov)))
    samples <- matrix(map$par, n_draws_laplace, d, byrow = TRUE) +
      t(Lr %*% matrix(rnorm(d * n_draws_laplace), d, n_draws_laplace))
  } else if (sampler == "mala") {
    fn <- function(p) rate_mod_eval(p, cr, B, hp, sigma2_a, nu, deriv = 1)
    mc <- mala_sample(fn, map$par, map$vcov, n_iter = n_iter, burn_in = burn_in, thin = thin, seed = seed + 2000)
    samples <- mc$samples; accept <- mc$accept_rate; eps <- mc$eps
  } else {
    log_post <- function(p) rate_mod_eval(p, cr, B, hp, sigma2_a, nu, deriv = 0)$value
    mc <- rwm_sample(log_post, map$par, map$vcov, n_iter = n_iter, burn_in = burn_in, thin = thin, seed = seed + 2000)
    samples <- mc$samples; accept <- mc$accept_rate
  }
  R <- cr$R
  xi_draws <- samples[, seq_len(R), drop = FALSE]; colnames(xi_draws) <- cr$rnames
  a_draws  <- samples[, R + seq_len(n), drop = FALSE]
  mu_draws <- samples[, R + n + 1]
  kap_draws <- samples[, R + n + 1 + seq_len(K), drop = FALSE]; colnames(kap_draws) <- colnames(B)

  cut <- NULL
  if (has_link && n_cut > 0) {
    set.seed(seed + 3000)
    ns1 <- dim(choice$U)[1]
    idx <- unique(round(seq(1, ns1, length.out = min(n_cut, ns1))))
    D <- choice$D
    cut_mu <- numeric(length(idx)); cut_kap <- matrix(NA_real_, length(idx), K, dimnames = list(NULL, colnames(B)))
    cut_xi <- matrix(NA_real_, length(idx), R, dimnames = list(NULL, cr$rnames))
    cut_sig <- numeric(length(idx)); cut_X <- array(NA_real_, c(length(idx), n, length(link)))
    for (k in seq_along(idx)) {
      ps_s <- position_summaries(matrix(choice$U[idx[k], , ], n, D), matrix(choice$V[idx[k], , ], n, D))
      X_s <- as.matrix(ps_s[, link, drop = FALSE])
      B_s <- eval_link_basis(basis$spec, X_s)
      if (length(exclude_from_link) > 0) B_s[exclude_from_link, ] <- 0
            f_s <- if (cut_refit_sigma && sigma_method == "eb") rate_mod_fit_eb(cr, B_s, hp, nu = nu, start = map$par) else
        rate_mod_fit(cr, B_s, hp, sigma2_a = sigma2_a, nu = nu, start = map$par)   # sigma_kappa stays at its plug-in value
      dr <- f_s$par + as.numeric(t(chol(make_pd(f_s$vcov))) %*% rnorm(d))
      cut_xi[k, ] <- dr[seq_len(R)]; cut_mu[k] <- dr[R + n + 1]; cut_kap[k, ] <- dr[R + n + 1 + seq_len(K)]
      cut_sig[k] <- sqrt(f_s$sigma2_a); cut_X[k, , ] <- X_s
    }
    cut <- list(mu_a = cut_mu, kappa_vec = cut_kap, xi = cut_xi, sigma_a = cut_sig, X = cut_X, idx = idx)
  }

  list(xi = xi_draws, a = a_draws, mu_a = mu_draws, kappa_vec = kap_draws,
       map = map, sigma_a = sqrt(sigma2_a), sigma_a_method = sigma_method, nu = nu,
       sigma_kappa = sqrt(hp$sigma2_kappa), sigma_kappa_method = if (eb_kappa) "eb" else "fixed",
       log_evidence = map$log_evidence, accept_rate = accept, eps = eps,
       inference = inference, sampler = if (inference == "mcmc") sampler else "laplace",
       link = list(vars = if (has_link) link else "none", type = link_type, df = df,
                   winsorize_q = winsorize_q, X = X, basis = basis, B = B,
                   exclude = exclude_from_link, plugin = plugin, scale_basis = scale_basis),
       cut = cut, cr = cr, priors = hp,
       seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
