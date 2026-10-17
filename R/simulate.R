#' Parameters of a latent space DyNAM with a known link
#'
#' Draws latent positions and activity effects for a simulation in which the
#' position-rate link is known. The sender positions are drawn from
#' \eqn{N(\text{mean\_u}, \text{sd\_u}^2 I)}{N(mean_u, sd_u^2 I)} and the
#' receiver positions from \eqn{N(0, \text{sd\_v}^2 I)}{N(0, sd_v^2 I)}; the
#' summaries `link_var` of the positions are computed with
#' [position_summaries()], and the activity effects are
#' \eqn{a_i = \mu_a + \sum_l g_l(x_{il}) + \varepsilon_i}{a_i = mu_a + sum_l g_l(x_il) + e_i}
#' with \eqn{\varepsilon_i \sim N(0, \sigma_a^2)}{e_i ~ N(0, sigma_a^2)}.
#' Each component \eqn{g_l} is a function of the centred summary, and its
#' average over the actors is removed, so that \eqn{\mu_a} is the average
#' level, as in the fitted model.
#'
#' With `mean_u` different from zero the senders share a direction of
#' preference and latent popularity has a variance comparable to that of
#' excess self-affinity; with `mean_u = 0` the mean sender position is close
#' to zero, and so is latent popularity.
#'
#' @param n number of actors.
#' @param D dimension of the latent space.
#' @param link_fun a function of the centred link variable or, for an
#'   additive link, a list with one function per element of `link_var`
#'   (`NULL` for a component without effect). `NULL`: no link.
#' @param link_var names of the summaries in the link, columns of
#'   [position_summaries()].
#' @param mu_a average activity effect.
#' @param sigma_a residual standard deviation of the activity effects.
#' @param xi named vector of coefficients of the rate statistics (names as in
#'   [compute_stats()]).
#' @param zeta named vector of coefficients of the choice statistics.
#' @param sd_u,sd_v standard deviations of the sender and receiver positions.
#' @param mean_u mean of the sender positions (a number or a vector of
#'   length `D`).
#' @param seed seed of the random number generator (`NULL`: use the current
#'   state). The function calls [set.seed()].
#'
#' @return A list with `theta` (a list with `xi`, `zeta`, `a`, `U`, `V`, the
#'   input of [simulate_lsdynam()]), the summaries of the positions
#'   `summaries`, the link variables `x`, the true components (`g`,
#'   `component_fun`), `level_fun(x_raw, component)`, the level of one
#'   component as a function of the raw summary, and the names of the
#'   statistics `stats`.
#' @seealso [simulate_lsdynam()].
#' @examples
#' f <- function(x) 1.5 * (plogis(3 * x) - 0.5)
#' tr <- make_truth(n = 20, D = 1, link_fun = list(f, NULL), link_var = c("pop", "rho"),
#'                  mean_u = 1, seed = 1)
#' str(tr$theta)
#' @export
make_truth <- function(n, D, link_fun = NULL, link_var = "r", mu_a = 0, sigma_a = 0.5,
                       xi = c(rate_indeg = 0.1), zeta = c(choice_recip = 0.5),
                       sd_u = 1, sd_v = 1, mean_u = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  U <- matrix(rnorm(n * D, 0, sd_u), n, D)
  V <- matrix(rnorm(n * D, 0, sd_v), n, D)
  if (any(mean_u != 0)) U <- sweep(U, 2, rep_len(mean_u, D), "+")
  ps <- position_summaries(U, V)
  L <- length(link_var)
  funs <- if (is.list(link_fun)) link_fun else rep(list(link_fun), L)
  if (length(funs) != L) stop("make_truth: give one link function per link variable")
  X <- as.matrix(ps[, link_var, drop = FALSE])
  centers <- colMeans(X)
  comp_raw <- lapply(seq_len(L), function(l)
    if (is.null(funs[[l]])) rep(0, n) else funs[[l]](X[, l] - centers[l]))
  g_means <- vapply(comp_raw, mean, numeric(1))
  g <- Reduce(`+`, lapply(seq_len(L), function(l) comp_raw[[l]] - g_means[l]))
  a <- mu_a + g + rnorm(n, 0, sigma_a)
  comp_fun <- lapply(seq_len(L), function(l) {
    f <- funs[[l]]; cl <- centers[l]; ml <- g_means[l]
    function(x_raw) (if (is.null(f)) rep(0, length(x_raw)) else f(x_raw - cl)) - ml
  })
  names(comp_fun) <- link_var
  # level of one component as a function of its raw summary (the other
  # components at their average contribution, which is zero)
  level_fun <- function(x_raw, component = 1) mu_a + comp_fun[[component]](x_raw)
  list(theta = list(xi = xi, zeta = zeta, a = a, U = U, V = V),
       n = n, D = D, summaries = ps, link_var = link_var,
       x = if (L == 1) X[, 1] else X, center = if (L == 1) unname(centers) else centers,
       level_fun = level_fun, component_fun = comp_fun, g = g, mu_a = mu_a, sigma_a = sigma_a,
       stats = c(names(xi), names(zeta)))
}

.rate_feats_now <- function(rate_names, outdeg, indeg) {
  n <- length(outdeg)
  if (length(rate_names) == 0) return(matrix(0, n, 0))
  feat <- matrix(0, n, length(rate_names), dimnames = list(NULL, rate_names))
  for (rn in rate_names) {
    feat[, rn] <- switch(rn, rate_outdeg = log1p(outdeg), rate_indeg = log1p(indeg),
                         stop("unrecognized rate stat: ", rn))
  }
  feat
}

.choice_feats_now <- function(choice_names, i, outdeg, indeg, inertia, recip) {
  n <- length(outdeg)
  if (length(choice_names) == 0) return(matrix(0, n, 0))
  feat <- matrix(0, n, length(choice_names), dimnames = list(NULL, choice_names))
  for (cn in choice_names) {
    feat[, cn] <- switch(cn,
      choice_outdeg = log1p(outdeg), choice_indeg = log1p(indeg),
      choice_inertia = log1p(inertia[i, ]), choice_recip = log1p(recip[, i]),
      stop("unrecognized choice stat: ", cn))
  }
  feat
}

#' Simulate relational events from a latent space DyNAM
#'
#' Generates events one at a time: the waiting time is exponential with rate
#' \eqn{\sum_i \lambda_i(t)}{sum_i lambda_i(t)}, the sender is drawn with
#' probabilities proportional to \eqn{\lambda_i(t)}{lambda_i(t)} and the
#' receiver from the choice model, with the statistics updated after every
#' event exactly as in [compute_stats()].
#'
#' @param theta list with the coefficients `xi` and `zeta` (named vectors,
#'   names as in [compute_stats()]), the activity effects `a` and the
#'   \eqn{n \times D}{n x D} position matrices `U` and `V`, for example
#'   `make_truth(...)$theta`.
#' @param n number of actors.
#' @param stats names of the statistics (default `names(theta$xi)` and
#'   `names(theta$zeta)`).
#' @param M number of events to simulate.
#' @param T_max alternatively, the end of the observation period.
#' @param seed seed of the random number generator (`NULL`: use the current
#'   state). The function calls [set.seed()].
#' @param max_events upper limit on the number of events when `T_max` is
#'   given.
#'
#' @return A data frame with columns `time`, `sender` and `receiver`.
#' @examples
#' tr <- make_truth(n = 15, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 200, seed = 2)
#' head(ev)
#' @export
simulate_lsdynam <- function(theta, n, stats = NULL, M = NULL, T_max = NULL, seed = NULL,
                             max_events = 500000) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(M) && is.null(T_max)) stop("give M or T_max")
  if (is.null(stats)) stats <- c(names(theta$xi), names(theta$zeta))
  rate_names   <- grep("^rate_",   stats, value = TRUE)
  choice_names <- grep("^choice_", stats, value = TRUE)
  cap <- if (!is.null(M)) M else 10000
  ev_time <- numeric(cap); ev_send <- integer(cap); ev_recv <- integer(cap)
  outdeg <- integer(n); indeg <- integer(n)
  inertia <- matrix(0, n, n); recip <- matrix(0, n, n)
  D <- ncol(theta$U)
  t <- 0; m <- 0
  repeat {
    rate_feat <- .rate_feats_now(rate_names, outdeg, indeg)
    lin_rate  <- if (length(rate_names) > 0) as.numeric(rate_feat %*% theta$xi[rate_names]) else 0
    lambda_i  <- exp(lin_rate + theta$a)
    tot <- sum(lambda_i)
    if (!is.finite(tot) || tot <= 0) { warning("simulate_lsdynam: non-finite total rate"); break }
    t <- t + rexp(1, rate = tot)
    if (!is.null(T_max) && t > T_max) break
    sender <- sample.int(n, 1, prob = lambda_i)

    choice_feat <- .choice_feats_now(choice_names, sender, outdeg, indeg, inertia, recip)
    eta <- if (length(choice_names) > 0) as.numeric(choice_feat %*% theta$zeta[choice_names]) else rep(0, n)
    if (D > 0) eta <- eta + as.numeric(theta$V %*% theta$U[sender, ])
    eta[sender] <- -Inf
    p <- exp(eta - max(eta))
    receiver <- sample.int(n, 1, prob = p)

    m <- m + 1
    if (m > cap) {
      cap <- 2 * cap
      length(ev_time) <- cap; length(ev_send) <- cap; length(ev_recv) <- cap
    }
    ev_time[m] <- t; ev_send[m] <- sender; ev_recv[m] <- receiver
    outdeg[sender] <- outdeg[sender] + 1L; indeg[receiver] <- indeg[receiver] + 1L
    inertia[sender, receiver] <- inertia[sender, receiver] + 1
    recip[sender, receiver]   <- recip[sender, receiver]   + 1
    if (!is.null(M) && m >= M) break
    if (m >= max_events) break
  }
  data.frame(time = ev_time[seq_len(m)], sender = ev_send[seq_len(m)], receiver = ev_recv[seq_len(m)])
}
