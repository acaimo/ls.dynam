#' Predictive scores on held-out events
#'
#' Logarithmic scores of a model fitted to the first `m_train` events, for
#' the events that follow, with the statistics computed from the full
#' observed history. The log-likelihood separates into a rate part (senders
#' and times) and a choice part (receivers given senders), and the rate part
#' separates further by actor.
#'
#' The scores are sums of log posterior predictive densities over units: the
#' actors for the rate part and the held-out events for the choice part.
#' `rate_joint` and `choice_joint` are the logs of the posterior means of the
#' likelihood of all held-out events together; they are given for comparison
#' only, since with a finite number of draws this Monte Carlo estimate is
#' dominated by the largest likelihoods and biased downwards.
#'
#' @param fit an object of class `lsdynam` fitted to events `1:m_train` of
#'   `events_full`.
#' @param events_full all the events (training and held-out).
#' @param n number of actors.
#' @param m_train number of training events.
#' @param n_draws number of posterior draws used.
#' @param rate_event as in [lsdynam_fit()], over all the events.
#'
#' @return A list with the scores `rate` and `choice`, their unit-level terms
#'   `rate_by_actor` and `choice_by_event`, `rate_joint`, `choice_joint`, the
#'   number of held-out events `n_test`, the training out-degrees
#'   `N_out_train` and the senders of the held-out events.
#' @seealso [score_difference()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' m_tr <- 400
#' f0 <- lsdynam_fit(ev[1:m_tr, ], n = 15, D = 0, link = "none", inference = "map", verbose = FALSE)
#' f1 <- lsdynam_fit(ev[1:m_tr, ], n = 15, D = 1, link = "none", inference = "map", verbose = FALSE)
#' s0 <- predictive_scores(f0, ev, n = 15, m_train = m_tr)
#' s1 <- predictive_scores(f1, ev, n = 15, m_train = m_tr)
#' score_difference(s1, s0)
#' @export
predictive_scores <- function(fit, events_full, n, m_train, n_draws = 200, rate_event = NULL) {
  M <- nrow(events_full)
  test <- (m_train + 1):M
  sl <- compute_stats(events_full, n, stats = fit$stats)
  senders <- events_full$sender; receivers <- events_full$receiver

  # ---- rate part
  cr <- compress_rate_stats(sl, senders, n, rows = test, rate_event = rate_event)
  nr <- nrow(fit$rate$a)
  ir <- unique(round(seq(1, nr, length.out = min(n_draws, nr))))
  # sum of the rate statistics over the held-out events sent by each actor
  # (the event term xi's(i_m;t_m) is attributed to the sending actor)
  ev_ok <- if (is.null(rate_event)) rep(TRUE, length(test)) else rate_event[test]
  S_ev_actor <- matrix(0, n, cr$R)
  for (r in seq_len(cr$R)) {
    st <- sl[[cr$rnames[r]]][cbind(test, senders[test])][ev_ok]
    S_ev_actor[, r] <- as.numeric(tapply(st, factor(senders[test][ev_ok], levels = seq_len(n)), sum, default = 0))
  }
  ll_actor <- matrix(NA_real_, length(ir), n)
  for (k in seq_along(ir)) {
    xi <- fit$rate$xi[ir[k], ]; a <- fit$rate$a[ir[k], ]
    lin <- a[cr$actor]
    if (cr$R > 0) lin <- lin + as.numeric(cr$S %*% xi)
    Lam_actor <- as.numeric(rowsum(cr$expo * exp(lin), cr$actor, reorder = TRUE))
    ev_term <- cr$N * a
    if (cr$R > 0) ev_term <- ev_term + as.numeric(S_ev_actor %*% xi)
    ll_actor[k, ] <- ev_term - Lam_actor
  }
  rate_joint <- log_mean_exp(rowSums(ll_actor))
  rate_actor <- apply(ll_actor, 2, log_mean_exp)

  # ---- choice part
  cnames <- choice_stat_names(sl)
  sl_test <- lapply(sl[cnames], function(m) m[test, , drop = FALSE])
  nc <- nrow(fit$choice$zeta)
  ic <- unique(round(seq(1, nc, length.out = min(n_draws, nc))))
  D <- fit$D
  ll_event <- matrix(NA_real_, length(ic), length(test))
  idx <- cbind(seq_along(test), receivers[test])
  for (k in seq_along(ic)) {
    U <- matrix(fit$choice$U[ic[k], , ], n, D); V <- matrix(fit$choice$V[ic[k], , ], n, D)
    ep <- choice_eta_pi(fit$choice$zeta[ic[k], ], U, V, sl_test, senders[test], cnames)
    ll_event[k, ] <- ep$eta[idx] - ep$log_psi
  }
  choice_joint <- log_mean_exp(rowSums(ll_event))
  choice_event <- apply(ll_event, 2, log_mean_exp)

  list(rate = sum(rate_actor), choice = sum(choice_event),
       rate_joint = rate_joint, choice_joint = choice_joint,
       rate_by_actor = rate_actor, choice_by_event = choice_event,
       n_test = length(test), m_train = m_train, n_draws = c(rate = length(ir), choice = length(ic)),
       N_out_train = tabulate(senders[seq_len(m_train)], nbins = n),
       test_sender = senders[test])
}

#' Differences of predictive scores
#'
#' Differences between the predictive scores of two models for the rate part
#' (sum over actors) and the choice part (sum over held-out events), with the
#' standard error of each sum of unit-level differences,
#' \eqn{\sqrt{k\,\mathrm{var}(d)}}{sqrt(k var(d))} for \eqn{k} units, as for
#' leave-one-out comparisons (Vehtari, Gelman and Gabry, 2017).
#'
#' @param sa,sb results of [predictive_scores()] for two models; positive
#'   differences favour `sa`.
#' @param groups optional factor over the actors (for example classes of
#'   training out-degree) that splits the rate part.
#'
#' @return A data frame with the part, the number of units, the difference
#'   and its standard error.
#' @references Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical
#'   Bayesian model evaluation using leave-one-out cross-validation and WAIC.
#'   *Statistics and Computing*, 27, 1413--1432.
#' @inherit predictive_scores examples
#' @export
score_difference <- function(sa, sb, groups = NULL) {
  se <- function(d) if (length(d) > 1) sqrt(length(d) * stats::var(d)) else NA_real_
  dr <- sa$rate_by_actor - sb$rate_by_actor
  dc <- sa$choice_by_event - sb$choice_by_event
  out <- data.frame(part = c("rate", "choice"), units = c(length(dr), length(dc)),
                    difference = c(sum(dr), sum(dc)), se = c(se(dr), se(dc)),
                    stringsAsFactors = FALSE)
  if (!is.null(groups)) {
    g <- do.call(rbind, lapply(levels(groups), function(l) {
      d <- dr[groups == l]
      data.frame(part = paste0("rate, ", l), units = length(d), difference = sum(d), se = se(d),
                 stringsAsFactors = FALSE)
    }))
    out <- rbind(out, g)
  }
  out
}
