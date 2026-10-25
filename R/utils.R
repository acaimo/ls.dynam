log_mean_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(mean(exp(x - m)))
}

#' Effective sample size of a Markov chain
#'
#' Effective sample size of one chain, with the autocorrelations summed by
#' Geyer's initial positive sequence estimator.
#'
#' @param x numeric vector of draws.
#' @return The effective sample size (`NA` for fewer than 10 draws or constant
#'   draws).
#' @seealso [mcmc_diagnostics()].
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = 0.8), n = 2000))
#' ess(x)
#' @export
ess <- function(x) {
  x <- as.numeric(x); n <- length(x)
  if (n < 10 || !is.finite(stats::var(x)) || stats::var(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(n - 2, 2000), plot = FALSE)$acf[-1]
  s <- 0; k <- 1
  while (k < length(ac)) {
    pr <- ac[k] + ac[k + 1]
    if (pr < 0) break
    s <- s + pr; k <- k + 2
  }
  min(n, n / (1 + 2 * s))
}
