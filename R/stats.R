#' Statistics of the rate and choice sub-processes
#'
#' Computes, for every event \eqn{m} and every actor, the statistics of the
#' rate sub-process \eqn{s(i;t_m)} and of the choice sub-process
#' \eqn{s(i_m,j;t_m)} from the history of the process before the event.
#' Counts start at zero at the first event of `events_df`.
#'
#' The available statistics, all of the form \eqn{\log(1 + \text{count})}{log(1 + count)}, are
#' * `rate_outdeg`, `rate_indeg`: number of events sent and received by
#'   actor \eqn{i} (rate sub-process);
#' * `choice_outdeg`, `choice_indeg`: number of events sent and received by
#'   the potential receiver \eqn{j};
#' * `choice_inertia`: number of events sent by the sender \eqn{i_m} to
#'   \eqn{j};
#' * `choice_recip`: number of events sent by \eqn{j} to the sender
#'   \eqn{i_m} (reciprocity).
#'
#' @param events_df data frame of events ordered by time, with columns
#'   `sender`, `receiver` (labels \eqn{1, \dots, n}) and `time`.
#' @param n number of actors.
#' @param stats names of the statistics to compute.
#'
#' @return A list with an \eqn{M \times n}{M x n} matrix for each statistic
#'   (row \eqn{m}: values before event \eqn{m}; for a choice statistic, column
#'   \eqn{j} is the value for receiver \eqn{j} given the sender of event
#'   \eqn{m}), the waiting times `delta_t` (with \eqn{t_0 = 0}), the number of
#'   events `M`, and the numbers of rate and choice statistics `R` and `S`.
#' @examples
#' ev <- data.frame(sender = c(1, 2, 1), receiver = c(2, 1, 2), time = c(0.5, 1, 2))
#' sl <- compute_stats(ev, n = 2, stats = c("rate_indeg", "choice_recip"))
#' sl$choice_recip
#' @export
compute_stats <- function(events_df, n,
                           stats = c("rate_outdeg", "rate_indeg",
                                     "choice_outdeg", "choice_indeg",
                                     "choice_inertia", "choice_recip")) {
  M   <- nrow(events_df)
  has <- function(x) x %in% stats
  stats <- unique(stats)

  if (has("rate_outdeg"))    rate_outdeg    <- matrix(0, M, n)
  if (has("rate_indeg"))     rate_indeg     <- matrix(0, M, n)
  if (has("choice_outdeg"))  choice_outdeg  <- matrix(0, M, n)
  if (has("choice_indeg"))   choice_indeg   <- matrix(0, M, n)
  if (has("choice_inertia")) choice_inertia <- matrix(0, M, n)
  if (has("choice_recip"))   choice_recip   <- matrix(0, M, n)

  outdeg <- integer(n)
  indeg  <- integer(n)
  if (has("choice_inertia")) inertia <- matrix(0, n, n)
  if (has("choice_recip"))   recip   <- matrix(0, n, n)

  for (m in seq_len(M)) {
    i <- events_df$sender[m]
    j <- events_df$receiver[m]

    if (has("rate_outdeg"))    rate_outdeg[m, ]    <- log1p(outdeg)
    if (has("rate_indeg"))     rate_indeg[m, ]     <- log1p(indeg)
    if (has("choice_outdeg"))  choice_outdeg[m, ]  <- log1p(outdeg)
    if (has("choice_indeg"))   choice_indeg[m, ]   <- log1p(indeg)
    if (has("choice_inertia")) choice_inertia[m, ] <- log1p(inertia[i, ])
    if (has("choice_recip"))   choice_recip[m, ]   <- log1p(recip[, i])

    outdeg[i] <- outdeg[i] + 1
    indeg[j]  <- indeg[j]  + 1
    if (has("choice_inertia")) inertia[i, j] <- inertia[i, j] + 1
    if (has("choice_recip"))   recip[i, j]   <- recip[i, j]   + 1
  }

  times   <- c(0, events_df$time)
  delta_t <- diff(times)

  out <- list(M = M, delta_t = delta_t)
  if (has("rate_outdeg"))    out$rate_outdeg    <- rate_outdeg
  if (has("rate_indeg"))     out$rate_indeg     <- rate_indeg
  if (has("choice_outdeg"))  out$choice_outdeg  <- choice_outdeg
  if (has("choice_indeg"))   out$choice_indeg   <- choice_indeg
  if (has("choice_inertia")) out$choice_inertia <- choice_inertia
  if (has("choice_recip"))   out$choice_recip   <- choice_recip
  out$R <- sum(grepl("rate_",   stats))
  out$S <- sum(grepl("choice_", stats))
  out
}

#' @noRd
rate_stat_names   <- function(sl) grep("^rate_",   names(sl), value = TRUE)

#' @noRd
choice_stat_names <- function(sl) grep("^choice_", names(sl), value = TRUE)
