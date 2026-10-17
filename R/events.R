#' Read time-stamped events from a KONECT file
#'
#' Reads a relational event file in the format of the KONECT network
#' collection: one event per line with sender, receiver, weight and time
#' stamp, and comment lines starting with `%`. This is the format of the
#' Radoslaw e-mail data (file `out.radoslaw_email_email`, available from
#' <https://konect.cc/networks/radoslaw_email>).
#'
#' @param file path of the file.
#' @param time_col column that holds the time stamp (4 in KONECT files).
#'
#' @return A data frame with integer columns `sender` and `receiver` and a
#'   numeric column `time`, ordered by time (events with the same time keep
#'   the order of the file).
#'
#' @seealso [prepare_events()] for the next step.
#' @examples
#' f <- tempfile()
#' writeLines(c("% sym unweighted", "1 2 1 100", "2 1 1 160", "1 3 1 130"), f)
#' read_konect_events(f)
#' @export
read_konect_events <- function(file, time_col = 4) {
  raw <- utils::read.table(file, header = FALSE, comment.char = "%")
  if (ncol(raw) < max(3, time_col)) stop("expected columns sender, receiver and a time stamp in column ", time_col)
  ev <- data.frame(sender = as.integer(raw[[1]]), receiver = as.integer(raw[[2]]),
                   time = as.numeric(raw[[time_col]]))
  ev <- ev[order(ev$time), , drop = FALSE]
  rownames(ev) <- NULL
  ev
}

#' Prepare relational events for model fitting
#'
#' Orders the events by time, removes self-loops (an actor is not in its own
#' choice set), identifies messages to several recipients, rescales the time
#' stamps and separates exact ties, as in the analysis of the Radoslaw e-mail
#' data.
#'
#' A message to \eqn{k} recipients is usually stored as \eqn{k} events with
#' the same sender and time stamp. Each of them is a choice of receiver; the
#' column `msg_id` identifies the message, so that only its first event can be
#' counted as a sending decision in the rate module, with
#' `rate_event = !duplicated(events$msg_id)` in [lsdynam_fit()] (recompute it
#' after any filtering of the events).
#'
#' @param events data frame with columns `sender`, `receiver` and `time`.
#' @param time_scale the time stamps are divided by this number after
#'   subtracting the first one; the default turns seconds into days.
#' @param break_ties if `TRUE`, events with exactly the same time are separated
#'   by adding an increasing sequence of at most `tie_eps` to the times, which
#'   keeps their order. The rate likelihood needs positive waiting times.
#' @param tie_eps total amount added to separate ties.
#' @param verbose print a short report.
#'
#' @return A data frame with columns `sender`, `receiver`, `time` and
#'   `msg_id`, ordered by time.
#' @seealso [select_time_window()], [filter_active_both()].
#' @examples
#' ev <- data.frame(sender = c(1, 1, 2, 3, 3), receiver = c(2, 3, 1, 3, 1),
#'                  time = c(0, 0, 3600, 7200, 9000))
#' prepare_events(ev)
#' @export
prepare_events <- function(events, time_scale = 86400, break_ties = TRUE, tie_eps = 1e-6,
                           verbose = TRUE) {
  if (!all(c("sender", "receiver", "time") %in% names(events)))
    stop("events needs columns sender, receiver and time")
  ev <- data.frame(sender = as.integer(events$sender), receiver = as.integer(events$receiver),
                   time = as.numeric(events$time))
  ev <- ev[order(ev$time), , drop = FALSE]
  n_self <- sum(ev$sender == ev$receiver)
  if (n_self > 0) ev <- ev[ev$sender != ev$receiver, , drop = FALSE]
  key <- paste(ev$sender, ev$time)
  ev$msg_id <- match(key, unique(key))
  ev$time <- (ev$time - min(ev$time)) / time_scale
  n_ties <- sum(duplicated(ev$time))
  if (break_ties && n_ties > 0) ev$time <- ev$time + seq(0, tie_eps, length.out = nrow(ev))
  rownames(ev) <- NULL
  if (verbose) {
    cat(sprintf("Events: %d (self-loops removed: %d); messages: %d; exact time ties: %d%s\n",
                nrow(ev), n_self, length(unique(ev$msg_id)), n_ties,
                if (break_ties && n_ties > 0) " (separated)" else ""))
  }
  ev
}

#' Select a time window
#'
#' Keeps the events within a time window, given as fractions of the
#' observation period or as times, and re-anchors the times at zero. The
#' statistics computed by [compute_stats()] then start from zero at the
#' beginning of the window.
#'
#' @param events_df data frame of events with a column `time`.
#' @param start_frac,end_frac start and end of the window as fractions of the
#'   observed time range (for example 0.2 and 0.4 for the second fifth).
#' @param start_time,end_time start and end of the window in the units of
#'   `events_df$time`, used when the fractions are not given.
#' @param verbose print a short report.
#'
#' @return The events in the window, with times measured from the first of
#'   them.
#' @examples
#' ev <- data.frame(sender = c(1, 2, 1, 3), receiver = c(2, 1, 3, 1), time = 1:4)
#' select_time_window(ev, start_frac = 0.25, end_frac = 1)
#' @export
select_time_window <- function(events_df, start_frac = NULL, end_frac = NULL,
                               start_time = NULL, end_time = NULL, verbose = TRUE) {
  full_range <- range(events_df$time)
  if (!is.null(start_frac)) start_time <- full_range[1] + start_frac * diff(full_range)
  if (!is.null(end_frac))   end_time   <- full_range[1] + end_frac   * diff(full_range)
  if (is.null(start_time)) start_time <- full_range[1]
  if (is.null(end_time))   end_time   <- full_range[2]
  sub <- events_df[events_df$time >= start_time & events_df$time <= end_time, , drop = FALSE]
  if (nrow(sub) == 0) stop("select_time_window: no events in the window")
  sub$time <- sub$time - min(sub$time)
  rownames(sub) <- NULL
  if (verbose) {
    cat(sprintf("Full period: [%.1f, %.1f] (%d events); window: [%.1f, %.1f] (%d events, %.0f%%)\n",
                full_range[1], full_range[2], nrow(events_df), start_time, end_time, nrow(sub),
                100 * nrow(sub) / nrow(events_df)))
  }
  sub
}

#' Exclude actors
#'
#' Removes the events that involve the given actors and relabels the others
#' as \eqn{1, \dots, n}.
#'
#' @param events_df data frame of events with columns `sender` and `receiver`
#'   holding labels \eqn{1, \dots, n}.
#' @param n number of actors.
#' @param actors_to_exclude labels of the actors to remove.
#' @param verbose print a short report.
#'
#' @return A list with the relabelled events `events_df`, the number of
#'   actors `n` and `old_to_new`, a vector whose names are the old labels of
#'   the actors kept and whose values are their new labels.
#' @examples
#' ev <- data.frame(sender = c(1, 2, 3, 4), receiver = c(2, 3, 4, 1), time = 1:4)
#' exclude_actors(ev, n = 4, actors_to_exclude = 2)
#' @export
exclude_actors <- function(events_df, n, actors_to_exclude, verbose = TRUE) {
  keep <- setdiff(seq_len(n), actors_to_exclude)
  remap <- setNames(seq_along(keep), keep)
  out <- events_df[!(events_df$sender %in% actors_to_exclude) &
                     !(events_df$receiver %in% actors_to_exclude), , drop = FALSE]
  out$sender   <- unname(remap[as.character(out$sender)])
  out$receiver <- unname(remap[as.character(out$receiver)])
  rownames(out) <- NULL
  if (verbose) {
    cat(sprintf("Excluded actor(s): %s; events %d -> %d; actors %d -> %d\n",
                paste(actors_to_exclude, collapse = ", "), nrow(events_df), nrow(out), n, length(keep)))
  }
  list(events_df = out, n = length(keep), old_to_new = remap)
}

#' Keep the actors that both send and receive
#'
#' An actor that only sends, or only receives, carries no information about
#' one of its two positions. This function keeps the actors that send and
#' receive at least one event. Removing an actor removes its events, which can
#' leave another actor active in one role only, so the filter is repeated
#' until nothing changes. The actors are relabelled \eqn{1, \dots, n}.
#'
#' @param events_df data frame of events with columns `sender` and
#'   `receiver`.
#' @param verbose print a short report.
#'
#' @return A list with the filtered and relabelled events `events_df`, the
#'   number of actors `n`, `old_id`, the original label of each new label,
#'   and `n_pass`, the number of passes.
#' @examples
#' ev <- data.frame(sender = c(1, 2, 3, 1), receiver = c(2, 1, 1, 4), time = 1:4)
#' filter_active_both(ev)
#' @export
filter_active_both <- function(events_df, verbose = TRUE) {
  n_pass <- 0
  repeat {
    n_pass <- n_pass + 1
    keep <- intersect(unique(events_df$sender), unique(events_df$receiver))
    ok <- events_df$sender %in% keep & events_df$receiver %in% keep
    if (all(ok)) break
    events_df <- events_df[ok, , drop = FALSE]
    if (nrow(events_df) == 0) stop("filter_active_both: no events left")
  }
  old_id <- sort(keep)
  events_df$sender   <- match(events_df$sender, old_id)
  events_df$receiver <- match(events_df$receiver, old_id)
  rownames(events_df) <- NULL
  if (verbose) cat(sprintf("Active in both roles: %d actors, %d events (%d pass%s)\n",
                           length(old_id), nrow(events_df), n_pass, if (n_pass > 1) "es" else ""))
  list(events_df = events_df, n = length(old_id), old_id = old_id, n_pass = n_pass)
}
