# =============================================================================
# radoslaw_application.R
#
# Application of the latent space DyNAM with position-rate link to the
# Radoslaw e-mail network (Section 8 of the paper).
#
# The central specification is the ADDITIVE link in latent popularity and
# excess self-affinity, LINK = c("pop", "rho"), with r = pop + rho. The link
# in the self-affinity alone, r, is fitted on the same positions for
# comparison: it constrains the two components to act through their sum.
#
#   1. Load and clean
#   2. Time window and actors active in both roles. All actors are kept by
#      default; the runs of 8-9 October excluded actor 6
#      (EXCLUDE_ACTORS <- 6 before sourcing reproduces them).
#   3. Fits: M2 (additive link), modes of the choice module, link in r on
#      the same positions, M1 (no link), M0 (no latent positions)
#   4. Posterior summaries and convergence
#   5. Figures: positions, activity effects, components of the link, link in r
#   6. Evidence for the link, and for each component given the other
#   7. Links in different summaries on the same positions
#   8. Diagnostics: shrinkage, leverage, plug-in agreement
#   9. Sensitivity of the additive link
#  10. Held-out predictive comparison of M0, M1, M2
#  11. Calibration by simulation: no association, and popularity only
#  12. Numbers for the paper (printed at the end and written to
#      radoslaw_paper_numbers.txt)
#
# Running time on the earlier sample (119 actors, 16,086 events): M2 about 16
# minutes, the held-out fits about 10, each calibration about 60 (20
# replicates). RUN_NULL <- FALSE and RUN_NULL_POP <- FALSE before sourcing
# skip the calibrations; QUICK <- TRUE runs everything with short chains, to
# check that the script runs. Results for D > 1 (D = 2 in detail, with its
# own calibration with popularity only): radoslaw_dimension.R, which reads the
# .rds file written here.
#
# Requires the package ls.dynam and, in the working directory, the data file
# out_radoslaw.txt (the file out.radoslaw_email_email of the KONECT
# collection, https://konect.cc/networks/radoslaw_email).
# =============================================================================

library(ls.dynam)

# -----------------------------------------------------------------------------
# Settings (single place to change the analysis)
# -----------------------------------------------------------------------------

if (!exists("DATA_FILE"))      DATA_FILE      <- "out_radoslaw.txt"
if (!exists("QUICK"))          QUICK          <- FALSE   # TRUE: short chains, to check that the script runs
if (!exists("EXCLUDE_ACTORS")) EXCLUDE_ACTORS <- NULL    # original labels of actors to drop; NULL = none
if (!exists("RUN_HOLDOUT"))    RUN_HOLDOUT    <- TRUE
if (!exists("RUN_NULL"))       RUN_NULL       <- TRUE    # calibration with no association (slow)
if (!exists("RUN_NULL_POP"))   RUN_NULL_POP   <- TRUE    # calibration with popularity only (slow)
if (!exists("N_NULL"))         N_NULL         <- if (QUICK) 2 else 20

START_FRAC <- 0.2; END_FRAC <- 0.4      # the second fifth of the observation period

D          <- 1
STATS      <- c("rate_indeg", "choice_recip")
LINK       <- c("pop", "rho")   # central specification: latent popularity + excess self-affinity
LINK_R     <- "r"               # self-affinity alone, for comparison
DF         <- 4                 # spline basis functions per component
SIGMA_A    <- "eb"              # maximum marginal likelihood; or a number
WINSORIZE_Q <- 1                # < 1 caps both tails of each link variable at that quantile
# The basis functions are centred and scaled to unit standard deviation
# across actors (scale_basis = TRUE, the default), so sigma2_kappa = 1 is the
# prior variance of the contribution of one basis function to a typical
# actor's activity effect, whatever the summary.

# One e-mail to k recipients is stored as k events with the same time stamp.
#   "event"  : every row is a sending decision (as in the earlier scripts)
#   "message": every row is a choice of receiver, but only the first row of
#              each message counts as a sending decision in the rate model
MULTI_RECIPIENT <- "event"

SAMPLER  <- "mala"
N_ITER   <- if (QUICK) 1500 else 12000
BURN_IN  <- if (QUICK) 500  else 2000
THIN     <- if (QUICK) 2    else 5
if (!exists("N_STARTS")) N_STARTS <- if (QUICK) 2 else 10   # random starts of the choice module (plus a spectral start)
N_CUT    <- if (QUICK) 60   else 300    # stage-1 draws used for the cut posterior
TRAIN_FRAC <- 0.8
OUT_PREFIX <- "radoslaw"
SK_GRID    <- c(0.1, 0.25, 0.5, 1, 2)   # prior standard deviations of kappa for the Bayes factors

# =============================================================================
# 1. Load and clean
# =============================================================================

raw <- read.table(DATA_FILE, header = FALSE, comment.char = "%",
                  col.names = c("sender", "receiver", "weight", "timestamp"))
cat("Rows read:", nrow(raw), "\n")

events_df <- data.frame(sender = raw$sender, receiver = raw$receiver, time = as.numeric(raw$timestamp))
events_df <- events_df[order(events_df$time), ]

# self-loops are not valid choices (the sender is excluded from its own choice set)
n_selfloop <- sum(events_df$sender == events_df$receiver)
if (n_selfloop > 0) {
  cat("Removing", n_selfloop, "self-loop event(s)\n")
  events_df <- events_df[events_df$sender != events_df$receiver, ]
}

# rows with the same sender and the same time stamp belong to one message
msg_key <- paste(events_df$sender, events_df$time)
events_df$msg_id <- match(msg_key, unique(msg_key))
n_msg <- length(unique(events_df$msg_id))
cat(sprintf("Messages: %d, recipient-events: %d (%.1f%% of the rows are additional recipients)\n",
            n_msg, nrow(events_df), 100 * (1 - n_msg / nrow(events_df))))

# seconds -> days; separate exact ties by a negligible amount, keeping the order
events_df$time <- (events_df$time - min(events_df$time)) / 86400
n_ties <- sum(duplicated(events_df$time))
cat("Exact time ties:", n_ties, sprintf("(%.1f%%)\n", 100 * n_ties / nrow(events_df)))
if (n_ties > 0) events_df$time <- events_df$time + seq(0, 1e-6, length.out = nrow(events_df))

# =============================================================================
# 2. Time window and actors active in both roles
# =============================================================================

sub <- select_time_window(events_df, start_frac = START_FRAC, end_frac = END_FRAC)
n_window_events <- nrow(sub)
n_sub <- max(c(sub$sender, sub$receiver))

orig_label <- seq_len(n_sub)
if (length(EXCLUDE_ACTORS) > 0) {
  res <- exclude_actors(sub, n_sub, actors_to_exclude = EXCLUDE_ACTORS)
  sub <- res$events_df
  orig_label <- as.integer(names(res$old_to_new))
}

# Repeated until every remaining actor both sends and receives (an actor
# active in one role only contributes no information to its other position).
fa <- filter_active_both(sub)
sub <- fa$events_df; n_sub <- fa$n
actor_label <- orig_label[fa$old_id]     # label of each actor in the original data file

N_out <- tabulate(sub$sender, nbins = n_sub); N_in <- tabulate(sub$receiver, nbins = n_sub)
cat(sprintf("Analysis sample: %d actors, %d events over %.1f days (window: %d events)\n",
            n_sub, nrow(sub), diff(range(sub$time)), n_window_events))
cat(sprintf("Events sent per actor: min %d, median %.0f, max %d\n", min(N_out), median(N_out), max(N_out)))

# first remaining row of each message (recomputed after the filters above)
rate_event <- if (MULTI_RECIPIENT == "message") !duplicated(sub$msg_id) else NULL
if (!is.null(rate_event)) cat(sprintf("Rate model on %d messages (choice model on %d recipient-events)\n", sum(rate_event), nrow(sub)))

# Statistics are computed within the window: counts start at zero at its beginning.
sl <- compute_stats(sub, n_sub, stats = STATS)

# =============================================================================
# 3. Model fits
#    M2 : latent positions + additive link in (pop, rho)      (main model)
#    M2r: latent positions + link in r, same positions          (comparison)
#    M1 : latent positions, no link       M0: no latent positions
# =============================================================================

t0 <- Sys.time()
fit2 <- lsdynam_fit(sub, n_sub, D = D, stats = STATS, link = LINK, df = DF, winsorize_q = WINSORIZE_Q,
                    sigma_a = SIGMA_A, sampler = SAMPLER, n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
                    n_starts = N_STARTS, n_cut = N_CUT, rate_event = rate_event, sl = sl, seed = 1)
minutes_M2 <- as.numeric(Sys.time() - t0, units = "mins")
cat("M2 took", round(minutes_M2, 1), "min\n")

# -----------------------------------------------------------------------------
# 3b. Starts of the choice module. convergence = 1 means that BFGS stopped at
#     its iteration limit; such a start is not a mode. Among converged
#     starts, check that the best objective is reached more than once and
#     that the link does not depend on the start.
# -----------------------------------------------------------------------------
cat("\n=== Starts of the choice module (best first): link at each solution ===\n")
pdf(paste0(OUT_PREFIX, "_modes.pdf"), width = 11, height = 5.5)
modes <- mode_stability_check(fit2)
dev.off()
n_best <- sum(modes$table$delta < 0.5 & modes$table$convergence == 0)
cat(sprintf("Starts within 0.5 of the best objective and converged: %d of %d\n", n_best, nrow(modes$table)))
if (n_best < 2) cat("The best solution was found only once: increase N_STARTS before relying on it.\n")

fit2r <- lsdynam_fit(sub, n_sub, D = D, stats = STATS, link = LINK_R, df = DF, winsorize_q = WINSORIZE_Q,
                     sigma_a = SIGMA_A, sampler = SAMPLER, n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
                     n_cut = N_CUT, rate_event = rate_event, choice = fit2$choice, sl = sl, seed = 1,
                     verbose = FALSE)
fit1 <- lsdynam_fit(sub, n_sub, D = D, stats = STATS, link = "none", sigma_a = SIGMA_A, sampler = SAMPLER,
                    n_iter = N_ITER, burn_in = BURN_IN, thin = THIN, rate_event = rate_event,
                    choice = fit2$choice, sl = sl, seed = 1, verbose = FALSE)
fit0 <- lsdynam_fit(sub, n_sub, D = 0, stats = STATS, link = "none", sigma_a = SIGMA_A, sampler = SAMPLER,
                    n_iter = N_ITER, burn_in = BURN_IN, thin = THIN, rate_event = rate_event,
                    sl = sl, seed = 1, verbose = FALSE)

# =============================================================================
# 4. Posterior summaries and convergence
# =============================================================================

cat("\n=== M2 (additive link) ===\n"); print(fit2)
cat("\n=== M2r (link in r, same positions) ===\n"); print(fit2r)
cat("\n=== M1 ===\n"); print(fit1)
cat("\n=== M0 ===\n"); print(fit0)

cat("\n=== Convergence (M2) ===\n")
dg <- mcmc_diagnostics(fit2)
print(summary(fit2))

pdf(paste0(OUT_PREFIX, "_traces.pdf"), width = 10, height = 7)
plot_traces(fit2)                                   # coefficients and log-posterior of the choice module
plot_traces(fit2, which = "activity", actors = order(N_out)[seq_len(min(12, n_sub))])   # the least active senders
dev.off()

# =============================================================================
# 5. Figures
# =============================================================================

ps <- fit2$summaries
a_M1 <- colMeans(fit1$rate$a)
cex_N <- 0.5 + 1.2 * sqrt(N_out / max(N_out))

pdf(paste0(OUT_PREFIX, "_fig_positions.pdf"), width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
if (D == 1) plot_positions(fit2, legend_pos = "bottomright") else plot(ps$pop, ps$rho, pch = ifelse(ps$r > 0, 19, 21), cex = cex_N,
                                           col = "black", bg = "white",
                                           xlab = "latent popularity", ylab = "excess self-affinity",
                                           main = "Invariant summaries (filled: r > 0; open: r < 0)")
plot_activity_effects(fit2, level = 0.95)
par(op); dev.off()

# page 1: activity effects of M1 (no link) against the two components
# page 2: levels of the components (shaded: plug-in band; dashed: cut posterior)
# page 3: local slopes of the components
pdf(paste0(OUT_PREFIX, "_fig_components.pdf"), width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
plot(ps$pop, a_M1, pch = 19, cex = cex_N, col = "grey30", xlab = "latent popularity",
     ylab = "activity effect, model without link", main = "Activity against popularity")
abline(v = 0, col = "grey40", lty = 2)
plot(ps$rho, a_M1, pch = 19, cex = cex_N, col = "grey30", xlab = "excess self-affinity",
     ylab = "activity effect, model without link", main = "Activity against excess self-affinity")
abline(v = 0, col = "grey40", lty = 2)
j_pop <- which(LINK == "pop"); j_rho <- which(LINK == "rho")
plot_link_level(fit2, component = j_pop, xlab = "latent popularity", main = sprintf("Component in latent popularity (D = %d)", D),
                legend_pos = NULL)
plot_link_level(fit2, component = j_rho, xlab = "excess self-affinity", main = sprintf("Component in excess self-affinity (D = %d)", D),
                legend_pos = NULL)
plot_link_slope(fit2, component = j_pop, xlab = "latent popularity", main = "Local slope, popularity")
plot_link_slope(fit2, component = j_rho, xlab = "excess self-affinity", main = "Local slope, excess self-affinity")
par(op); dev.off()

pdf(paste0(OUT_PREFIX, "_fig_link_r.pdf"), width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
plot_link_level(fit2r, xlab = "self-affinity", main = "Level of the link in r", legend_pos = NULL)
plot_link_slope(fit2r, xlab = "self-affinity", main = "Local slope of the link in r")
par(op); dev.off()

cat("\n=== Invariant summaries of the positions (plug-in) ===\n")
sv <- c("pop", "rho", "r", "z", "phi", "delta")
print(round(t(sapply(ps[sv], function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x)))), 3))
cat("Correlations:\n"); print(round(cor(ps[sv]), 2))
cat(sprintf("Actors with positive self-affinity: %d of %d\n", sum(ps$r > 0), n_sub))
quad <- NULL
if (D == 1) {
  # common sign of the positions as displayed by plot_positions (mean sender position positive)
  sg <- if (mean(fit2$choice$U_mean[, 1]) < 0) -1 else 1
  u <- sg * fit2$choice$U_mean[, 1]; dv <- sg * (fit2$choice$V_mean[, 1] - mean(fit2$choice$V_mean[, 1]))
  quad <- c(I = sum(u > 0 & dv > 0), II = sum(u < 0 & dv > 0), III = sum(u < 0 & dv < 0), IV = sum(u > 0 & dv < 0))
  cat(sprintf("Sender positions (sign as displayed): u_i > 0 for %d of %d actors (ubar = %.3f)\n", sum(u > 0), n_sub, mean(u)))
  cat("Quadrants (u, v - vbar):\n"); print(quad)
}

# =============================================================================
# 6. Evidence for the link and for each component given the other
# =============================================================================

evd   <- link_evidence(fit2, sigma_kappa_grid = SK_GRID)
evd_c <- lapply(setNames(LINK, LINK), function(v) link_evidence(fit2, components = v, sigma_kappa_grid = SK_GRID))
evd_r <- link_evidence(fit2r, sigma_kappa_grid = SK_GRID)
cat("\n=== Evidence (log Bayes factors; sigma_a at its value under each link) ===\n")
ev_tab <- rbind(
  "whole additive link (M2 against M1)" = c(sd = evd$log_bf10_sd, laplace = evd$log_bf10_laplace),
  "excess self-affinity given popularity" = c(evd_c$rho$log_bf10_sd, evd_c$rho$log_bf10_laplace),
  "popularity given excess self-affinity" = c(evd_c$pop$log_bf10_sd, evd_c$pop$log_bf10_laplace),
  "link in r (M2r against M1)" = c(evd_r$log_bf10_sd, evd_r$log_bf10_laplace))
print(round(ev_tab, 2))
cat("Laplace log Bayes factors for other prior standard deviations of kappa:\n")
ev_grid <- rbind(whole = evd$grid, rho_given_pop = evd_c$rho$grid, pop_given_rho = evd_c$pop$grid, r = evd_r$grid)
print(round(ev_grid, 2))
cat(sprintf("Residual sd of the activity effects: M1 %.3f, M2 %.3f, M2r %.3f ; R2_a: M2 %.3f, M2r %.3f\n",
            fit1$rate$sigma_a, fit2$rate$sigma_a, fit2r$rate$sigma_a, link_r2(fit2), link_r2(fit2r)))

contrast_rows <- function(f, label) {
  vars <- f$rate$link$vars
  do.call(rbind, lapply(seq_along(vars), function(j) {
    x <- f$rate$link$X[, j]
    out <- list(link_contrast(f, component = j, source = "plugin"),
                if (!is.null(f$rate$cut)) link_contrast(f, component = j, source = "cut"),
                link_contrast(f, component = j, x_lo = mean(x) - sd(x), x_hi = mean(x) + sd(x)))
    data.frame(link = label, component = vars[j],
               what = c("10-90, plug-in", "10-90, cut", "mean -/+ 1 sd, plug-in")[!vapply(out, is.null, logical(1))],
               mean = vapply(Filter(Negate(is.null), out), function(o) o$mean, numeric(1)),
               lo = vapply(Filter(Negate(is.null), out), function(o) o$ci[1], numeric(1)),
               hi = vapply(Filter(Negate(is.null), out), function(o) o$ci[2], numeric(1)),
               p_positive = vapply(Filter(Negate(is.null), out), function(o) o$prob_positive, numeric(1)),
               stringsAsFactors = FALSE)
  }))
}
cat("\nContrasts of the level (level at the higher value minus level at the lower value):\n")
contr <- rbind(contrast_rows(fit2, "M2"), contrast_rows(fit2r, "M2r"))
print(format(contr, digits = 3), row.names = FALSE)

# =============================================================================
# 7. Links in different summaries, on the same positions (Laplace version of
#    the rate module; sigma_a by maximum marginal likelihood in each model)
# =============================================================================

specs <- list("pop + rho" = LINK, "pop" = "pop", "rho" = "rho", "r" = "r",
              "pop + rho, linear" = list(link = LINK, link_type = "linear"),
              "r, linear" = list(link = "r", link_type = "linear"),
              "z" = "z", "delta" = "delta", "phi" = "phi")
if (D >= 2) specs[["log phi + log delta + align"]] <- c("log_phi", "log_delta", "align")
cmp <- link_comparison_table(fit2, specs)
cat(sprintf("\n=== Links in different summaries (log ML relative to M1, for which sigma_a = %.3f) ===\n", cmp$sigma_a_M1))
print(format(cmp$table, digits = 3), row.names = FALSE)
ml <- setNames(cmp$table$log_ml_vs_M1, cmp$table$link)
cmp_derived <- c(rho_given_pop = unname(ml["pop + rho"] - ml["pop"]),
                 pop_given_rho = unname(ml["pop + rho"] - ml["rho"]),
                 additive_minus_r = unname(ml["pop + rho"] - ml["r"]),
                 spline_minus_linear_additive = unname(ml["pop + rho"] - ml["pop + rho, linear"]),
                 spline_minus_linear_r = unname(ml["r"] - ml["r, linear"]))
print(round(cmp_derived, 2))

# =============================================================================
# 8. Diagnostics
# =============================================================================

pdf(paste0(OUT_PREFIX, "_diagnostics.pdf"), width = 12, height = 6)
cat("\n=== Summaries against the number of events sent ===\n")
sh <- link_shrinkage_check(fit2, vars = c("pop", "rho", "r", "z"))
cat("\n=== Leverage of each component ===\n")
lev <- lapply(setNames(seq_along(LINK), LINK), function(j) link_leverage_check(fit2, sub, n_exclude = 5, component = j))
lev_r <- link_leverage_check(fit2r, sub, n_exclude = 5)
cat("\n=== Plug-in summaries against the posterior means of the summaries ===\n")
ag <- summary_agreement_check(fit2, vars = c("pop", "rho", "r", "norm_u"))
dev.off()

# =============================================================================
# 9. Sensitivity of the additive link (same positions, Laplace version of
#    stage 2)
# =============================================================================

sens_specs <- list(
  "pop + rho, spline (main)"              = list(link = LINK),
  "pop + rho, linear"                     = list(link = LINK, link_type = "linear"),
  "pop + rho, K = 3 per component"        = list(link = LINK, df = 3),
  "pop + rho, K = 6 per component"        = list(link = LINK, df = 6),
  "pop + rho, winsorised at 5%/95%"       = list(link = LINK, winsorize_q = 0.95),
  "pop + rho, Student-t (nu = 4)"         = list(link = LINK, nu = 4),
  "pop + rho, sigma_a = 1"                = list(link = LINK, sigma_a = 1),
  "pop + rho, unscaled basis"             = list(link = LINK, scale_basis = FALSE),
  "pop + rho, sigma_kappa by marg. lik."  = list(link = LINK, sigma_kappa = "eb"),
  "pop + rho, posterior-mean summaries"   = list(link = LINK, plugin = "summary"),
  "r, spline"                             = list(link = "r"),
  "uv = u'v (uncentred)"                  = list(link = "uv"),
  "norm_u = ||u||^2"                      = list(link = "norm_u")
)
sens <- lapply(sens_specs, function(sp) {
  args <- c(list(events_df = sub, n = n_sub, D = D, stats = STATS, df = DF, winsorize_q = WINSORIZE_Q,
                 sigma_a = SIGMA_A, inference = "map", choice = fit2$choice, rate_event = rate_event, sl = sl,
                 verbose = FALSE), sp)
  args <- args[!duplicated(names(args), fromLast = TRUE)]
  do.call(lsdynam_fit, args)
})
sens_tab <- do.call(rbind, lapply(names(sens), function(nm) {
  f <- sens[[nm]]; vars <- f$rate$link$vars; e <- link_evidence(f)
  ct <- lapply(seq_along(vars), function(k) link_contrast(f, component = k))
  data.frame(link = nm, K = e$K, sigma_a = f$rate$sigma_a, sigma_kappa = f$rate$sigma_kappa,
             log_bf10 = e$log_bf10_laplace,
             log_bf10_rho_given_pop = if (all(LINK %in% vars)) link_evidence(f, components = "rho")$log_bf10_laplace else NA_real_,
             R2 = link_r2(f),
             contrast = paste(vapply(seq_along(vars), function(k)
               sprintf("%s %.2f [%.2f, %.2f]", vars[k], ct[[k]]$mean, ct[[k]]$ci[1], ct[[k]]$ci[2]), character(1)),
               collapse = "; "),
             stringsAsFactors = FALSE)
}))
cat("\n=== Sensitivity of the link (log Bayes factors at the sigma_a of each fit) ===\n")
print(format(sens_tab, digits = 3), row.names = FALSE)

pdf(paste0(OUT_PREFIX, "_sensitivity.pdf"), width = 14, height = 7)
show <- c("pop + rho, spline (main)", "pop + rho, winsorised at 5%/95%", "pop + rho, Student-t (nu = 4)",
          "pop + rho, unscaled basis")
op <- par(mfrow = c(2, length(show)))
for (j in seq_along(LINK)) for (nm in show) plot_link_level(sens[[nm]], component = j, show_cut = FALSE,
                                                             main = paste0(nm, ": ", LINK[j]), cex.main = 0.8)
par(op); dev.off()

# =============================================================================
# 10. Held-out predictive comparison (fit on the first part of the window)
# =============================================================================

holdout <- NULL
if (RUN_HOLDOUT) {
  m_tr <- round(TRAIN_FRAC * nrow(sub)); sub_tr <- sub[seq_len(m_tr), ]
  re_tr <- if (is.null(rate_event)) NULL else rate_event[seq_len(m_tr)]
  cat(sprintf("\n=== Held-out scores: %d training events, %d test events ===\n", m_tr, nrow(sub) - m_tr))
  h2 <- lsdynam_fit(sub_tr, n_sub, D = D, stats = STATS, link = LINK, df = DF, winsorize_q = WINSORIZE_Q,
                    sigma_a = SIGMA_A, inference = "map", n_starts = N_STARTS, rate_event = re_tr, seed = 1, verbose = FALSE)
  h2r <- lsdynam_fit(sub_tr, n_sub, D = D, stats = STATS, link = LINK_R, df = DF, winsorize_q = WINSORIZE_Q,
                     sigma_a = SIGMA_A, inference = "map", choice = h2$choice, rate_event = re_tr, seed = 1, verbose = FALSE)
  h1 <- lsdynam_fit(sub_tr, n_sub, D = D, stats = STATS, link = "none", sigma_a = SIGMA_A, inference = "map",
                    choice = h2$choice, rate_event = re_tr, seed = 1, verbose = FALSE)
  h0 <- lsdynam_fit(sub_tr, n_sub, D = 0, stats = STATS, link = "none", sigma_a = SIGMA_A, inference = "map",
                    rate_event = re_tr, seed = 1, verbose = FALSE)
  sc <- lapply(list(M0 = h0, M1 = h1, M2 = h2, M2r = h2r), predictive_scores, events_full = sub, n = n_sub,
               m_train = m_tr, n_draws = 200, rate_event = rate_event)
  # rate: sum over actors of log predictive densities; choice: sum over held-out events
  print(round(sapply(sc, function(s) c(rate = s$rate, choice = s$choice, choice_per_event = s$choice / s$n_test)), 4))
  N_tr <- sc$M1$N_out_train
  grp <- cut(N_tr, breaks = unique(quantile(N_tr, c(0, 0.25, 0.5, 0.75, 1))), include.lowest = TRUE)
  cat("\nDifferences with standard errors (rate: over actors, also by training out-degree; choice: over events):\n")
  hd <- rbind(cbind(comparison = "M1 - M0", score_difference(sc$M1, sc$M0)),
              cbind(comparison = "M2 - M1", score_difference(sc$M2, sc$M1, groups = grp)),
              cbind(comparison = "M2r - M1", score_difference(sc$M2r, sc$M1, groups = grp)))
  print(hd, digits = 3, row.names = FALSE)
  holdout <- list(scores = sc, groups = grp, m_train = m_tr, differences = hd,
                  uniform_choice_per_event = -log(n_sub - 1))
}

settings <- list(DATA_FILE = DATA_FILE, START_FRAC = START_FRAC, END_FRAC = END_FRAC,
                 EXCLUDE_ACTORS = EXCLUDE_ACTORS, D = D, STATS = STATS, LINK = LINK, DF = DF,
                 SIGMA_A = SIGMA_A, WINSORIZE_Q = WINSORIZE_Q, MULTI_RECIPIENT = MULTI_RECIPIENT,
                 N_ITER = N_ITER, BURN_IN = BURN_IN, THIN = THIN, N_STARTS = N_STARTS,
                 N_CUT = N_CUT, N_NULL = N_NULL)
# interim save: the calibrations below are slow, and the fits should not be
# lost if they are interrupted (the final save overwrites this file)
saveRDS(list(fit2 = fit2, fit2r = fit2r, fit1 = fit1, fit0 = fit0, comparison = cmp$table,
             sensitivity = sens_tab, modes = modes$table, holdout = holdout, actor_label = actor_label,
             settings = settings), paste0(OUT_PREFIX, "_results.rds"))

# =============================================================================
# 11. Calibration by simulation (Laplace version, one random start plus the
#     spectral start per replicate)
#   no association : activity effects randomly reassigned to the actors
#   popularity only: activity effects = fitted link in popularity plus
#                    reassigned residuals; calibrates the evidence for excess
#                    self-affinity given popularity
# The replicates use the summaries at the mode of the choice module, so the
# observed values to compare with are those of the best start in the modes
# table (also printed).
# =============================================================================

calib_print <- function(nc, label) {
  cat(sprintf("\n%s (%d replicates):\n", label, length(nc$r2)))
  q <- function(x) sprintf("median %6.2f  min %6.2f  max %6.2f%s", median(x, na.rm = TRUE), min(x, na.rm = TRUE),
                           max(x, na.rm = TRUE), if (anyNA(x)) sprintf("  (%d failed)", sum(is.na(x))) else "")
  cat("  R2_a                   ", q(nc$r2), "\n")
  cat("  log BF10, whole link   ", q(nc$log_bf10), "\n")
  if (!is.null(nc$log_bf10_components)) for (v in colnames(nc$log_bf10_components))
    cat(sprintf("  log BF10, %-4s | rest   %s\n", v, q(nc$log_bf10_components[, v])))
  for (v in colnames(nc$contrast)) cat(sprintf("  contrast 10-90, %-6s %s\n", v, q(nc$contrast[, v])))
}
obs_mode <- modes$table[1, ]
null_cal <- NULL; null_pop <- NULL
if (RUN_NULL) {
  cat("\n=== Calibration: no association ===\n")
  null_cal <- null_link_calibration(fit2, n_rep = N_NULL, n_starts = 1)
  calib_print(null_cal, "No association")
}
if (RUN_NULL_POP) {
  cat("\n=== Calibration: association with popularity only ===\n")
  null_pop <- null_link_calibration(fit2, n_rep = N_NULL, n_starts = 1, keep = "pop", seed = 1001)
  calib_print(null_pop, "Popularity only")
}
if (RUN_NULL || RUN_NULL_POP) {
  cat("\nObserved, at the best mode (comparable with the replicates):\n")
  print(format(obs_mode, digits = 3), row.names = FALSE)
  # observed link with the summaries at the mode of the choice module, as in the replicates
  ch_mode <- fit2$choice; ch_mode$U_mean <- ch_mode$map$U; ch_mode$V_mean <- ch_mode$map$V
  fit_mode <- fit2
  fit_mode$rate <- lsdynam_rate_stage(NULL, n_sub, sl = NULL, choice = ch_mode, link = LINK, df = DF,
                                      winsorize_q = WINSORIZE_Q, priors = fit2$priors, sigma_a = SIGMA_A,
                                      inference = "map", cr = fit2$rate$cr, seed = 1)
  pdf(paste0(OUT_PREFIX, "_calibration.pdf"), width = 11, height = 5.5)
  for (nc in Filter(Negate(is.null), list(null_cal, null_pop))) {
    op <- par(mfrow = c(1, length(LINK)))
    for (j in seq_along(LINK)) {
      obs <- link_curve(fit_mode, j, grid = nc$grid[[j]])
      matplot(nc$grid[[j]], t(nc$levels[[j]]), type = "l", lty = 1, col = "grey70",
              ylim = range(nc$levels[[j]], obs$level_mean - obs$mu_a_mean, na.rm = TRUE), xlab = LINK[j],
              ylab = "component, relative to mu_a",
              main = sprintf("%s: observed at the mode against replicates (%s)", LINK[j],
                             if (is.null(nc$keep)) "no association" else "popularity only"))
      lines(nc$grid[[j]], obs$level_mean - obs$mu_a_mean, lwd = 3, col = "black")
      abline(v = 0, col = "grey40", lty = 2)
    }
    par(op)
  }
  dev.off()
}

# =============================================================================
# 12. Numbers for the paper
# =============================================================================

q3 <- function(v) c(mean = mean(v), lo = quantile(v, 0.025, names = FALSE), hi = quantile(v, 0.975, names = FALSE))
paper <- list(
  sample = c(window_events = n_window_events, n = n_sub, M = nrow(sub), days = diff(range(sub$time)),
             N_out_min = min(N_out), N_out_median = median(N_out), N_out_max = max(N_out),
             excluded = length(EXCLUDE_ACTORS)),
  mcmc = c(accept_choice = fit2$choice$accept_rate, accept_rate = fit2$rate$accept_rate,
           ess_zeta = min(dg$ess$zeta), ess_pop_min = min(dg$ess$pop), ess_pop_median = median(dg$ess$pop),
           ess_rho_min = min(dg$ess$rho), ess_rho_median = median(dg$ess$rho),
           ess_r_min = min(dg$ess$r), ess_kappa_min = min(dg$ess$kappa), ess_a_min = min(dg$ess$a),
           minutes_M2 = minutes_M2, seconds_choice = fit2$choice$seconds, seconds_rate = fit2$rate$seconds,
           starts_at_best = n_best, starts = nrow(modes$table)),
  modes = modes$table,
  coefficients = rbind(zeta_M2 = q3(fit2$choice$zeta[, 1]), zeta_M0 = q3(fit0$choice$zeta[, 1]),
                       xi_M2 = q3(fit2$rate$xi[, 1])),
  positions = c(u_positive = if (D == 1) sum(sg * fit2$choice$U_mean[, 1] > 0) else NA,
                quad, r_positive = sum(ps$r > 0), r_min = min(ps$r), r_max = max(ps$r), r_sd = sd(ps$r),
                pop_sd = sd(ps$pop), rho_sd = sd(ps$rho), cor_pop_rho = cor(ps$pop, ps$rho),
                cor_pop_delta = cor(ps$pop, ps$delta), cor_r_pop = cor(ps$r, ps$pop)),
  sigma_a = c(M1 = fit1$rate$sigma_a, M2 = fit2$rate$sigma_a, M2r = fit2r$rate$sigma_a),
  R2 = c(M2 = link_r2(fit2), M2r = link_r2(fit2r)),
  log_bf = ev_tab, log_bf_grid = ev_grid,
  contrasts = contr,
  comparison = cmp$table, comparison_derived = cmp_derived,
  shrinkage = sh,
  leverage_max_change = c(sapply(lev, function(l) l$max_diff), r = lev_r$max_diff),
  agreement = ag$table,
  sensitivity = sens_tab,
  holdout = if (!is.null(holdout)) list(per_event = sapply(holdout$scores, function(s) s$choice / s$n_test),
                                        uniform = holdout$uniform_choice_per_event,
                                        differences = holdout$differences),
  calibration_none = if (!is.null(null_cal)) list(r2 = summary(null_cal$r2), log_bf10 = summary(null_cal$log_bf10),
                                                  rho_given_pop = summary(null_cal$log_bf10_components[, "rho"]),
                                                  pop_given_rho = summary(null_cal$log_bf10_components[, "pop"]),
                                                  contrast_pop = summary(null_cal$contrast[, "pop"]),
                                                  contrast_rho = summary(null_cal$contrast[, "rho"])),
  calibration_pop_only = if (!is.null(null_pop)) list(rho_given_pop = summary(null_pop$log_bf10_components[, "rho"]),
                                                      contrast_rho = summary(null_pop$contrast[, "rho"]),
                                                      pop_given_rho = summary(null_pop$log_bf10_components[, "pop"])),
  observed_at_mode = obs_mode
)
print_paper <- function(p) {
  for (nm in names(p)) {
    if (is.null(p[[nm]])) next
    cat(sprintf("\n--- %s ---\n", nm))
    x <- p[[nm]]
    if (is.data.frame(x)) print(format(x, digits = 3), row.names = FALSE)
    else if (is.list(x)) for (k in names(x)) { cat(sprintf("%s:\n", k)); print(if (is.numeric(x[[k]])) round(x[[k]], 3) else x[[k]]) }
    else print(round(x, 3))
  }
}
cat("\n================ NUMBERS FOR THE PAPER ================\n")
print_paper(paper)
capture.output(print_paper(paper), file = paste0(OUT_PREFIX, "_paper_numbers.txt"))

saveRDS(list(fit2 = fit2, fit2r = fit2r, fit1 = fit1, fit0 = fit0, comparison = cmp$table,
             sensitivity = sens_tab, modes = modes$table, holdout = holdout,
             null_cal = null_cal, null_pop = null_pop, paper = paper, actor_label = actor_label,
             settings = settings),
        paste0(OUT_PREFIX, "_results.rds"))
cat("\nDone.\n")
