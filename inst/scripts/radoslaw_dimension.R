# =============================================================================
# radoslaw_dimension.R
#
# Latent spaces of dimension D > 1 (Sections 8.7 and 8.8 of the paper), with
# the additive link in latent popularity and excess self-affinity of the main
# analysis. D = 2 (D_DETAIL) is analysed in detail, since with two dimensions
# the senders can favour different receivers and excess self-affinity has
# room to matter; D = 3 and D = 4 give the essential results only.
#
#   F1. Held-out choice score for D = 2, 3, 4 on the split of the main
#       analysis, with the differences between successive dimensions and
#       their standard errors (models without link, Laplace version).
#   F2. M2 with the additive link for D = 2 and D = 3 on all events, with the
#       settings of the main analysis: evidence for the link and for each
#       component given the other, contrasts of the two components, links in
#       pop, rho, pop + rho, r and (log phi, log delta, align) on the same
#       positions, correlation of the summaries with those for D = 1, starts
#       of the choice module. For D = D_DETAIL also the cut posterior, the
#       Bayes factors for several prior scales, the model without link (M1)
#       and a three-page figure (components, activity against the summaries,
#       local slopes).
#   F4. Calibration for D = D_DETAIL with the association with popularity
#       only: reference distribution of the evidence for excess
#       self-affinity given popularity, including what the shrinkage of the
#       sender positions produces on its own (and, optionally, with no
#       association at all).
#   F3. Table of the essential results for D = 1, 2, 3 and the detailed
#       numbers for D = D_DETAIL (written to radoslaw_dimension_essentials.txt).
#
# Reads radoslaw_results.rds written by radoslaw_application.R (same sample
# and settings: run that script first). Results in
# radoslaw_dimension_results.rds, figures in radoslaw_dimension_*.pdf.
# Running time on the earlier sample: F1 about 40 minutes, F2 about 25
# minutes for D = 2 and 30 for D = 3, F4 about 2 hours (20 replicates with
# D = 2). RUN_F4 <- FALSE before sourcing skips F4.
# =============================================================================

library(ls.dynam)

if (!exists("RESULTS_FILE")) RESULTS_FILE <- "radoslaw_results.rds"
saved <- readRDS(RESULTS_FILE)
st    <- saved$settings
dflt  <- function(x, d) if (is.null(x)) d else x

# defaults: the settings of the main analysis, saved in RESULTS_FILE
if (!exists("DATA_FILE"))    DATA_FILE    <- dflt(st$DATA_FILE, "out_radoslaw.txt")
if (!exists("RUN_F1"))   RUN_F1   <- TRUE          # held-out choice score for D_SCORE
if (!exists("RUN_F2"))   RUN_F2   <- TRUE          # M2 with the additive link for D_LINK
if (!exists("D_SCORE"))  D_SCORE  <- c(2, 3, 4)
if (!exists("D_LINK"))   D_LINK   <- c(2, 3)
if (!exists("D_DETAIL")) D_DETAIL <- 2             # dimension analysed in detail
if (!exists("RUN_F4"))   RUN_F4   <- TRUE          # calibration with popularity only, D = D_DETAIL
if (!exists("RUN_F4_NONE")) RUN_F4_NONE <- FALSE   # calibration with no association, D = D_DETAIL
if (!exists("N_NULL"))   N_NULL   <- dflt(st$N_NULL, 20)
if (!exists("N_CUT"))    N_CUT    <- dflt(st$N_CUT, 300)
if (!exists("N_STARTS")) N_STARTS <- dflt(st$N_STARTS, 10)
if (!exists("N_ITER"))   N_ITER   <- dflt(st$N_ITER, 12000)
if (!exists("BURN_IN"))  BURN_IN  <- dflt(st$BURN_IN, 2000)
if (!exists("THIN"))     THIN     <- dflt(st$THIN, 5)
OUT_PREFIX <- "radoslaw_dimension"
DF    <- st$DF
WINSORIZE_Q <- dflt(st$WINSORIZE_Q, 1)
LINK  <- if (!is.null(st$LINK)) st$LINK else c("pop", "rho")
if (!all(c("pop", "rho") %in% LINK)) stop("radoslaw_results.rds was not produced with the additive link in pop and rho")
out   <- list()

# ---- the sample of the main analysis ------------------------------------------

rebuild_sample <- function() {
  raw <- read.table(DATA_FILE, header = FALSE, comment.char = "%",
                    col.names = c("sender", "receiver", "weight", "timestamp"))
  ev <- data.frame(sender = raw$sender, receiver = raw$receiver, time = as.numeric(raw$timestamp))
  ev <- ev[order(ev$time), ]
  ev <- ev[ev$sender != ev$receiver, ]
  key <- paste(ev$sender, ev$time); ev$msg_id <- match(key, unique(key))
  ev$time <- (ev$time - min(ev$time)) / 86400
  if (sum(duplicated(ev$time)) > 0) ev$time <- ev$time + seq(0, 1e-6, length.out = nrow(ev))
  s <- select_time_window(ev, start_frac = st$START_FRAC, end_frac = st$END_FRAC)
  ns <- max(c(s$sender, s$receiver))
  if (length(st$EXCLUDE_ACTORS) > 0) s <- exclude_actors(s, ns, actors_to_exclude = st$EXCLUDE_ACTORS)$events_df
  fa <- filter_active_both(s)
  list(sub = fa$events_df, n = fa$n,
       rate_event = if (st$MULTI_RECIPIENT == "message") !duplicated(fa$events_df$msg_id) else NULL)
}

smp <- rebuild_sample()
stopifnot(smp$n == saved$fit2$n, nrow(smp$sub) == saved$fit2$M)

specs <- list("pop + rho" = LINK, "pop" = "pop", "rho" = "rho", "r" = "r")
SK_GRID <- c(0.1, 0.25, 0.5, 1, 2)

contrast_rows <- function(f) {
  vars <- f$rate$link$vars
  do.call(rbind, lapply(seq_along(vars), function(j) {
    x <- f$rate$link$X[, j]
    out <- Filter(Negate(is.null), list(
      plugin = link_contrast(f, component = j, source = "plugin"),
      cut = if (!is.null(f$rate$cut)) link_contrast(f, component = j, source = "cut"),
      sd = link_contrast(f, component = j, x_lo = mean(x) - sd(x), x_hi = mean(x) + sd(x))))
    data.frame(component = vars[j],
               what = c(plugin = "10-90, plug-in", cut = "10-90, cut", sd = "mean -/+ 1 sd, plug-in")[names(out)],
               mean = vapply(out, function(o) o$mean, numeric(1)), lo = vapply(out, function(o) o$ci[1], numeric(1)),
               hi = vapply(out, function(o) o$ci[2], numeric(1)),
               p_positive = vapply(out, function(o) o$prob_positive, numeric(1)), stringsAsFactors = FALSE, row.names = NULL)
  }))
}
calib_print <- function(nc, label) {
  cat(sprintf("\n%s (%d replicates):\n", label, length(nc$r2)))
  q <- function(x) sprintf("median %6.2f  min %6.2f  max %6.2f%s", median(x, na.rm = TRUE), min(x, na.rm = TRUE),
                           max(x, na.rm = TRUE), if (anyNA(x)) sprintf("  (%d failed)", sum(is.na(x))) else "")
  cat("  R2_a                   ", q(nc$r2), "\n")
  cat("  log BF10, whole link   ", q(nc$log_bf10), "\n")
  for (v in colnames(nc$log_bf10_components)) cat(sprintf("  log BF10, %-4s | rest   %s\n", v, q(nc$log_bf10_components[, v])))
  for (v in colnames(nc$contrast)) cat(sprintf("  contrast 10-90, %-6s %s\n", v, q(nc$contrast[, v])))
}

# Essential quantities of an M2 fit with the additive link
essentials <- function(f, cmp) {
  ml <- setNames(cmp$table$log_ml_vs_M1, cmp$table$link)
  e_all <- link_evidence(f); e_rho <- link_evidence(f, components = "rho")
  c_pop <- link_contrast(f, component = which(f$rate$link$vars == "pop"))
  c_rho <- link_contrast(f, component = which(f$rate$link$vars == "rho"))
  c(zeta = mean(f$choice$zeta[, 1]), sigma_a_M1 = cmp$sigma_a_M1, sigma_a_M2 = f$rate$sigma_a,
    R2 = link_r2(f), R2_pop_alone = cmp$table$R2[cmp$table$link == "pop"], log_bf_link = e_all$log_bf10_laplace,
    log_bf_rho_given_pop_sd = e_rho$log_bf10_sd, log_bf_rho_given_pop = e_rho$log_bf10_laplace,
    logml_pop = unname(ml["pop"]), logml_rho = unname(ml["rho"]), logml_pop_rho = unname(ml["pop + rho"]),
    logml_r = unname(ml["r"]), rho_given_pop = unname(ml["pop + rho"] - ml["pop"]),
    pop_given_rho = unname(ml["pop + rho"] - ml["rho"]),
    contrast_pop = c_pop$mean, contrast_pop_lo = c_pop$ci[1], contrast_pop_hi = c_pop$ci[2],
    contrast_rho = c_rho$mean, contrast_rho_lo = c_rho$ci[1], contrast_rho_hi = c_rho$ci[2])
}

# =============================================================================
# F1. Held-out choice score for higher dimensions
# =============================================================================

if (RUN_F1 && !is.null(saved$holdout)) {
  cat(sprintf("\n=== F1. Held-out choice score for D = %s ===\n", paste(D_SCORE, collapse = ", ")))
  m_tr <- saved$holdout$m_train
  sub_tr <- smp$sub[seq_len(m_tr), ]
  re_tr <- if (is.null(smp$rate_event)) NULL else smp$rate_event[seq_len(m_tr)]
  scD <- list(D0 = saved$holdout$scores$M0, D1 = saved$holdout$scores$M1)
  for (Dk in D_SCORE) {
    t0 <- Sys.time()
    h <- lsdynam_fit(sub_tr, smp$n, D = Dk, stats = st$STATS, link = "none", sigma_a = "eb", inference = "map",
                     n_starts = N_STARTS, rate_event = re_tr, seed = 1, verbose = TRUE)
    scD[[paste0("D", Dk)]] <- predictive_scores(h, smp$sub, smp$n, m_tr, n_draws = 200, rate_event = smp$rate_event)
    cat(sprintf("  D = %d done in %.1f min\n", Dk, as.numeric(Sys.time() - t0, units = "mins")))
  }
  dims <- as.integer(sub("D", "", names(scD)))
  scD <- scD[order(dims)]; dims <- sort(dims)
  per_event <- vapply(scD, function(s) sum(s$choice_by_event) / s$n_test, numeric(1))
  cat("\nChoice score per held-out event:\n"); print(round(per_event, 4))
  dF1 <- do.call(rbind, lapply(seq_along(dims)[-1], function(k) {
    x <- score_difference(scD[[k]], scD[[k - 1]])[2, ]
    cbind(comparison = sprintf("D = %d minus D = %d", dims[k], dims[k - 1]), x)
  }))
  print(dF1, digits = 3, row.names = FALSE)
  out$F1 <- list(per_event = per_event, differences = dF1, dims = dims)
}

# =============================================================================
# F2. M2 with the additive link for D = D_LINK on all events
# =============================================================================

ps1 <- saved$fit2$summaries
if (RUN_F2) {
  out$F2 <- list()
  for (Dk in D_LINK) {
    cat(sprintf("\n=== F2. M2 with the additive link, D = %d, all events ===\n", Dk))
    t0 <- Sys.time()
    detail <- Dk == D_DETAIL
    fitD <- lsdynam_fit(smp$sub, smp$n, D = Dk, stats = st$STATS, link = LINK, df = DF,
                        winsorize_q = WINSORIZE_Q, sigma_a = st$SIGMA_A, sampler = "mala", n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
                        n_starts = N_STARTS, n_cut = if (detail) N_CUT else 0, rate_event = smp$rate_event, seed = 1)
    minutes <- as.numeric(Sys.time() - t0, units = "mins")
    cat(sprintf("(%.1f min)\n", minutes))
    print(fitD); dgD <- mcmc_diagnostics(fitD)

    cat("\nStarts of the choice module and the link at each of them:\n")
    pdf(paste0(OUT_PREFIX, "_modes_D", Dk, ".pdf"), width = 11, height = 5.5)
    modesD <- mode_stability_check(fitD)
    dev.off()

    eD <- link_evidence(fitD); eD_rho <- link_evidence(fitD, components = "rho"); eD_pop <- link_evidence(fitD, components = "pop")
    cat(sprintf("log BF10: whole link %.2f (SD %.2f) ; rho given pop %.2f (SD %.2f) ; pop given rho %.2f (SD %.2f) ; R2_a %.3f\n",
                eD$log_bf10_laplace, eD$log_bf10_sd, eD_rho$log_bf10_laplace, eD_rho$log_bf10_sd,
                eD_pop$log_bf10_laplace, eD_pop$log_bf10_sd, link_r2(fitD)))

    psD <- fitD$summaries
    cr1 <- vapply(c("pop", "rho", "r"), function(v) cor(psD[[v]], ps1[[v]]), numeric(1))
    cat("Correlation of the summaries with D = 1:\n"); print(round(cr1, 3))
    cat(sprintf("Correlation of pop and rho: %.3f ; alignment: median %.2f\n", cor(psD$pop, psD$rho), median(psD$align)))

    specsD <- specs
    specsD[["log phi + log delta + align"]] <- c("log_phi", "log_delta", "align")
    cmpD <- link_comparison_table(fitD, specsD)
    cat(sprintf("\nLinks on the same positions (log ML relative to M1, for which sigma_a = %.3f):\n", cmpD$sigma_a_M1))
    print(format(cmpD$table, digits = 3), row.names = FALSE)

    # page 1: components (shaded: plug-in band; dashed: cut posterior, for D_DETAIL)
    # for D_DETAIL also page 2: activity effects without link against the
    # summaries, and page 3: local slopes
    j_pop <- which(LINK == "pop"); j_rho <- which(LINK == "rho")
    fit1D <- NULL
    if (detail) fit1D <- lsdynam_fit(smp$sub, smp$n, D = Dk, stats = st$STATS, link = "none", sigma_a = st$SIGMA_A,
                                     sampler = "mala", n_iter = N_ITER, burn_in = BURN_IN, thin = THIN,
                                     choice = fitD$choice, rate_event = smp$rate_event, seed = 1, verbose = FALSE)
    pdf(paste0(OUT_PREFIX, "_components_D", Dk, ".pdf"), width = 11, height = 5.5)
    op <- par(mfrow = c(1, 2))
    plot_link_level(fitD, component = j_pop, show_cut = detail, xlab = "latent popularity",
                    main = sprintf("Component in latent popularity (D = %d)", Dk), legend_pos = NULL)
    plot_link_level(fitD, component = j_rho, show_cut = detail, xlab = "excess self-affinity",
                    main = sprintf("Component in excess self-affinity (D = %d)", Dk), legend_pos = NULL)
    if (detail) {
      a_M1 <- colMeans(fit1D$rate$a); cexN <- 0.5 + 1.2 * sqrt(fitD$N_out / max(fitD$N_out))
      plot(psD$pop, a_M1, pch = 19, cex = cexN, col = "grey30", xlab = "latent popularity",
           ylab = "activity effect, model without link", main = sprintf("Activity against popularity (D = %d)", Dk))
      abline(v = 0, col = "grey40", lty = 2)
      plot(psD$rho, a_M1, pch = 19, cex = cexN, col = "grey30", xlab = "excess self-affinity",
           ylab = "activity effect, model without link", main = sprintf("Activity against excess self-affinity (D = %d)", Dk))
      abline(v = 0, col = "grey40", lty = 2)
      plot_link_slope(fitD, component = j_pop, xlab = "latent popularity", main = "Local slope, popularity")
      plot_link_slope(fitD, component = j_rho, xlab = "excess self-affinity", main = "Local slope, excess self-affinity")
    }
    par(op); dev.off()

    res_D <- list(fit = fitD, minutes = minutes, ess = dgD$ess, modes = modesD$table,
                  comparison = cmpD$table, cor_with_D1 = cr1, essentials = essentials(fitD, cmpD))
    if (detail) {
      ev_tab <- rbind(whole = c(sd = eD$log_bf10_sd, laplace = eD$log_bf10_laplace),
                      rho_given_pop = c(eD_rho$log_bf10_sd, eD_rho$log_bf10_laplace),
                      pop_given_rho = c(eD_pop$log_bf10_sd, eD_pop$log_bf10_laplace))
      ev_grid <- rbind(whole = link_evidence(fitD, sigma_kappa_grid = SK_GRID)$grid,
                       rho_given_pop = link_evidence(fitD, components = "rho", sigma_kappa_grid = SK_GRID)$grid,
                       pop_given_rho = link_evidence(fitD, components = "pop", sigma_kappa_grid = SK_GRID)$grid)
      contr <- contrast_rows(fitD)
      cat("\nEvidence (log Bayes factors, sigma_a at its value under the link):\n"); print(round(ev_tab, 2))
      cat("For other prior standard deviations of kappa (Laplace):\n"); print(round(ev_grid, 2))
      cat("Contrasts:\n"); print(format(contr, digits = 3), row.names = FALSE)
      cat(sprintf("Residual sd of the activity effects: M1 %.3f, M2 %.3f ; R2_a %.3f\n",
                  fit1D$rate$sigma_a, fitD$rate$sigma_a, link_r2(fitD)))
      res_D$detail <- list(evidence = ev_tab, evidence_grid = ev_grid, contrasts = contr,
                           sigma_a = c(M1 = fit1D$rate$sigma_a, M2 = fitD$rate$sigma_a), R2 = link_r2(fitD),
                           cor_pop_rho = cor(psD$pop, psD$rho), cor_pop_align = cor(psD$pop, psD$align),
                           align_median = median(psD$align),
                           summaries_sd = sapply(psD[c("pop", "rho", "r")], sd), fit1 = fit1D)
    }
    out$F2[[paste0("D", Dk)]] <- res_D
  }
}

# =============================================================================
# F4. Calibration for D = D_DETAIL: association with popularity only (and,
#     optionally, no association). Laplace version of both stages, one random
#     start plus the spectral start per replicate. The observed values to
#     compare with are those at the best mode (first row of the modes table),
#     where, as in the replicates, the summaries are evaluated at the mode.
# =============================================================================

key_D <- paste0("D", D_DETAIL)
if (RUN_F4 || RUN_F4_NONE) {
  if (is.null(out$F2[[key_D]])) stop("F4 needs the fit for D = D_DETAIL from F2 (RUN_F2 <- TRUE and D_DETAIL in D_LINK)")
  saveRDS(out, paste0(OUT_PREFIX, "_results.rds"))      # interim save: F4 is slow
  fitDD <- out$F2[[key_D]]$fit
  out$F4 <- list(observed_at_mode = out$F2[[key_D]]$modes[1, ])
  if (RUN_F4) {
    cat(sprintf("\n=== F4. Calibration with popularity only, D = %d ===\n", D_DETAIL))
    out$F4$pop_only <- null_link_calibration(fitDD, n_rep = N_NULL, n_starts = 1, keep = "pop", seed = 3001)
    calib_print(out$F4$pop_only, "Popularity only")
  }
  if (RUN_F4_NONE) {
    cat(sprintf("\n=== F4. Calibration with no association, D = %d ===\n", D_DETAIL))
    out$F4$none <- null_link_calibration(fitDD, n_rep = N_NULL, n_starts = 1, seed = 4001)
    calib_print(out$F4$none, "No association")
  }
  cat("\nObserved, at the best mode (comparable with the replicates):\n")
  print(format(out$F4$observed_at_mode, digits = 3), row.names = FALSE)
  pdf(paste0(OUT_PREFIX, "_calibration_D", D_DETAIL, ".pdf"), width = 11, height = 5.5)
  ch_mode <- fitDD$choice; ch_mode$U_mean <- ch_mode$map$U; ch_mode$V_mean <- ch_mode$map$V
  fit_mode <- fitDD
  fit_mode$rate <- lsdynam_rate_stage(NULL, smp$n, sl = NULL, choice = ch_mode, link = LINK, df = DF,
                                      winsorize_q = WINSORIZE_Q, priors = fitDD$priors, sigma_a = st$SIGMA_A,
                                      inference = "map", cr = fitDD$rate$cr, seed = 1)
  for (nc in Filter(Negate(is.null), list(out$F4$pop_only, out$F4$none))) {
    op <- par(mfrow = c(1, length(LINK)))
    for (j in seq_along(LINK)) {
      obs <- link_curve(fit_mode, j, grid = nc$grid[[j]])
      matplot(nc$grid[[j]], t(nc$levels[[j]]), type = "l", lty = 1, col = "grey70",
              ylim = range(nc$levels[[j]], obs$level_mean - obs$mu_a_mean, na.rm = TRUE), xlab = LINK[j],
              ylab = "component, relative to mu_a",
              main = sprintf("%s, D = %d: observed at the mode against replicates (%s)", LINK[j], D_DETAIL,
                             if (is.null(nc$keep)) "no association" else "popularity only"), cex.main = 0.85)
      lines(nc$grid[[j]], obs$level_mean - obs$mu_a_mean, lwd = 3, col = "black")
      abline(v = 0, col = "grey40", lty = 2)
    }
    par(op)
  }
  dev.off()
}

# =============================================================================
# F3. Essentials across dimensions, and details for D = D_DETAIL
# =============================================================================

cmp1 <- link_comparison_table(saved$fit2, specs)
ess_tab <- cbind(D1 = essentials(saved$fit2, cmp1),
                 do.call(cbind, lapply(out$F2, function(x) x$essentials)))
if (!is.null(out$F1)) {
  pe <- out$F1$per_event; gain <- c(NA, diff(pe) * (if (!is.null(saved$holdout)) saved$holdout$scores$M1$n_test else NA))
  names(gain) <- names(pe)
  extra <- rbind(choice_score_per_event = pe[colnames(ess_tab)],
                 choice_gain_over_previous_D = gain[colnames(ess_tab)])
  ess_tab <- rbind(ess_tab, extra)
}
cat("\n=== F3. Essentials for D = 1, 2, 3 (additive link; log ML relative to M1 with sigma_a set in each model) ===\n")
print(round(ess_tab, 3))
if (!is.null(out$F1)) { cat("\nHeld-out choice score, differences between successive dimensions:\n"); print(out$F1$differences, digits = 3, row.names = FALSE) }
print_detail <- function() {
  dd <- out$F2[[key_D]]
  if (is.null(dd$detail)) return(invisible(NULL))
  cat(sprintf("\n================ D = %d IN DETAIL ================\n", D_DETAIL))
  cat(sprintf("Fit: %.1f min ; acceptance choice %.3f, rate %.3f\n", dd$minutes, dd$fit$choice$accept_rate, dd$fit$rate$accept_rate))
  cat("ESS (min, median):\n"); print(round(t(sapply(dd$ess, function(e) c(min = min(e, na.rm = TRUE), median = median(e, na.rm = TRUE)))), 0))
  cat("Modes of the choice module:\n"); print(format(dd$modes, digits = 3), row.names = FALSE)
  cat(sprintf("zeta: %.3f [%.3f, %.3f]\n", mean(dd$fit$choice$zeta[, 1]), quantile(dd$fit$choice$zeta[, 1], 0.025),
              quantile(dd$fit$choice$zeta[, 1], 0.975)))
  cat("sigma_a:\n"); print(round(dd$detail$sigma_a, 3)); cat(sprintf("R2_a: %.3f\n", dd$detail$R2))
  cat("Evidence:\n"); print(round(dd$detail$evidence, 2)); print(round(dd$detail$evidence_grid, 2))
  cat("Contrasts:\n"); print(format(dd$detail$contrasts, digits = 3), row.names = FALSE)
  cat("Links on the same positions:\n"); print(format(dd$comparison, digits = 3), row.names = FALSE)
  cat("Correlation of the summaries with D = 1:\n"); print(round(dd$cor_with_D1, 3))
  cat(sprintf("cor(pop, rho) %.3f ; cor(pop, align) %.3f ; median alignment %.2f\n",
              dd$detail$cor_pop_rho, dd$detail$cor_pop_align, dd$detail$align_median))
  cat("sd of the summaries:\n"); print(round(dd$detail$summaries_sd, 3))
  if (!is.null(out$F4$pop_only)) calib_print(out$F4$pop_only, "Calibration, popularity only")
  if (!is.null(out$F4$none)) calib_print(out$F4$none, "Calibration, no association")
  if (!is.null(out$F4)) { cat("Observed at the best mode:\n"); print(format(out$F4$observed_at_mode, digits = 3), row.names = FALSE) }
}
print_detail()
capture.output({
  print(round(ess_tab, 3))
  if (!is.null(out$F1)) print(out$F1$differences, digits = 3, row.names = FALSE)
  print_detail()
}, file = paste0(OUT_PREFIX, "_essentials.txt"))
out$essentials <- ess_tab

saveRDS(out, paste0(OUT_PREFIX, "_results.rds"))
cat("\nDone.\n")
