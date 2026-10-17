# =============================================================================
# simulation_study.R
#
# Simulation study for the latent space DyNAM with position-rate link, with
# the additive link in latent popularity (pop) and excess self-affinity (rho)
# as the fitted model, as in the paper.
#
#   Study A  Recovery and link estimation, for true links in
#              none : no association
#              pop  : association with latent popularity only
#              rho  : association with excess self-affinity only
#              both : association with both
#            - coefficients (xi, zeta), activity effects, the summaries
#              pop, rho and r
#            - each component of the link: error and coverage of the 95%
#              band under the plug-in and the cut posterior, 10-90 contrast
#              (bias, attenuation, coverage), also between the quantiles of
#              the estimated summaries as in the application ("_q"), and the
#              error of the component of each actor
#            - evidence: whole link, each component given the other (in
#              particular rho given pop, whose false-positive rate under
#              'pop' is the main check of the decomposition), and the link in
#              the self-affinity r alone on the same positions
#            - association between the estimated summaries and the number of
#              events sent (shrinkage pulls rho towards -pop for actors with
#              few sent events)
#   Study B  Predictive comparison on held-out events
#            - choice score for D = 0 (M0), 1, 2, 3
#            - rate score with (M2, additive link) and without (M1) the link,
#              overall and for the actors with the fewest training events
#
# The true sender positions are shifted by mean_u in the first coordinate,
# so that the senders share a direction of preference and latent popularity
# has a variance comparable to that of excess self-affinity (as in the
# application, where most sender positions have the same sign). Every
# replicate draws new positions and activity effects from the truth of its
# scenario, simulates the events and applies the two-stage procedure.
#
# Requires the package ls.dynam. Results are written to OUT_DIR as .rds,
# .csv and .pdf.
# =============================================================================

library(ls.dynam)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

if (!exists("MODE"))    MODE    <- "quick"   # "quick": a few minutes, checks that everything runs
                                             # "full" : the study for the paper
if (!exists("CORES"))   CORES   <- 1         # > 1 uses parallel::mclapply (not on Windows)
if (!exists("OUT_DIR")) OUT_DIR <- "sim_results"
if (!exists("RUN_A"))   RUN_A   <- TRUE
if (!exists("RUN_B"))   RUN_B   <- TRUE

cfg <- if (MODE == "quick") {
  # inference = "map" (Laplace approximation in both stages) is used here for
  # speed only. It plugs in the summaries at the MODE of the positions, which
  # are noisier than at the posterior means and attenuate the link when the
  # data are sparse; the full study uses posterior simulation.
  list(n = 30, M = c(1200), D = c(1), truths = c("none", "both"), n_rep = 2,
       inference = "map", n_iter = 1500, burn_in = 500, thin = 2, n_cut = 60, n_starts = 1,
       D_fit = 0:2, n_rep_B = 2, scen_B = data.frame(D = 1, M = 1200, truth = "both"))
} else {
  list(n = 60, M = c(2000, 8000), D = c(1, 2), truths = c("none", "pop", "rho", "both"), n_rep = 100,
       inference = "mcmc", n_iter = 5000, burn_in = 1000, thin = 4, n_cut = 200, n_starts = 2,
       D_fit = 0:3, n_rep_B = 50,
       scen_B = expand.grid(D = c(1, 2), M = c(2000, 8000), truth = c("none", "both"),
                            stringsAsFactors = FALSE))
}
if (exists("CFG_OVERRIDE")) cfg[names(CFG_OVERRIDE)] <- CFG_OVERRIDE   # e.g. CFG_OVERRIDE <- list(n_rep = 20)
cfg$sigma_a <- 0.5                         # residual sd of the activity effects
cfg$xi      <- c(rate_indeg = 0.1)
cfg$zeta    <- c(choice_recip = 0.5)
cfg$mean_u  <- 1                           # shift of the first coordinate of the sender positions
cfg$df      <- 4                           # basis functions per component of the fitted link
cfg$train_frac <- 0.8                      # Study B
cfg$grid_fix <- seq(-1.5, 1.5, length.out = 31)   # raw values of pop and rho, for the figures

LINK <- c("pop", "rho")
f_thr <- function(x) 1.5 * (plogis(3 * x) - 0.5)
# True components, as functions of the centred summaries. Their average over
# the actors is removed in make_truth(), so mu_a = 0 is the average level.
TRUTHS <- list(
  none = list(pop = NULL,  rho = NULL),
  pop  = list(pop = f_thr, rho = NULL),
  rho  = list(pop = NULL,  rho = f_thr),
  both = list(pop = f_thr, rho = f_thr)
)

truth_for <- function(sc, seed) {
  make_truth(n = cfg$n, D = sc$D, link_fun = TRUTHS[[sc$truth]], link_var = LINK, mu_a = 0,
             sigma_a = cfg$sigma_a, xi = cfg$xi, zeta = cfg$zeta,
             mean_u = c(cfg$mean_u, rep(0, sc$D - 1)), seed = seed)
}

dir.create(OUT_DIR, showWarnings = FALSE)

run_jobs <- function(jobs, FUN, cores) {
  if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(jobs, FUN, mc.cores = cores, mc.preschedule = FALSE)
  } else lapply(jobs, FUN)
}

covers <- function(draws, truth) {
  q <- quantile(draws, c(0.025, 0.975), names = FALSE)
  truth >= q[1] && truth <= q[2]
}

# -----------------------------------------------------------------------------
# Study A: one replicate
# -----------------------------------------------------------------------------

component_metrics <- function(fit, tr, j) {
  v <- LINK[j]
  x_true <- tr$x[, v]
  grid  <- quantile(x_true, seq(0.05, 0.95, by = 0.05), names = FALSE)
  truth <- tr$level_fun(grid, j)
  lp <- link_curve(fit, j, grid = grid, source = "plugin")
  lc <- link_curve(fit, j, grid = grid, source = "cut")
  xq <- quantile(x_true, c(0.10, 0.90), names = FALSE)
  cp <- link_contrast(fit, xq[1], xq[2], component = j, source = "plugin")
  cc <- link_contrast(fit, xq[1], xq[2], component = j, source = "cut")
  c_true <- diff(tr$level_fun(xq, j))
  # contrast as reported in the application: between the 10% and 90% quantiles
  # of the ESTIMATED (plug-in) summary, which does not depend on the scale of
  # the estimated summaries, against the true contrast between the same
  # quantiles of the true summary
  cpq <- link_contrast(fit, component = j, source = "plugin")
  ccq <- link_contrast(fit, component = j, source = "cut")
  # component of each actor, centred over the actors: estimated at the plug-in
  # summary against the truth at the true summary
  g_hat <- link_curve(fit, j, grid = fit$rate$link$X[, j], source = "plugin")$level_mean
  g_tru <- tr$level_fun(x_true, j)
  g_hat <- g_hat - mean(g_hat); g_tru <- g_tru - mean(g_tru)
  e <- link_evidence(fit, components = j)
  out <- c(level_rmse = sqrt(mean((lp$level_mean - truth)^2)),
           cover_plugin = mean(truth >= lp$level_ci[, 1] & truth <= lp$level_ci[, 2]),
           cover_cut    = mean(truth >= lc$level_ci[, 1] & truth <= lc$level_ci[, 2]),
           width_plugin = mean(lp$level_ci[, 2] - lp$level_ci[, 1]),
           width_cut    = mean(lc$level_ci[, 2] - lc$level_ci[, 1]),
           contrast_true = c_true, contrast_plugin = cp$mean, contrast_cut = cc$mean,
           contrast_cover_plugin = as.numeric(c_true >= cp$ci[1] && c_true <= cp$ci[2]),
           contrast_cover_cut    = as.numeric(c_true >= cc$ci[1] && c_true <= cc$ci[2]),
           excl0_plugin = as.numeric(cp$ci[1] > 0 || cp$ci[2] < 0),
           excl0_cut    = as.numeric(cc$ci[1] > 0 || cc$ci[2] < 0),
           contrast_q_plugin = cpq$mean, contrast_q_cut = ccq$mean,
           contrast_q_cover_plugin = as.numeric(c_true >= cpq$ci[1] && c_true <= cpq$ci[2]),
           contrast_q_cover_cut    = as.numeric(c_true >= ccq$ci[1] && c_true <= ccq$ci[2]),
           excl0_q_plugin = as.numeric(cpq$ci[1] > 0 || cpq$ci[2] < 0),
           excl0_q_cut    = as.numeric(ccq$ci[1] > 0 || ccq$ci[2] < 0),
           actor_rmse = sqrt(mean((g_hat - g_tru)^2)),
           actor_cor = if (sd(g_tru) > 0) cor(g_hat, g_tru) else NA_real_,
           logbf_given_rest_sd = e$log_bf10_sd, logbf_given_rest = e$log_bf10_laplace)
  names(out) <- paste0(v, "_", names(out))
  out
}

run_rep_A <- function(job) {
  sc <- job$sc; rep_id <- job$rep
  seed <- 100000 * sc$id + rep_id
  out <- tryCatch({
    tr <- truth_for(sc, seed)
    ev <- simulate_lsdynam(tr$theta, cfg$n, M = sc$M, seed = seed + 1)
    N  <- tabulate(ev$sender, nbins = cfg$n)
    t0 <- Sys.time()
    fit <- lsdynam_fit(ev, cfg$n, D = sc$D, stats = tr$stats, link = LINK, df = cfg$df, sigma_a = "eb",
                       inference = cfg$inference, n_iter = cfg$n_iter, burn_in = cfg$burn_in,
                       thin = cfg$thin, n_cut = cfg$n_cut, n_starts = cfg$n_starts,
                       seed = seed + 2, verbose = FALSE)
    secs <- as.numeric(Sys.time() - t0, units = "secs")

    a_hat <- colMeans(fit$rate$a)
    a_cov <- mean(vapply(seq_len(cfg$n), function(i) covers(fit$rate$a[, i], tr$theta$a[i]), logical(1)))
    comp <- c(component_metrics(fit, tr, 1), component_metrics(fit, tr, 2))
    evd <- link_evidence(fit)

    # the link in the self-affinity alone, on the same positions (Laplace version of stage 2)
    fr <- lsdynam_fit(ev, cfg$n, D = sc$D, stats = tr$stats, link = "r", df = cfg$df, sigma_a = "eb",
                      inference = "map", choice = fit$choice, seed = seed + 3, verbose = FALSE)
    cr_ <- link_contrast(fr)

    ps <- fit$summaries; pt <- tr$summaries
    sp <- function(x) suppressWarnings(cor(x, N, method = "spearman"))
    row <- data.frame(
      scenario = sc$id, D = sc$D, M = sc$M, truth = sc$truth, rep = rep_id,
      N_min = min(N), N_median = median(N),
      zeta_hat = mean(fit$choice$zeta[, 1]), zeta_true = unname(cfg$zeta[1]),
      zeta_cover = covers(fit$choice$zeta[, 1], cfg$zeta[1]),
      xi_hat = mean(fit$rate$xi[, 1]), xi_true = unname(cfg$xi[1]),
      xi_cover = covers(fit$rate$xi[, 1], cfg$xi[1]),
      sigma_a_hat = fit$rate$sigma_a, mu_a_hat = mean(fit$rate$mu_a),
      a_rmse = sqrt(mean((a_hat - tr$theta$a)^2)), a_cover = a_cov,
      pop_cor = cor(ps$pop, pt$pop), rho_cor = cor(ps$rho, pt$rho), r_cor = cor(ps$r, pt$r),
      pop_slope = unname(coef(lm(ps$pop ~ pt$pop))[2]), rho_slope = unname(coef(lm(ps$rho ~ pt$rho))[2]),
      cor_pop_rho_true = cor(pt$pop, pt$rho),
      t(comp),
      logbf_sd = evd$log_bf10_sd, logbf_laplace = evd$log_bf10_laplace, r2 = link_r2(fit),
      logbf_r = link_evidence(fr)$log_bf10_laplace, excl0_r = as.numeric(cr_$ci[1] > 0 || cr_$ci[2] < 0),
      r2_r = link_r2(fr),
      rho_pop_N = sp(ps$pop), rho_absrho_N = sp(abs(ps$rho)), rho_absr_N = sp(abs(ps$r)),
      accept_choice = fit$choice$accept_rate, accept_rate = fit$rate$accept_rate,
      seconds = secs, error = NA_character_, stringsAsFactors = FALSE, row.names = NULL)
    pq <- seq(0.05, 0.95, by = 0.05)
    curves <- lapply(seq_along(LINK), function(j) {
      lf <- suppressWarnings(link_curve(fit, j, grid = cfg$grid_fix, source = "plugin"))
      # component against the quantiles of the summary (estimated at the
      # quantiles of the plug-in summary, truth at those of the true summary),
      # relative to the average level, which does not depend on the scale of
      # the estimated summaries
      xq_hat <- quantile(fit$rate$link$X[, j], pq, names = FALSE)
      lq  <- link_curve(fit, j, grid = xq_hat, source = "plugin")
      lqc <- link_curve(fit, j, grid = xq_hat, source = "cut")
      list(est = lf$level_mean, truth = tr$level_fun(cfg$grid_fix, j),
           est_q = lq$level_mean - lq$mu_a_mean, est_q_cut = lqc$level_mean - lqc$mu_a_mean,
           truth_q = tr$level_fun(quantile(tr$x[, LINK[j]], pq, names = FALSE), j) - tr$mu_a)
    })
    names(curves) <- LINK
    list(row = row, curve = curves)
  }, error = function(e) {
    list(row = data.frame(scenario = sc$id, D = sc$D, M = sc$M, truth = sc$truth, rep = rep_id,
                          error = conditionMessage(e), stringsAsFactors = FALSE), curve = NULL)
  })
  out
}

# -----------------------------------------------------------------------------
# Study B: one replicate
# -----------------------------------------------------------------------------

run_rep_B <- function(job) {
  sc <- job$sc; rep_id <- job$rep
  seed <- 500000 + 100000 * sc$id + rep_id
  tryCatch({
    tr <- truth_for(sc, seed)
    ev <- simulate_lsdynam(tr$theta, cfg$n, M = sc$M, seed = seed + 1)
    m_tr <- round(cfg$train_frac * nrow(ev))
    ev_tr <- ev[seq_len(m_tr), ]
    N_tr <- tabulate(ev_tr$sender, nbins = cfg$n)
    low <- N_tr <= quantile(N_tr, 0.25)                 # actors with the fewest training events

    res <- list()
    for (Df in cfg$D_fit) {
      f1 <- lsdynam_fit(ev_tr, cfg$n, D = Df, stats = tr$stats, link = "none", sigma_a = "eb",
                        inference = "map", n_starts = cfg$n_starts, seed = seed + 2, verbose = FALSE)
      s1 <- predictive_scores(f1, ev, cfg$n, m_tr, n_draws = 100)
      res[[paste0("choice_D", Df)]] <- s1$choice
      if (Df == 0) {                                     # the rate module of M0 and M1 is the same
        res$rate_M1 <- s1$rate
        res$rate_M1_low <- sum(s1$rate_by_actor[low]); res$rate_M1_high <- sum(s1$rate_by_actor[!low])
      }
      if (Df == sc$D) {
        f2 <- lsdynam_fit(ev_tr, cfg$n, D = Df, stats = tr$stats, link = LINK, df = cfg$df, sigma_a = "eb",
                          inference = "map", choice = f1$choice, seed = seed + 2, verbose = FALSE)
        s2 <- predictive_scores(f2, ev, cfg$n, m_tr, n_draws = 100)
        res$rate_M2 <- s2$rate
        res$rate_M2_low <- sum(s2$rate_by_actor[low]); res$rate_M2_high <- sum(s2$rate_by_actor[!low])
      }
    }
    ch <- unlist(res[grep("^choice_D", names(res))])
    data.frame(scenario = sc$id, D = sc$D, M = sc$M, truth = sc$truth, rep = rep_id,
               as.data.frame(res), D_best = cfg$D_fit[which.max(ch)], n_test = nrow(ev) - m_tr,
               n_low = sum(low), N_low_median = median(N_tr[low]), error = NA_character_,
               stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(scenario = sc$id, D = sc$D, M = sc$M, truth = sc$truth, rep = rep_id,
               error = conditionMessage(e), stringsAsFactors = FALSE)
  })
}

rbind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}

# -----------------------------------------------------------------------------
# Run Study A
# -----------------------------------------------------------------------------

if (RUN_A) {
  scen_A <- expand.grid(D = cfg$D, M = cfg$M, truth = cfg$truths, stringsAsFactors = FALSE)
  scen_A$id <- seq_len(nrow(scen_A))
  jobs <- unlist(lapply(seq_len(nrow(scen_A)), function(s)
    lapply(seq_len(cfg$n_rep), function(r) list(sc = scen_A[s, ], rep = r))), recursive = FALSE)
  cat(sprintf("Study A: %d scenarios x %d replicates, n = %d, inference = %s\n",
              nrow(scen_A), cfg$n_rep, cfg$n, cfg$inference))
  t0 <- Sys.time()
  out_A <- run_jobs(jobs, run_rep_A, CORES)
  cat(sprintf("Study A took %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
  res_A <- rbind_fill(lapply(out_A, function(o) o$row))
  curves_A <- lapply(out_A, function(o) o$curve)
  saveRDS(list(results = res_A, curves = curves_A, scenarios = scen_A, cfg = cfg), file.path(OUT_DIR, "study_A.rds"))
  write.csv(res_A, file.path(OUT_DIR, "study_A_replicates.csv"), row.names = FALSE)

  n_err <- sum(!is.na(res_A$error))
  if (n_err > 0) { cat(n_err, "replicate(s) failed:\n"); print(table(res_A$error)) }
  ok <- res_A[is.na(res_A$error), ]
  by <- list(D = ok$D, M = ok$M, truth = ok$truth)
  agg <- function(vars, FUN = mean) aggregate(ok[vars], by, FUN, na.rm = TRUE)

  cat("\n=== A1. Coefficients, activity effects, summaries (means over replicates) ===\n")
  ok$zeta_bias <- ok$zeta_hat - ok$zeta_true; ok$xi_bias <- ok$xi_hat - ok$xi_true
  tab_A1 <- agg(c("N_median", "zeta_bias", "zeta_cover", "xi_bias", "xi_cover", "sigma_a_hat",
                  "a_rmse", "a_cover", "pop_cor", "rho_cor", "r_cor", "pop_slope", "rho_slope", "cor_pop_rho_true"))
  print(format(tab_A1, digits = 3), row.names = FALSE)

  cat("\n=== A2. Components of the link: error, coverage of the 95% band, 10-90 contrast ===\n")
  tab_A2 <- do.call(rbind, lapply(LINK, function(v) {
    t <- agg(paste0(v, "_", c("level_rmse", "cover_plugin", "cover_cut", "width_plugin", "width_cut",
                              "contrast_true", "contrast_plugin", "contrast_cut",
                              "contrast_cover_plugin", "contrast_cover_cut",
                              "contrast_q_plugin", "contrast_q_cut", "contrast_q_cover_plugin",
                              "contrast_q_cover_cut", "actor_rmse", "actor_cor")))
    names(t) <- sub(paste0("^", v, "_"), "", names(t))
    cbind(component = v, t)
  }))
  print(format(tab_A2, digits = 3), row.names = FALSE)

  cat("\n=== A3. Evidence: proportion of replicates in which an association is declared ===\n")
  cat("(BF10 > 3 by Laplace evidence; 'excl0' = 95% interval of the 10-90 contrast excludes zero.\n")
  cat(" Under truth = none every column is a false-positive rate; under truth = pop the columns for rho\n")
  cat(" are the false-positive rates of the decomposition, under truth = rho those for pop.)\n")
  ok$bf3_link <- ok$logbf_laplace > log(3)
  ok$bf3_pop_given_rho <- ok$pop_logbf_given_rest > log(3)
  ok$bf3_rho_given_pop <- ok$rho_logbf_given_rest > log(3)
  ok$bf3_r <- ok$logbf_r > log(3)
  tab_A3 <- agg(c("bf3_link", "bf3_pop_given_rho", "bf3_rho_given_pop", "pop_excl0_plugin", "rho_excl0_plugin",
                  "pop_excl0_cut", "rho_excl0_cut", "r2", "bf3_r", "excl0_r", "r2_r"))
  print(format(tab_A3, digits = 3), row.names = FALSE)
  cat("\nMedian log Bayes factors:\n")
  tab_A3m <- agg(c("logbf_laplace", "pop_logbf_given_rest", "rho_logbf_given_rest", "logbf_r"), FUN = median)
  print(format(tab_A3m, digits = 3), row.names = FALSE)

  cat("\n=== A4. Rank correlation of the estimated summaries with the number of events sent ===\n")
  tab_A4 <- agg(c("rho_pop_N", "rho_absrho_N", "rho_absr_N"))
  print(format(tab_A4, digits = 3), row.names = FALSE)

  cat(sprintf("\nAgreement of the two Bayes factors for the whole link (Savage-Dickey vs Laplace): correlation %.3f, mean abs difference %.2f\n",
              cor(ok$logbf_sd, ok$logbf_laplace), mean(abs(ok$logbf_sd - ok$logbf_laplace))))
  cat(sprintf("Seconds per replicate: median %.0f\n", median(ok$seconds)))
  for (nm in c("A1", "A2", "A3", "A3m", "A4")) {
    write.csv(get(paste0("tab_", nm)), file.path(OUT_DIR, paste0("study_A_table_", nm, ".csv")), row.names = FALSE)
  }

  # Figure: estimated components across replicates against the truth
  dm <- unique(scen_A[c("D", "M")])
  pdf(file.path(OUT_DIR, "study_A_components.pdf"), width = 3.4 * length(cfg$truths), height = 3.2 * 2 * nrow(dm))
  op <- par(mfrow = c(2 * nrow(dm), length(cfg$truths)), mar = c(4, 4, 2.5, 1))
  for (k in seq_len(nrow(dm))) for (v in LINK) for (tt in cfg$truths) {
    sid <- scen_A$id[scen_A$D == dm$D[k] & scen_A$M == dm$M[k] & scen_A$truth == tt]
    idx <- which(res_A$scenario == sid & is.na(res_A$error))
    est <- do.call(rbind, lapply(curves_A[idx], function(cv) cv[[v]]$est))
    tru <- do.call(rbind, lapply(curves_A[idx], function(cv) cv[[v]]$truth))
    if (is.null(est)) { plot.new(); next }
    q <- apply(est, 2, quantile, probs = c(0.05, 0.95))
    plot(cfg$grid_fix, colMeans(est), type = "n", ylim = range(q, tru), xlab = v,
         ylab = "component of the link", main = sprintf("%s | truth %s, D = %d, M = %d", v, tt, dm$D[k], dm$M[k]),
         cex.main = 0.85)
    polygon(c(cfg$grid_fix, rev(cfg$grid_fix)), c(q[1, ], rev(q[2, ])), col = "grey80", border = NA)
    lines(cfg$grid_fix, colMeans(est), lwd = 2, col = "black")
    lines(cfg$grid_fix, colMeans(tru), lwd = 2, lty = 4, col = "grey40")
    abline(v = 0, col = "grey50", lty = 2)
    if (k == 1 && v == LINK[1] && tt == cfg$truths[1]) legend("topleft", c("mean estimate, 5-95% of replicates", "truth"),
                                                              col = c("black", "grey40"), lty = c(1, 4), lwd = 2, bty = "n", cex = 0.7)
  }
  par(op); dev.off()
}

# -----------------------------------------------------------------------------
# Run Study B
# -----------------------------------------------------------------------------

if (RUN_B) {
  scen_B <- cfg$scen_B
  scen_B$id <- seq_len(nrow(scen_B))
  jobs <- unlist(lapply(seq_len(nrow(scen_B)), function(s)
    lapply(seq_len(cfg$n_rep_B), function(r) list(sc = scen_B[s, ], rep = r))), recursive = FALSE)
  cat(sprintf("\nStudy B: %d scenarios x %d replicates, D fitted = %s\n",
              nrow(scen_B), cfg$n_rep_B, paste(cfg$D_fit, collapse = ", ")))
  t0 <- Sys.time()
  res_B <- rbind_fill(run_jobs(jobs, run_rep_B, CORES))
  cat(sprintf("Study B took %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
  saveRDS(list(results = res_B, scenarios = scen_B, cfg = cfg), file.path(OUT_DIR, "study_B.rds"))
  write.csv(res_B, file.path(OUT_DIR, "study_B_replicates.csv"), row.names = FALSE)

  n_err <- sum(!is.na(res_B$error))
  if (n_err > 0) { cat(n_err, "replicate(s) failed:\n"); print(table(res_B$error)) }
  okB <- res_B[is.na(res_B$error), ]
  byB <- list(D_true = okB$D, M = okB$M, truth = okB$truth)

  cat("\n=== B1. Held-out choice score per event, relative to D = 0 (M0); share of replicates choosing each D ===\n")
  ch_cols <- paste0("choice_D", cfg$D_fit)
  gain <- (okB[ch_cols] - okB$choice_D0) / okB$n_test
  names(gain) <- paste0("gain_D", cfg$D_fit)
  tab_B1 <- aggregate(gain, byB, mean)
  pk <- as.data.frame(lapply(cfg$D_fit, function(d) okB$D_best == d))
  names(pk) <- paste0("picked_D", cfg$D_fit)
  pick <- aggregate(pk, byB, mean)
  tab_B1 <- merge(tab_B1, pick)
  print(format(tab_B1, digits = 3), row.names = FALSE)

  cat("\n=== B2. Held-out rate score: model with the additive link (M2) minus model without (M1) ===\n")
  cat("(sums over actors of the actor-level scores; 'low' = quarter of the actors with the fewest training events)\n")
  d <- data.frame(diff_total = okB$rate_M2 - okB$rate_M1,          # = diff_low + diff_high
                  diff_low = okB$rate_M2_low - okB$rate_M1_low,
                  diff_high = okB$rate_M2_high - okB$rate_M1_high,
                  M2_better = okB$rate_M2 > okB$rate_M1,
                  N_low_median = okB$N_low_median)
  tab_B2 <- aggregate(d, byB, mean)
  print(format(tab_B2, digits = 3), row.names = FALSE)
  write.csv(tab_B1, file.path(OUT_DIR, "study_B_table_B1.csv"), row.names = FALSE)
  write.csv(tab_B2, file.path(OUT_DIR, "study_B_table_B2.csv"), row.names = FALSE)
}

cat("\nDone. Output in", OUT_DIR, "\n")
