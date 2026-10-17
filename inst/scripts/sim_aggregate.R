# Combines the per-replicate files written by sim_driver.R and prints the
# summaries used in the appendix; writes CSV tables and a greyscale figure.
OUT_DIR <- "sim_full"
rbind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
fa <- sort(list.files(file.path(OUT_DIR, "reps"), pattern = "^A_.*\\.rds$", full.names = TRUE))
fb <- sort(list.files(file.path(OUT_DIR, "reps"), pattern = "^B_.*\\.rds$", full.names = TRUE))
cat(length(fa), "replicates of study A,", length(fb), "of study B\n")
options(width = 200)
if (length(fa) > 0) {
  outA <- lapply(fa, readRDS)
  resA <- rbind_fill(lapply(outA, function(o) o$row))
  curvesA <- lapply(outA, function(o) o$curve)
  cat("errors:", sum(!is.na(resA$error)), "\n"); if (any(!is.na(resA$error))) print(table(resA$error))
  ok <- resA[is.na(resA$error), ]
  ok$truth <- factor(ok$truth, levels = c("none", "pop", "rho", "both"))
  by <- list(D = ok$D, M = ok$M, truth = ok$truth)
  agg <- function(vars, FUN = mean) aggregate(ok[vars], by, FUN, na.rm = TRUE)
  cnt <- aggregate(list(reps = ok$rep), by, length)
  L2 <- log(3)
  ok$bf3_link <- ok$logbf_laplace > L2
  ok$bf3_pop_given_rho <- ok$pop_logbf_given_rest > L2
  ok$bf3_rho_given_pop <- ok$rho_logbf_given_rest > L2
  ok$bf3_r <- ok$logbf_r > L2
  ok$pop_att_plugin <- ok$pop_contrast_plugin / ok$pop_contrast_true
  ok$pop_att_cut <- ok$pop_contrast_cut / ok$pop_contrast_true
  ok$rho_att_plugin <- ok$rho_contrast_plugin / ok$rho_contrast_true
  ok$rho_att_cut <- ok$rho_contrast_cut / ok$rho_contrast_true
  for (v in c("pop", "rho")) for (src in c("plugin", "cut"))
    ok[[sprintf("%s_attq_%s", v, src)]] <- ok[[sprintf("%s_contrast_q_%s", v, src)]] / ok[[sprintf("%s_contrast_true", v)]]
  ok$zeta_bias <- ok$zeta_hat - ok$zeta_true; ok$xi_bias <- ok$xi_hat - ok$xi_true
  tab1 <- merge(cnt, agg(c("N_median", "pop_cor", "rho_cor", "r_cor", "cor_pop_rho_true", "zeta_bias", "zeta_cover",
                           "xi_bias", "xi_cover", "sigma_a_hat", "a_rmse", "a_cover")))
  tab2 <- agg(c("pop_contrast_true", "pop_contrast_plugin", "pop_contrast_cut", "pop_contrast_cover_plugin",
                "pop_contrast_cover_cut", "pop_cover_plugin", "pop_cover_cut", "pop_level_rmse",
                "rho_contrast_true", "rho_contrast_plugin", "rho_contrast_cut", "rho_contrast_cover_plugin",
                "rho_contrast_cover_cut", "rho_cover_plugin", "rho_cover_cut", "rho_level_rmse"))
  qv <- function(v) paste0(v, "_", c("contrast_true", "contrast_q_plugin", "contrast_q_cut", "contrast_q_cover_plugin",
                                     "contrast_q_cover_cut", "excl0_q_plugin", "excl0_q_cut", "actor_cor", "actor_rmse"))
  tab2q <- agg(c(qv("pop"), qv("rho")))
  tab3 <- agg(c("bf3_link", "bf3_pop_given_rho", "bf3_rho_given_pop", "pop_excl0_plugin", "rho_excl0_plugin",
                "pop_excl0_cut", "rho_excl0_cut", "r2", "bf3_r", "r2_r"))
  tab3m <- agg(c("logbf_laplace", "pop_logbf_given_rest", "rho_logbf_given_rest", "logbf_r"), FUN = median)
  tab4 <- agg(c("rho_pop_N", "rho_absrho_N", "rho_absr_N", "accept_choice", "accept_rate", "seconds"))
  for (nm in c("tab1", "tab2", "tab2q", "tab3", "tab3m", "tab4")) {
    cat("\n=== ", nm, " ===\n"); print(format(get(nm), digits = 3), row.names = FALSE)
    write.csv(get(nm), file.path(OUT_DIR, paste0("A_", nm, ".csv")), row.names = FALSE)
  }
  cat(sprintf("\nSavage-Dickey vs Laplace (whole link): cor %.3f, mean abs diff %.2f\n",
              cor(ok$logbf_sd, ok$logbf_laplace), mean(abs(ok$logbf_sd - ok$logbf_laplace))))
  # attenuation where the component is active
  act <- ok
  cat("\n=== ratio estimated / true contrast, active components (mean over replicates) ===\n")
  print(format(aggregate(act[act$truth %in% c("pop", "both"), c("pop_att_plugin", "pop_att_cut", "pop_attq_plugin", "pop_attq_cut")],
                         list(D = act$D[act$truth %in% c("pop", "both")], M = act$M[act$truth %in% c("pop", "both")]), mean), digits = 3))
  print(format(aggregate(act[act$truth %in% c("rho", "both"), c("rho_att_plugin", "rho_att_cut", "rho_attq_plugin", "rho_attq_cut")],
                         list(D = act$D[act$truth %in% c("rho", "both")], M = act$M[act$truth %in% c("rho", "both")]), mean), digits = 3))
  saveRDS(list(results = resA, curves = curvesA), file.path(OUT_DIR, "study_A_combined.rds"))

  # Figure: components against the percentile of the summary, across
  # replicates, for one (D, M): truth (grey, dashed), mean of the plug-in
  # estimates (black) with 5-95% of the replicates (shaded), mean of the
  # cut-posterior estimates (black, dotted)
  pq <- 100 * seq(0.05, 0.95, by = 0.05)
  fig <- function(Dsel, Msel, file) {
    pdf(file, width = 9, height = 4.8)
    op <- par(mfrow = c(2, 4), mar = c(3.6, 3.8, 2, 0.6), mgp = c(2.2, 0.7, 0), cex.lab = 1.05)
    for (v in c("pop", "rho")) for (tt in c("none", "pop", "rho", "both")) {
      idx <- which(resA$D == Dsel & resA$M == Msel & resA$truth == tt & is.na(resA$error))
      idx <- idx[vapply(curvesA[idx], function(cv) !is.null(cv[[v]]$est_q), logical(1))]
      if (length(idx) == 0) { plot.new(); next }
      est <- do.call(rbind, lapply(curvesA[idx], function(cv) cv[[v]]$est_q))
      estc <- do.call(rbind, lapply(curvesA[idx], function(cv) cv[[v]]$est_q_cut))
      tru <- do.call(rbind, lapply(curvesA[idx], function(cv) cv[[v]]$truth_q))
      q <- apply(est, 2, quantile, probs = c(0.05, 0.95), na.rm = TRUE)
      xl <- if (v == "pop") "percentile of latent popularity" else "percentile of excess self-affinity"
      ttl <- c(none = "No association", pop = "Popularity only", rho = "Excess self-affinity only", both = "Both")[tt]
      plot(pq, colMeans(est), type = "n", ylim = c(-1.3, 1.3), xlab = xl,
           ylab = if (tt == "none") (if (v == "pop") expression("component in " * pi) else expression("component in " * rho)) else "",
           main = if (v == "pop") ttl else "", cex.main = 1, font.main = 1, xaxt = "n")
      axis(1, at = c(5, 25, 50, 75, 95))
      polygon(c(pq, rev(pq)), c(q[1, ], rev(q[2, ])), col = "grey85", border = NA)
      abline(h = 0, col = "grey60", lty = 3)
      lines(pq, colMeans(tru), lwd = 2.2, lty = 2, col = "grey45")
      lines(pq, colMeans(estc), lwd = 1.5, lty = 3, col = "black")
      lines(pq, colMeans(est), lwd = 1.8, col = "black")
    }
    par(op); dev.off()
  }
  for (Dsel in unique(resA$D)) for (Msel in unique(resA$M))
    fig(Dsel, Msel, file.path(OUT_DIR, sprintf("sim_components_D%d_M%d.pdf", Dsel, Msel)))
}
if (length(fb) > 0) {
  resB <- rbind_fill(lapply(fb, readRDS))
  cat("\nStudy B errors:", sum(!is.na(resB$error)), "\n")
  okB <- resB[is.na(resB$error), ]
  byB <- list(D_true = okB$D, M = okB$M, truth = okB$truth)
  chc <- grep("^choice_D", names(okB), value = TRUE)
  gain <- (okB[chc] - okB$choice_D0) / okB$n_test
  names(gain) <- sub("choice_", "gain_", chc)
  tb1 <- merge(aggregate(list(reps = okB$rep), byB, length), aggregate(gain, byB, mean))
  pk <- as.data.frame(lapply(as.integer(sub("choice_D", "", chc)), function(d) okB$D_best == d))
  names(pk) <- paste0("picked_D", sub("choice_D", "", chc))
  tb1 <- merge(tb1, aggregate(pk, byB, mean))
  gain_tot <- okB[chc] - okB$choice_D0
  names(gain_tot) <- sub("choice_", "total_gain_", chc)
  tb1 <- merge(tb1, aggregate(gain_tot, byB, mean))
  d <- data.frame(diff_total = okB$rate_M2 - okB$rate_M1, diff_low = okB$rate_M2_low - okB$rate_M1_low,
                  diff_high = okB$rate_M2_high - okB$rate_M1_high, M2_better = okB$rate_M2 > okB$rate_M1,
                  low_better = okB$rate_M2_low > okB$rate_M1_low, N_low_median = okB$N_low_median, n_test = okB$n_test)
  tb2 <- aggregate(d, byB, mean)
  se <- aggregate(d[c("diff_total", "diff_low", "diff_high")], byB, function(x) sd(x) / sqrt(length(x)))
  names(se)[4:6] <- paste0("se_", names(se)[4:6])
  tb2 <- merge(tb2, se)
  cat("\n=== B1 ===\n"); print(format(tb1, digits = 3), row.names = FALSE)
  cat("\n=== B2 ===\n"); print(format(tb2, digits = 3), row.names = FALSE)
  write.csv(tb1, file.path(OUT_DIR, "B_tab1.csv"), row.names = FALSE)
  write.csv(tb2, file.path(OUT_DIR, "B_tab2.csv"), row.names = FALSE)
}
