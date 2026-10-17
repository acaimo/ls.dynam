#' Convergence diagnostics
#'
#' Effective sample sizes and acceptance rates of the two samplers. The
#' invariant summaries of the positions (the self-affinity and the summaries
#' in the link) are monitored instead of individual coordinates, since they do
#' not depend on the alignment of the draws.
#'
#' @param fit an object of class `lsdynam` fitted with `inference = "mcmc"`.
#' @param verbose print the minimum and median effective sample size of each
#'   group of quantities.
#' @return A list with the effective sample sizes `ess` (one vector per group
#'   of quantities), the acceptance rates `accept` and the numbers of draws
#'   kept `n_draws`, invisibly.
#' @seealso [ess()], [plot_traces()].
#' @examples
#' tr <- make_truth(n = 10, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 10, M = 200, seed = 2)
#' fit <- lsdynam_fit(ev, n = 10, link = "none", n_iter = 400, burn_in = 200, thin = 1,
#'                    verbose = FALSE)
#' mcmc_diagnostics(fit)
#' @export
mcmc_diagnostics <- function(fit, verbose = TRUE) {
  ch <- fit$choice; rt <- fit$rate
  e <- list()
  if (ncol(ch$zeta) > 0) e$zeta <- apply(ch$zeta, 2, ess)
  if (fit$D > 0) {
    vars <- unique(c("r", if (!identical(rt$link$vars, "none")) rt$link$vars))
    sdr <- summaries_draws(ch$U, ch$V, vars)
    for (v in vars) e[[v]] <- apply(matrix(sdr[, , v], nrow = dim(sdr)[1]), 2, ess)
  }
  if (ncol(rt$xi) > 0) e$xi <- apply(rt$xi, 2, ess)
  e$mu_a <- ess(rt$mu_a)
  if (ncol(rt$kappa_vec) > 0) e$kappa <- apply(rt$kappa_vec, 2, ess)
  e$a <- apply(rt$a, 2, ess)
  out <- list(ess = e, accept = c(choice = ch$accept_rate, rate = rt$accept_rate),
              n_draws = c(choice = nrow(ch$zeta), rate = nrow(rt$a)))
  if (verbose) {
    cat(sprintf("Draws kept: choice %d, rate %d | acceptance: choice %.3f, rate %.3f\n",
                out$n_draws[1], out$n_draws[2], out$accept[1], out$accept[2]))
    for (nm in names(e)) cat(sprintf("  ESS %-6s min %7.0f  median %7.0f\n", nm, min(e[[nm]], na.rm = TRUE), median(e[[nm]], na.rm = TRUE)))
  }
  invisible(out)
}

#' Shrinkage check for the summaries
#'
#' Is a summary of the positions related to the amount of data on each
#' actor? Shrinkage pulls the sender position of an actor with few events
#' sent towards the origin, and with it the self-affinity towards zero from
#' either side, so the absolute value of a summary is checked as well as its
#' sign: the function reports Spearman correlations of each summary, and of
#' its absolute value, with the number of events sent.
#'
#' @param fit an object of class `lsdynam` with \eqn{D \ge 1}{D >= 1}.
#' @param vars summaries to check (columns of [position_summaries()]).
#' @param plot plot each summary against the number of events sent.
#' @param verbose print the table.
#' @return A data frame of correlations and p-values, invisibly.
#' @examples
#' tr <- make_truth(n = 15, D = 1, mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, link = "none", inference = "map", verbose = FALSE)
#' link_shrinkage_check(fit, vars = c("pop", "rho"), plot = FALSE)
#' @export
link_shrinkage_check <- function(fit, vars = c("r", "z"), plot = TRUE, verbose = TRUE) {
  ps <- fit$summaries; N <- fit$N_out
  res <- do.call(rbind, lapply(vars, function(v) {
    x <- ps[[v]]
    c1 <- suppressWarnings(cor.test(x, N, method = "spearman"))
    c2 <- suppressWarnings(cor.test(abs(x), N, method = "spearman"))
    data.frame(var = v, rho_signed = unname(c1$estimate), p_signed = c1$p.value,
               rho_abs = unname(c2$estimate), p_abs = c2$p.value)
  }))
  if (plot) {
    op <- par(mfrow = c(1, length(vars))); on.exit(par(op))
    for (v in vars) {
      plot(N, ps[[v]], log = "x", pch = 19, col = "grey30",
           xlab = "events sent (log scale)", ylab = v,
           main = sprintf("%s: rho = %.2f, rho(|%s|) = %.2f", v, res$rho_signed[res$var == v], v, res$rho_abs[res$var == v]))
      abline(h = 0, col = "grey50", lty = 2)
    }
  }
  if (verbose) print(res, digits = 3, row.names = FALSE)
  invisible(res)
}

#' Leverage check for the link
#'
#' Refits the rate module with the `n_exclude` actors with the most extreme
#' values of one summary removed from the link (they keep the exchangeable
#' prior of the activity effects) and compares the level of the link over the
#' range of the other actors.
#'
#' @param fit an object of class `lsdynam` with a link.
#' @param events_df not needed (kept for compatibility); the statistics are
#'   taken from `fit`.
#' @param n_exclude number of actors removed from the link.
#' @param component index or name of the component.
#' @param plot plot the two levels.
#' @param seed seed of the random number generator.
#' @param verbose print the actors removed and the largest change.
#' @return A list with the actors removed `excluded`, the two levels
#'   (`full`, `reduced`, see [link_curve()]) and the largest change in the
#'   level `max_diff`, invisibly.
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' link_leverage_check(fit, n_exclude = 2, component = "pop", plot = FALSE)$max_diff
#' @export
link_leverage_check <- function(fit, events_df = NULL, n_exclude = 5, component = 1, plot = TRUE, seed = 1,
                                verbose = TRUE) {
  lk <- fit$rate$link
  x <- lk$X[, component]
  excl <- order(abs(x - mean(x)), decreasing = TRUE)[seq_len(n_exclude)]
  rt2 <- lsdynam_rate_stage(events_df, fit$n, sl = NULL, choice = fit$choice, link = lk$vars,
                            link_type = lk$type, df = lk$df, winsorize_q = lk$winsorize_q,
                            priors = fit$priors, sigma_a = fit$rate$sigma_a, nu = fit$rate$nu,
                            inference = "map", exclude_from_link = excl, plugin = lk$plugin,
                            scale_basis = lk$scale_basis, sigma_kappa = fit$rate$sigma_kappa,
                            seed = seed, cr = fit$rate$cr)
  fit2 <- fit; fit2$rate <- rt2
  grid <- seq(min(x[-excl]), max(x[-excl]), length.out = 50)
  l1 <- link_curve(fit, component, grid = grid); l2 <- link_curve(fit2, component, grid = grid)
  max_diff <- max(abs(l1$level_mean - l2$level_mean))
  if (plot) {
    yl <- range(l1$level_ci, l2$level_ci)
    plot(grid, l1$level_mean, type = "n", lwd = 2, col = "black", ylim = yl, xlab = l1$var,
         ylab = "expected activity effect g",
         main = sprintf("Leverage: link with and without the %d most extreme actors", n_exclude))
    polygon(c(grid, rev(grid)), c(l1$level_ci[, 1], rev(l1$level_ci[, 2])), col = "grey85", border = NA)
    lines(grid, l1$level_mean, lwd = 2, col = "black")
    lines(grid, l2$level_mean, lwd = 2, col = "grey35", lty = 2)
    lines(grid, l2$level_ci[, 1], col = "grey35", lty = 3); lines(grid, l2$level_ci[, 2], col = "grey35", lty = 3)
    rug(x[-excl]); abline(v = 0, col = "grey40", lty = 2)
    legend("topleft", c("all actors", "extreme actors removed from the link"),
           col = c("black", "grey35"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  }
  if (verbose) {
    cat(sprintf("Actors removed from the link: %s\n", paste(excl, collapse = ", ")))
    cat(sprintf("  their link variable: %s\n", paste(round(x[excl], 2), collapse = ", ")))
    cat(sprintf("  their events sent  : %s  (median over all actors: %.0f)\n",
                paste(fit$N_out[excl], collapse = ", "), median(fit$N_out)))
    cat(sprintf("Largest change in the level over the common range: %.3f\n", max_diff))
  }
  invisible(list(excluded = excl, full = l1, reduced = l2, max_diff = max_diff))
}

#' Agreement of the plug-in summaries with their posterior means
#'
#' Compares the summaries evaluated at the posterior means of the aligned
#' positions (the plug-in values) with the posterior means of the summaries
#' computed draw by draw. For the self-affinity they differ by the posterior
#' covariance of \eqn{u_i} and \eqn{d_i}, which should be small and unrelated
#' to the number of events sent.
#'
#' @param fit an object of class `lsdynam` with \eqn{D \ge 1}{D >= 1}.
#' @param vars summaries to compare.
#' @param plot plot one against the other.
#' @param verbose print the table.
#' @return A list with the table of mean absolute differences and the
#'   posterior means, invisibly.
#' @examples
#' tr <- make_truth(n = 10, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 10, M = 200, seed = 2)
#' fit <- lsdynam_fit(ev, n = 10, link = "none", inference = "map", verbose = FALSE)
#' summary_agreement_check(fit, vars = c("pop", "rho"), plot = FALSE)
#' @export
summary_agreement_check <- function(fit, vars = c("r", "norm_u"), plot = TRUE, verbose = TRUE) {
  sd_draws <- summaries_draws(fit$choice$U, fit$choice$V, vars)
  pm <- apply(sd_draws, c(2, 3), mean)
  res <- do.call(rbind, lapply(seq_along(vars), function(j) {
    gap <- pm[, j] - fit$summaries[[vars[j]]]
    data.frame(var = vars[j], mean_abs_gap = mean(abs(gap)), sd_plugin = sd(fit$summaries[[vars[j]]]),
               rho_gap_outdeg = suppressWarnings(cor(gap, fit$N_out, method = "spearman")))
  }))
  if (plot) {
    op <- par(mfrow = c(1, length(vars))); on.exit(par(op))
    for (j in seq_along(vars)) {
      plot(fit$summaries[[vars[j]]], pm[, j], pch = 19, col = "grey30",
           xlab = paste(vars[j], "at posterior-mean positions"), ylab = paste("posterior mean of", vars[j]),
           main = vars[j])
      abline(0, 1, col = "black", lty = 2)
    }
  }
  if (verbose) print(res, digits = 3, row.names = FALSE)
  invisible(list(table = res, posterior_mean = pm))
}

#' Stability of the link across the modes of the choice module
#'
#' The posterior of the choice module can have local modes. For every
#' starting configuration of the posterior-mode search (`n_starts` in
#' [lsdynam_fit()]), this function reports the objective (negative
#' log-posterior) and its distance from the best start, the convergence code
#' (0: converged, 1: iteration limit), how similar the fitted structure is to
#' that of the best start (correlations of the self-affinities, of the
#' summaries in the link and of the centred latent affinities
#' \eqn{u_i^\top(v_j - \bar v)}{u_i'(v_j - vbar)}, all invariant), and the link
#' fitted on that start's positions (Laplace approximation): its evidence,
#' \eqn{R^2_a}, the contrast of each component and, for an additive link, the
#' evidence for each component given the others.
#'
#' Starts with convergence code 1 are not modes: increase `maxit`. If
#' converged starts with a clearly worse objective give a different link, the
#' conclusions rest on having found the best mode, and more starts are
#' needed.
#'
#' @param fit an object of class `lsdynam` fitted with at least two starts
#'   (`n_starts > 1`, or `n_starts = 1` and `spectral_start = TRUE`).
#' @param plot plot the components of the link fitted at each start.
#' @param verbose print the table.
#' @return A list with the table, ordered by objective, and the fitted
#'   components at each start, invisibly.
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", n_starts = 3, verbose = FALSE)
#' ms <- mode_stability_check(fit, plot = FALSE, verbose = FALSE)
#' ms$table[, c("start", "delta", "cor_r", "log_bf10")]
#' @export
mode_stability_check <- function(fit, plot = TRUE, verbose = TRUE) {
  mp <- fit$choice$map
  if (is.null(mp$start_par) || nrow(mp$start_par) < 2) stop("refit with n_starts > 1 (or spectral_start = TRUE)")
  lk <- fit$rate$link
  has_link <- !identical(lk$vars, "none")
  nc <- if (has_link) length(lk$vars) else 0
  n <- fit$n; D <- fit$D; S <- length(mp$zeta)
  best <- which.min(mp$starts$value)
  unpack <- function(k) choice_unpack(mp$start_par[k, ], n, D, S)
  aff <- function(pk) as.vector(pk$U %*% t(sweep(pk$V, 2, colMeans(pk$V))))
  pb <- unpack(best); ps_best <- position_summaries(pb$U, pb$V); aff_best <- aff(pb)
  rows <- vector("list", nrow(mp$starts)); curves <- vector("list", nrow(mp$starts))
  for (k in seq_len(nrow(mp$starts))) {
    pk <- unpack(k); ps <- position_summaries(pk$U, pk$V)
    row <- data.frame(start = mp$starts$start[k], objective = mp$starts$value[k],
                      delta = mp$starts$value[k] - mp$starts$value[best],
                      convergence = mp$starts$convergence[k],
                      cor_r = cor(ps$r, ps_best$r), cor_affinity = cor(aff(pk), aff_best),
                      zeta1 = if (S > 0) pk$zeta[1] else NA_real_, stringsAsFactors = FALSE)
    if (has_link) {
      for (v in setdiff(lk$vars, "r")) row[[paste0("cor_", v)]] <- cor(ps[[v]], ps_best[[v]])
      ch <- list(D = D, n = n, U_mean = pk$U, V_mean = pk$V, U = array(pk$U, c(1, n, D)),
                 V = array(pk$V, c(1, n, D)), zeta = matrix(pk$zeta, 1, S))
      rt <- lsdynam_rate_stage(NULL, n, sl = NULL, choice = ch, link = lk$vars, link_type = lk$type, df = lk$df,
                               winsorize_q = lk$winsorize_q, priors = fit$priors,
                               sigma_a = if (fit$rate$sigma_a_method == "eb") "eb" else fit$rate$sigma_a,
                               sigma_kappa = if (fit$rate$sigma_kappa_method == "eb") "eb" else fit$rate$sigma_kappa,
                               nu = fit$rate$nu, inference = "map", scale_basis = lk$scale_basis,
                               cr = fit$rate$cr, seed = 1)
      fk <- fit; fk$rate <- rt; fk$summaries <- ps
      row$sigma_a <- rt$sigma_a
      row$log_bf10 <- link_evidence(fk)$log_bf10_laplace
      row$R2 <- link_r2(fk)
      for (j in seq_len(nc)) {
        ct <- link_contrast(fk, component = j)
        row[[paste0("contrast_", lk$vars[j])]] <- ct$mean
        row[[paste0("lo_", lk$vars[j])]] <- ct$ci[1]
        row[[paste0("hi_", lk$vars[j])]] <- ct$ci[2]
        if (nc > 1) row[[paste0("log_bf10_", lk$vars[j], "_given_rest")]] <- link_evidence(fk, components = j)$log_bf10_laplace
      }
      curves[[k]] <- lapply(seq_len(nc), function(j) link_curve(fk, j))
    }
    rows[[k]] <- row
  }
  tab <- do.call(rbind, rows)
  ord <- order(tab$objective)
  if (verbose) print(format(tab[ord, ], digits = 3), row.names = FALSE)
  if (plot && has_link) {
    if (nc > 1) { op <- par(mfrow = c(1, nc)); on.exit(par(op), add = TRUE) }
    cols <- grDevices::gray(seq(0.35, 0.8, length.out = max(2, length(ord))))[seq_along(ord)]
    for (j in seq_len(nc)) {
      yl <- range(unlist(lapply(curves, function(cv) cv[[j]]$level_mean)))
      xl <- range(unlist(lapply(curves, function(cv) cv[[j]]$x)))
      vlab <- summary_label(lk$vars[j])
      plot(NA, xlim = xl, ylim = yl, xlab = vlab, ylab = "expected activity effect g",
           main = if (nc > 1) sprintf("Component in %s at each solution", vlab) else
             "Level of the link at each solution of the choice module")
      for (jj in rev(seq_along(ord))) {
        k <- ord[jj]
        lines(curves[[k]][[j]]$x, curves[[k]][[j]]$level_mean, col = if (jj == 1) "black" else cols[jj],
              lwd = if (jj == 1) 3 else 1.5)
      }
      abline(v = 0, col = "grey40", lty = 2)
      if (j == 1) legend("topleft", c("best start", "other starts (darker = closer to the best)"),
                         col = c("black", cols[2]), lwd = c(3, 1.5), bty = "n", cex = 0.8)
    }
  }
  invisible(list(table = tab[ord, ], curves = curves[ord]))
}

#' Calibration of the link by simulation
#'
#' Reference distribution of the link summaries when there is no
#' association: events are simulated from the fitted model, with the
#' estimated positions and coefficients and with the estimated activity
#' effects randomly reassigned to the actors (so that they are independent of
#' the positions but keep their distribution), and the two-stage procedure is
#' applied to each simulated sequence.
#'
#' With `keep` (for example `"pop"`), the association with these summaries is
#' kept and the rest removed: the link in `keep` is refitted on the positions
#' of `fit`, and the simulated activity effects are its fitted level plus the
#' residuals of the estimated activity effects, randomly reassigned. With
#' `keep = "pop"` this gives the reference distribution of the evidence for
#' excess self-affinity given popularity when only popularity is associated
#' with activity, which includes whatever the shrinkage of the sender
#' positions produces.
#'
#' @param fit an object of class `lsdynam` with a link.
#' @param n_rep number of simulated data sets.
#' @param grid optional grid (or list of grids, one per component) at which
#'   the levels are recorded.
#' @param inference inference for the replicates (`"map"` is much faster).
#' @param keep summaries whose association with activity is kept.
#' @param n_starts number of random starts for the replicates.
#' @param seed seed of the random number generator; replicate \eqn{b} uses
#'   `seed + b`.
#' @param verbose print a line per replicate.
#' @param ... further arguments for [lsdynam_fit()].
#' @return A list with, for every replicate, \eqn{R^2_a} `r2`, the log Bayes
#'   factor for the whole link `log_bf10`, for each component given the
#'   others `log_bf10_components`, the contrasts `contrast` and the levels of
#'   each component on its grid `levels`, together with the grids and the
#'   messages of failed replicates.
#' @examples
#' \donttest{
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' cal <- null_link_calibration(fit, n_rep = 3, verbose = FALSE)
#' cal$log_bf10
#' cal_pop <- null_link_calibration(fit, n_rep = 3, keep = "pop", verbose = FALSE)
#' cal_pop$log_bf10_components[, "rho"]
#' }
#' @export
null_link_calibration <- function(fit, n_rep = 20, grid = NULL, inference = "map", keep = NULL,
                                  n_starts = 1, seed = 1, verbose = TRUE, ...) {
  lk <- fit$rate$link
  nc <- length(lk$vars)
  grids <- lapply(seq_len(nc), function(j) {
    if (!is.null(grid)) { if (is.list(grid)) grid[[j]] else grid } else
      seq(quantile(lk$X[, j], 0.02), quantile(lk$X[, j], 0.98), length.out = 40)
  })
  names(grids) <- lk$vars
  a_hat <- colMeans(fit$rate$a)
  if (is.null(keep)) {
    base <- rep(0, fit$n); res <- a_hat
  } else {
    rk <- lsdynam_rate_stage(NULL, fit$n, sl = NULL, choice = fit$choice, link = keep, link_type = lk$type,
                             df = lk$df, winsorize_q = lk$winsorize_q, priors = fit$priors,
                             sigma_a = if (fit$rate$sigma_a_method == "eb") "eb" else fit$rate$sigma_a,
                             sigma_kappa = if (fit$rate$sigma_kappa_method == "eb") "eb" else fit$rate$sigma_kappa,
                             nu = fit$rate$nu, inference = "map", scale_basis = lk$scale_basis,
                             plugin = lk$plugin, cr = fit$rate$cr, seed = seed)
    base <- rk$map$mu_a + as.numeric(rk$link$B %*% rk$map$kappa_vec)    # fitted level of the kept link
    res <- rk$map$a - base                                              # residuals of the activity effects
  }
  theta0 <- list(xi = if (is.null(keep)) colMeans(fit$rate$xi) else rk$map$xi, zeta = colMeans(fit$choice$zeta),
                 U = fit$choice$U_mean, V = fit$choice$V_mean)
  levels <- lapply(grids, function(g) matrix(NA_real_, n_rep, length(g)))
  r2 <- numeric(n_rep); lbf <- numeric(n_rep)
  contrast <- matrix(NA_real_, n_rep, nc, dimnames = list(NULL, lk$vars))
  lbf_comp <- if (nc > 1) matrix(NA_real_, n_rep, nc, dimnames = list(NULL, lk$vars)) else NULL
  failed <- character(0)
  for (b in seq_len(n_rep)) {
    # a failing replicate is recorded as NA instead of stopping the calibration
    ok <- tryCatch({
      set.seed(seed + b)
      theta0$a <- base + sample(res)
      ev <- simulate_lsdynam(theta0, fit$n, stats = fit$stats, M = fit$M)
      fb <- lsdynam_fit(ev, fit$n, D = fit$D, stats = fit$stats, link = lk$vars, link_type = lk$type,
                        df = lk$df, winsorize_q = lk$winsorize_q, priors = fit$priors, scale_basis = lk$scale_basis,
                        plugin = lk$plugin, nu = fit$rate$nu,
                        sigma_kappa = if (fit$rate$sigma_kappa_method == "eb") "eb" else fit$rate$sigma_kappa,
                        sigma_a = if (fit$rate$sigma_a_method == "eb") "eb" else fit$rate$sigma_a,
                        inference = inference, n_starts = n_starts, seed = seed + b, verbose = FALSE, ...)
      for (j in seq_len(nc)) {
        lc <- suppressWarnings(link_curve(fb, j, grid = grids[[j]]))
        levels[[j]][b, ] <- lc$level_mean - lc$mu_a_mean
        contrast[b, j] <- link_contrast(fb, component = j)$mean
        if (nc > 1) lbf_comp[b, j] <- link_evidence(fb, components = j)$log_bf10_laplace
      }
      r2[b] <- link_r2(fb); lbf[b] <- link_evidence(fb)$log_bf10_laplace
      TRUE
    }, error = function(e) { failed <<- c(failed, sprintf("replicate %d: %s", b, conditionMessage(e))); FALSE })
    if (!ok) { r2[b] <- NA; lbf[b] <- NA; contrast[b, ] <- NA; if (nc > 1) lbf_comp[b, ] <- NA }
    if (verbose) cat(if (!ok) sprintf("  replicate %d / %d FAILED: %s\n", b, n_rep, failed[length(failed)]) else
      sprintf("  replicate %d / %d: R2 = %.3f, log BF10 = %.2f%s\n", b, n_rep, r2[b], lbf[b],
              if (nc > 1) paste0(", given the rest: ", paste(sprintf("%s %.2f", lk$vars, lbf_comp[b, ]), collapse = ", ")) else ""))
  }
  list(grid = grids, levels = levels, r2 = r2, log_bf10 = lbf, log_bf10_components = lbf_comp,
       contrast = contrast, keep = keep, vars = lk$vars, failed = failed)
}
