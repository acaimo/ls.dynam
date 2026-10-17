#' Plot a component of the link
#'
#' `plot_link_level()` plots the posterior mean and 95\% band of the level of
#' one component of the link under the plug-in posterior (shaded) and, if
#' available, the cut posterior (dashed), with a rug of the actors' values and
#' a vertical line at zero. `plot_link_slope()` plots its local slope.
#'
#' @inheritParams link_curve
#' @param show_cut add the cut posterior if the fit has it.
#' @param truth optional function of the raw summary giving the true level
#'   (simulations).
#' @param col colour of the posterior mean.
#' @param band colour of the plug-in band.
#' @param main,xlab title and label of the horizontal axis.
#' @param centred_axis plot against the centred summary.
#' @param legend_pos position of the legend (`NULL`: no legend).
#' @param ... further graphical parameters for [plot()].
#'
#' @return The output of [link_curve()] for the plug-in and the cut
#'   posterior, invisibly.
#' @seealso [plot.lsdynam()], [link_curve()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' plot_link_level(fit, "pop")
#' plot_link_slope(fit, "pop")
#' @export
plot_link_level <- function(fit, component = 1, show_cut = TRUE, truth = NULL,
                            col = "black", band = "grey80", main = NULL, centred_axis = FALSE, xlab = NULL,
                            legend_pos = "topleft", ...) {
  lc <- link_curve(fit, component, source = "plugin")
  lcc <- if (show_cut && !is.null(fit$rate$cut)) link_curve(fit, component, grid = lc$x, source = "cut") else NULL
  x <- if (centred_axis) lc$x - lc$center else lc$x
  yl <- range(lc$level_ci, if (!is.null(lcc)) lcc$level_ci, if (!is.null(truth)) truth(lc$x))
  if (is.null(main)) main <- sprintf("Level of the position-rate link (%s)", lc$var)
  plot(x, lc$level_mean, type = "n", ylim = yl,
       xlab = if (!is.null(xlab)) xlab else if (centred_axis) paste0(lc$var, " (centred)") else lc$var,
       ylab = "expected activity effect g", main = main, ...)
  polygon(c(x, rev(x)), c(lc$level_ci[, 1], rev(lc$level_ci[, 2])), col = band, border = NA)
  if (!is.null(lcc)) {
    lines(x, lcc$level_ci[, 1], lty = 2, col = "grey30"); lines(x, lcc$level_ci[, 2], lty = 2, col = "grey30")
    lines(x, lcc$level_mean, lty = 3, col = "grey30")
  }
  lines(x, lc$level_mean, lwd = 2, col = col)
  if (!is.null(truth)) lines(x, truth(lc$x), col = "black", lwd = 2, lty = 4)
  abline(h = lc$mu_a_mean, col = "grey60", lty = 3)
  abline(v = if (centred_axis) -lc$center else 0, col = "grey40", lty = 2)
  rug(if (centred_axis) lc$x_actors - lc$center else lc$x_actors, col = "grey40")
  leg <- c("posterior mean (plug-in)", "95% band (plug-in)"); lt <- c(1, NA); lw <- c(2, NA); cl <- c(col, band)
  fl <- c(NA, band); bd <- c(NA, NA)
  if (!is.null(lcc)) { leg <- c(leg, "cut posterior: mean and 95% band"); lt <- c(lt, 2); lw <- c(lw, 1); cl <- c(cl, "grey30"); fl <- c(fl, NA); bd <- c(bd, NA) }
  if (!is.null(truth)) { leg <- c(leg, "truth"); lt <- c(lt, 4); lw <- c(lw, 2); cl <- c(cl, "black"); fl <- c(fl, NA); bd <- c(bd, NA) }
  if (!is.null(legend_pos)) legend(legend_pos, legend = leg, lty = lt, lwd = lw, col = cl, fill = fl, border = bd, bty = "n", cex = 0.8)
  invisible(list(plugin = lc, cut = lcc))
}

#' @rdname plot_link_level
#' @export
plot_link_slope <- function(fit, component = 1, show_cut = TRUE, col = "black", band = "grey80", main = NULL, xlab = NULL, ...) {
  lc <- link_curve(fit, component, source = "plugin")
  lcc <- if (show_cut && !is.null(fit$rate$cut)) link_curve(fit, component, grid = lc$x, source = "cut") else NULL
  yl <- range(lc$slope_ci, if (!is.null(lcc)) lcc$slope_ci, 0)
  if (is.null(main)) main <- sprintf("Local slope of the link (%s)", lc$var)
  plot(lc$x, lc$slope_mean, type = "n", ylim = yl, xlab = if (is.null(xlab)) lc$var else xlab,
       ylab = "local slope g'", main = main, ...)
  polygon(c(lc$x, rev(lc$x)), c(lc$slope_ci[, 1], rev(lc$slope_ci[, 2])), col = band, border = NA)
  if (!is.null(lcc)) { lines(lc$x, lcc$slope_ci[, 1], lty = 2, col = "grey30"); lines(lc$x, lcc$slope_ci[, 2], lty = 2, col = "grey30") }
  lines(lc$x, lc$slope_mean, lwd = 2, col = col)
  abline(h = 0, col = "grey50", lty = 2); abline(v = 0, col = "grey40", lty = 2)
  rug(lc$x_actors, col = "grey40")
  invisible(list(plugin = lc, cut = lcc))
}

#' Plot the latent positions for one dimension
#'
#' For \eqn{D = 1}, plots the sender position \eqn{u_i} against the centred
#' receiver position \eqn{v_i - \bar v}{v_i - vbar} of each actor, with 95\%
#' posterior regions. Actors in the first and third quadrants have positive
#' self-affinity (filled points), those in the second and fourth negative
#' (open points); the grey hyperbolas are level sets of the self-affinity.
#' Point size increases with the number of events sent. The common sign of
#' the positions is not identified; it is chosen so that the mean sender
#' position is positive.
#'
#' @param fit an object of class `lsdynam` with \eqn{D = 1}.
#' @param labels label the points with the actor numbers.
#' @param ellipses draw the posterior regions.
#' @param level probability of the posterior regions.
#' @param contours levels of the self-affinity for the hyperbolas.
#' @param legend_pos position of the legend (`NULL`: no legend).
#' @param ... further graphical parameters for [plot()].
#' @return `NULL`, invisibly.
#' @examples
#' tr <- make_truth(n = 15, D = 1, mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, link = "none", inference = "map", verbose = FALSE)
#' plot_positions(fit)
#' @export
plot_positions <- function(fit, labels = FALSE, ellipses = TRUE, level = 0.95, contours = c(-2, -1, -0.5, 0.5, 1, 2),
                           legend_pos = "topleft", ...) {
  if (fit$D != 1) stop("plot_positions is for D = 1; for D >= 2 plot the invariant summaries")
  ch <- fit$choice
  sg <- if (mean(ch$U_mean[, 1]) < 0) -1 else 1
  u <- sg * ch$U_mean[, 1]; vbar <- mean(sg * ch$V_mean[, 1]); dv <- sg * ch$V_mean[, 1] - vbar
  Ud <- sg * ch$U[, , 1]; Vd <- sg * (ch$V[, , 1] - rowMeans(ch$V[, , 1]))
  xl <- range(u, if (ellipses) quantile(Ud, c(0.005, 0.995))); yl <- range(dv, if (ellipses) quantile(Vd, c(0.005, 0.995)))
  plot(u, dv, type = "n", xlim = xl, ylim = yl, xlab = expression(u[i] ~ "(sender position)"),
       ylab = expression(v[i] - bar(v) ~ "(receiver position, centred)"),
       main = "Latent positions (D = 1)", ...)
  ug <- seq(xl[1], xl[2], length.out = 400)
  for (cc in contours) {
    for (sgn in c(-1, 1)) {
      uu <- ug[sgn * ug > 1e-3]
      lines(uu, cc / uu, col = "grey88", lty = 1)
    }
  }
  if (ellipses) {
    th <- seq(0, 2 * pi, length.out = 60); circ <- cbind(cos(th), sin(th)); q <- qchisq(level, 2)
    for (i in seq_len(fit$n)) {
      S <- cov(cbind(Ud[, i], Vd[, i])); eg <- eigen(S, symmetric = TRUE)
      ell <- circ %*% diag(sqrt(pmax(eg$values, 0) * q), 2) %*% t(eg$vectors)
      lines(ell[, 1] + u[i], ell[, 2] + dv[i], col = "grey75")
    }
  }
  abline(h = 0, v = 0, col = "grey25", lty = 2, lwd = 1.2)
  cex <- 0.5 + 1.2 * sqrt(fit$N_out / max(fit$N_out))
  points(u, dv, pch = ifelse(u * dv > 0, 19, 21), cex = cex, col = "black", bg = "white")
  if (labels) text(u, dv, labels = seq_len(fit$n), pos = 1, cex = 0.5)
  if (!is.null(legend_pos)) legend(legend_pos, c("r > 0 (quadrants I, III)", "r < 0 (quadrants II, IV)"), pch = c(19, 21),
         col = "black", pt.bg = "white", bty = "n", cex = 0.8)
  invisible(NULL)
}

#' Plot the activity effects
#'
#' Posterior means and credible intervals of the activity effects
#' \eqn{a_i}, sorted, or against their true values in a simulation.
#'
#' @param fit an object of class `lsdynam` (or its rate module `fit$rate`).
#' @param true_a optional true activity effects (simulations).
#' @param level probability of the credible intervals.
#' @param labels label the points with the actor numbers (with `true_a`).
#' @param main title.
#' @param ... further graphical parameters for [plot()].
#' @return A list with the posterior means `mean` and the credible intervals
#'   `ci` (2 \eqn{\times n}{x n}), invisibly.
#' @examples
#' tr <- make_truth(n = 15, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, link = "none", inference = "map", verbose = FALSE)
#' plot_activity_effects(fit, true_a = tr$theta$a)
#' @export
plot_activity_effects <- function(fit, true_a = NULL, level = 0.95, labels = FALSE, main = NULL, ...) {
  rate_post <- if (inherits(fit, "lsdynam")) fit$rate else fit
  a_draws <- rate_post$a
  n <- ncol(a_draws)
  a_mean <- colMeans(a_draws)
  q_lo <- (1 - level) / 2
  a_ci <- apply(a_draws, 2, quantile, probs = c(q_lo, 1 - q_lo))
  if (!is.null(true_a)) {
    plot(true_a, a_mean, pch = 19, col = "grey30", xlim = range(true_a), ylim = range(c(a_ci, true_a)),
         xlab = expression("true " * a[i]), ylab = expression("posterior mean of " * a[i]),
         main = if (is.null(main)) "Recovery of the activity effects" else main, ...)
    segments(true_a, a_ci[1, ], true_a, a_ci[2, ], col = "grey70", lwd = 1.5)
    points(true_a, a_mean, pch = 19, col = "grey30")
    abline(0, 1, col = "black", lty = 2)
    if (labels) text(true_a, a_mean, labels = seq_len(n), pos = 3, cex = 0.5)
    covered <- true_a >= a_ci[1, ] & true_a <= a_ci[2, ]
    legend("topleft", legend = sprintf("%d of %d %.0f%% intervals cover the truth", sum(covered), n, 100 * level),
           bty = "n", cex = 0.8)
  } else {
    ord <- order(a_mean)
    plot(seq_len(n), a_mean[ord], pch = 19, col = "grey30", ylim = range(a_ci),
         xlab = "actor (sorted by posterior mean)", ylab = expression(a[i]),
         main = if (is.null(main)) "Activity effects" else main, ...)
    segments(seq_len(n), a_ci[1, ord], seq_len(n), a_ci[2, ord], col = "grey70", lwd = 1.5)
    points(seq_len(n), a_mean[ord], pch = 19, col = "grey30")
  }
  invisible(list(mean = a_mean, ci = a_ci))
}

#' Trace plots
#'
#' Trace plots of the draws of the coefficients (\eqn{\xi}, \eqn{\zeta},
#' \eqn{\mu_a}, the link coefficients and the log-posterior of the choice
#' module) or of the activity effects of selected actors.
#'
#' @param fit an object of class `lsdynam`.
#' @param which `"parameters"` or `"activity"`.
#' @param actors actors whose activity effects are plotted (default: the
#'   `max_actors` actors with the fewest events sent).
#' @param max_actors number of actors plotted by default.
#' @return `NULL`, invisibly.
#' @examples
#' tr <- make_truth(n = 10, D = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 10, M = 200, seed = 2)
#' fit <- lsdynam_fit(ev, n = 10, link = "none", n_iter = 400, burn_in = 200, thin = 1,
#'                    verbose = FALSE)
#' plot_traces(fit)
#' @export
plot_traces <- function(fit, which = c("parameters", "activity"), actors = NULL, max_actors = 12) {
  which <- match.arg(which)
  series <- list()
  if (which == "parameters") {
    for (nm in colnames(fit$rate$xi)) series[[paste0("xi[", nm, "]")]] <- fit$rate$xi[, nm]
    for (nm in colnames(fit$choice$zeta)) series[[paste0("zeta[", nm, "]")]] <- fit$choice$zeta[, nm]
    series[["mu_a"]] <- fit$rate$mu_a
    for (nm in colnames(fit$rate$kappa_vec)) series[[paste0("kappa[", nm, "]")]] <- fit$rate$kappa_vec[, nm]
    if (!is.null(fit$choice$log_post)) series[["log-posterior (choice)"]] <- fit$choice$log_post
  } else {
    if (is.null(actors)) actors <- order(fit$N_out)[seq_len(min(max_actors, fit$n))]
    for (i in actors) series[[paste0("a[", i, "]")]] <- fit$rate$a[, i]
  }
  k <- length(series)
  nc <- ceiling(sqrt(k)); nr <- ceiling(k / nc)
  op <- par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0)); on.exit(par(op))
  for (nm in names(series)) plot(series[[nm]], type = "l", xlab = "draw", ylab = "", main = nm, col = "grey20")
  invisible(NULL)
}
