# Index of a component of the link, given by position or by name.
.component_index <- function(fit, component) {
  vars <- fit$rate$link$vars
  if (is.character(component)) {
    j <- match(component, vars)
    if (anyNA(j)) stop("unknown component '", component[is.na(j)][1], "'; the link has ",
                       paste(vars, collapse = ", "))
    return(j)
  }
  if (length(component) != 1 || component < 1 || component > length(vars)) stop("component out of range")
  as.integer(component)
}

#' Level and slope of a component of the link
#'
#' Posterior distribution of the level
#' \eqn{\mu_a + B_l(x)^\top \kappa_l}{mu_a + B_l(x)' kappa_l} of one component
#' of the link, and of its local slope, on a grid of raw (uncentred) values
#' \eqn{x} of the summary, under the plug-in or the cut posterior.
#'
#' @param fit an object of class `lsdynam` with a link.
#' @param component index or name of the component (for example `"pop"`).
#' @param grid values of the summary (default: `n_grid` values over the range
#'   of the actors' values).
#' @param n_grid number of grid values when `grid` is not given.
#' @param source `"plugin"` for the plug-in posterior, `"cut"` for the cut
#'   posterior (needs `n_cut > 0` in [lsdynam_fit()]).
#' @param level probability of the credible intervals.
#'
#' @return A list with the name of the summary `var`, the grid `x`, the
#'   actors' values `x_actors`, the posterior mean and credible intervals of
#'   the level (`level_mean`, `level_ci`) and of the slope (`slope_mean`,
#'   `slope_ci`), the draws of both (grid \eqn{\times}{x} draws), and the
#'   posterior mean of \eqn{\mu_a}.
#' @seealso [plot_link_level()], [link_contrast()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' lc <- link_curve(fit, "pop", n_grid = 5)
#' cbind(lc$x, lc$level_mean, lc$level_ci)
#' @export
link_curve <- function(fit, component = 1, grid = NULL, n_grid = 60,
                       source = c("plugin", "cut"), level = 0.95) {
  source <- match.arg(source)
  lk <- fit$rate$link
  if (identical(lk$vars, "none")) stop("this fit has no link")
  component <- .component_index(fit, component)
  comp <- lk$basis$spec$comps[[component]]
  cols <- lk$basis$spec$cols[[component]]
  x_act <- lk$X[, component]
  if (is.null(grid)) grid <- seq(min(x_act), max(x_act), length.out = n_grid)
  if (source == "cut") {
    if (is.null(fit$rate$cut)) stop("no cut-posterior draws in this fit (refit with n_cut > 0)")
    mu <- fit$rate$cut$mu_a; kap <- fit$rate$cut$kappa_vec[, cols, drop = FALSE]
  } else {
    mu <- fit$rate$mu_a; kap <- fit$rate$kappa_vec[, cols, drop = FALSE]
  }
  h <- 1e-4 * max(1e-8, diff(range(x_act)))
  Bg <- eval_link_component(comp, grid)
  dB <- (eval_link_component(comp, grid + h) - eval_link_component(comp, grid - h)) / (2 * h)
  lev <- sweep(Bg %*% t(kap), 2, mu, "+")          # n_grid x n_draws
  slo <- dB %*% t(kap)
  pr <- c((1 - level) / 2, 1 - (1 - level) / 2)
  list(var = comp$var, x = grid, center = comp$center, x_actors = x_act,
       level_mean = rowMeans(lev), level_ci = t(apply(lev, 1, quantile, probs = pr)),
       slope_mean = rowMeans(slo), slope_ci = t(apply(slo, 1, quantile, probs = pr)),
       level_draws = lev, slope_draws = slo, mu_a_mean = mean(mu), source = source, level = level)
}

#' Contrast of a component of the link
#'
#' Posterior distribution of the difference between the levels of one
#' component of the link at two values of its summary, by default the 10\% and
#' 90\% quantiles of the actors' values. For latent popularity this is how
#' much higher, on the log scale, the typical rate of an actor at the upper
#' end is than that of an actor at the lower end, other things being equal.
#'
#' @inheritParams link_curve
#' @param x_lo,x_hi the two values of the summary (default: its 10\% and 90\%
#'   quantiles over the actors).
#'
#' @return A list with the two values `x`, the posterior mean `mean`, the
#'   credible interval `ci`, the posterior probability that the contrast is
#'   positive `prob_positive` and the draws `draws`.
#' @seealso [link_curve()], [summary.lsdynam()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' ct <- link_contrast(fit, component = "pop")
#' c(ct$mean, ct$ci)
#' @export
link_contrast <- function(fit, component = 1, x_lo = NULL, x_hi = NULL, source = c("plugin", "cut"),
                          level = 0.95) {
  source <- match.arg(source)
  component <- .component_index(fit, component)
  x_act <- fit$rate$link$X[, component]
  if (is.null(x_lo)) x_lo <- quantile(x_act, 0.10, names = FALSE)
  if (is.null(x_hi)) x_hi <- quantile(x_act, 0.90, names = FALSE)
  lc <- link_curve(fit, component, grid = c(x_lo, x_hi), source = source)
  dd <- lc$level_draws[2, ] - lc$level_draws[1, ]
  pr <- c((1 - level) / 2, 1 - (1 - level) / 2)
  list(x = c(x_lo, x_hi), mean = mean(dd), ci = quantile(dd, pr, names = FALSE),
       prob_positive = mean(dd > 0), draws = dd)
}

#' Evidence for the link
#'
#' Bayes factor for the model with the link against the same model with the
#' link coefficients set to zero, given the plug-in summaries and
#' \eqn{\sigma_a}. It is computed in two ways: by the Savage-Dickey density
#' ratio, with a Normal approximation to the posterior of the coefficients
#' computed from the draws, and by the ratio of Laplace approximations to the
#' marginal likelihoods. Their agreement is a check on both.
#'
#' With `components`, only those components are tested and the others stay in
#' both models: for the additive link `c("pop", "rho")`,
#' `components = "rho"` gives the evidence for excess self-affinity given
#' latent popularity. A Bayes factor for a point hypothesis depends on the
#' prior scale of the coefficients: `sigma_kappa_grid` returns the Laplace
#' version for several prior standard deviations. If `sigma_kappa` was
#' estimated (`"eb"`) the ratio is maximised over the prior scale and is not a
#' Bayes factor; calibrate it with [null_link_calibration()].
#'
#' @param fit an object of class `lsdynam` with a link.
#' @param sigma_kappa_grid optional prior standard deviations of the link
#'   coefficients at which to recompute the Laplace version.
#' @param components indices or names of the components to test (default:
#'   the whole link).
#'
#' @return A list with the log Bayes factors `log_bf10_sd` (Savage-Dickey)
#'   and `log_bf10_laplace` (positive values favour the link), the number of
#'   coefficients tested `K`, the components tested, the prior variance and,
#'   if requested, the log Bayes factors on the grid of prior scales `grid`.
#' @seealso [link_comparison_table()], [null_link_calibration()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' link_evidence(fit)$log_bf10_laplace
#' # excess self-affinity given popularity
#' link_evidence(fit, components = "rho", sigma_kappa_grid = c(0.5, 1, 2))$grid
#' @export
link_evidence <- function(fit, sigma_kappa_grid = NULL, components = NULL) {
  rt <- fit$rate
  K <- ncol(rt$kappa_vec)
  if (K == 0) stop("this fit has no link")
  spec <- rt$link$basis$spec
  if (is.null(components)) {
    cols <- seq_len(K); comp <- spec$vars
  } else {
    ci <- if (is.character(components)) match(components, spec$vars) else as.integer(components)
    if (anyNA(ci) || any(ci < 1 | ci > length(spec$vars))) stop("link_evidence: unknown component(s)")
    cols <- unlist(spec$cols[ci]); comp <- spec$vars[ci]
  }
  Kt <- length(cols)
  kd <- rt$kappa_vec[, cols, drop = FALSE]
  m <- colMeans(kd); C <- stats::cov(kd)
  Cc <- chol(make_pd(C))
  log_post0  <- -0.5 * Kt * log(2 * pi) - sum(log(diag(Cc))) - 0.5 * sum(backsolve(Cc, m, transpose = TRUE)^2)
  log_prior0 <- -0.5 * Kt * log(2 * pi * rt$priors$sigma2_kappa)
  B0 <- rt$link$B[, -cols, drop = FALSE]               # the link without the tested components
  fit0 <- rate_mod_fit(rt$cr, B0, rt$priors, sigma2_a = rt$sigma_a^2, nu = rt$nu)
  grid <- NULL
  if (!is.null(sigma_kappa_grid)) {
    grid <- vapply(sigma_kappa_grid, function(sk) {
      hp <- rt$priors; hp$sigma2_kappa <- sk^2
      e1 <- rate_mod_fit(rt$cr, rt$link$B, hp, sigma2_a = rt$sigma_a^2, nu = rt$nu, start = rt$map$par)$log_evidence
      e0 <- if (ncol(B0) == 0) fit0$log_evidence else
        rate_mod_fit(rt$cr, B0, hp, sigma2_a = rt$sigma_a^2, nu = rt$nu)$log_evidence
      e1 - e0
    }, numeric(1))
    names(grid) <- paste0("sigma_kappa=", sigma_kappa_grid)
  }
  list(log_bf10_sd = log_prior0 - log_post0,
       log_bf10_laplace = rt$map$log_evidence - fit0$log_evidence, K = Kt, components = comp,
       sigma2_kappa = rt$priors$sigma2_kappa, sigma_kappa_method = rt$sigma_kappa_method, grid = grid)
}

#' Proportion of variation explained by the link
#'
#' The proportion \eqn{R^2_a} of the variation of the activity effects
#' between actors that the link accounts for, averaged over the posterior
#' draws (Gelman and Pardoe, 2006).
#'
#' @param fit an object of class `lsdynam`.
#' @return A number (0 for a model without link).
#' @references Gelman, A. and Pardoe, I. (2006). Bayesian measures of
#'   explained variance and pooling in multilevel (hierarchical) models.
#'   *Technometrics*, 48, 241--251.
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' link_r2(fit)
#' @export
link_r2 <- function(fit) {
  rt <- fit$rate
  if (ncol(rt$kappa_vec) == 0) return(0)
  fitted <- sweep(rt$kappa_vec %*% t(rt$link$B), 1, rt$mu_a, "+")   # n_draws x n
  1 - mean(apply(rt$a - fitted, 1, var)) / mean(apply(rt$a, 1, var))
}

#' Compare links on the same positions
#'
#' Fits the rate module with several links on the positions of `fit` (Laplace
#' approximation, \eqn{\sigma_a} set by maximum marginal likelihood in each
#' model) and compares them by their log marginal likelihoods relative to the
#' model without link. Differences between rows compare links directly; for
#' nested additive links they are log Bayes factors for the components that
#' differ, for example `"pop + rho"` minus `"pop"` for excess self-affinity
#' given popularity.
#'
#' @param fit an object of class `lsdynam` with \eqn{D \ge 1}{D >= 1}.
#' @param specs named list; each element is a link (a character vector of
#'   summaries) or a list of arguments of [lsdynam_rate_stage()], such as
#'   `list(link = "r", link_type = "linear")`. Settings that a specification
#'   does not give are those of `fit`.
#' @param df number of basis functions per component (default: that of
#'   `fit`).
#' @param level probability of the credible intervals of the contrasts.
#'
#' @return A list with `table` (one row per link: number of coefficients,
#'   \eqn{\sigma_a}, log marginal likelihood relative to the model without
#'   link, \eqn{R^2_a} and contrasts), the fitted models `fits`, and
#'   \eqn{\sigma_a} and the log marginal likelihood of the model without link.
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' link_comparison_table(fit, list("pop + rho" = c("pop", "rho"), pop = "pop",
#'                                 rho = "rho", r = "r"))$table
#' @export
link_comparison_table <- function(fit, specs, df = NULL, level = 0.95) {
  lk <- fit$rate$link
  if (is.null(df)) df <- lk$df
  dflt <- function(x, d) if (is.null(x)) d else x
  refit <- function(sp) {
    args <- list(events_df = NULL, n = fit$n, sl = NULL, choice = fit$choice, link = "none", df = df,
                 winsorize_q = dflt(lk$winsorize_q, 1), scale_basis = dflt(lk$scale_basis, TRUE),
                 plugin = dflt(lk$plugin, "positions"), nu = dflt(fit$rate$nu, Inf),
                 sigma_kappa = if (identical(fit$rate$sigma_kappa_method, "eb")) "eb" else fit$rate$sigma_kappa,
                 priors = fit$priors, sigma_a = "eb", inference = "map", cr = fit$rate$cr, seed = 1)
    if (is.list(sp)) args[names(sp)] <- sp else args$link <- sp
    rt <- do.call(lsdynam_rate_stage, args)
    f <- fit; f$rate <- rt; f
  }
  f0 <- refit("none")
  fits <- lapply(specs, refit)
  tab <- do.call(rbind, lapply(names(fits), function(nm) {
    f <- fits[[nm]]; vars <- f$rate$link$vars
    ct <- lapply(seq_along(vars), function(k) link_contrast(f, component = k, level = level))
    data.frame(link = nm, K = ncol(f$rate$kappa_vec), sigma_a = f$rate$sigma_a,
               log_ml_vs_M1 = f$rate$log_evidence - f0$rate$log_evidence, R2 = link_r2(f),
               contrast_10_90 = paste(vapply(seq_along(vars), function(k)
                 sprintf("%s %.2f [%.2f, %.2f]", vars[k], ct[[k]]$mean, ct[[k]]$ci[1], ct[[k]]$ci[2]),
                 character(1)), collapse = "; "),
               stringsAsFactors = FALSE)
  }))
  list(table = tab, fits = fits, sigma_a_M1 = f0$rate$sigma_a, log_evidence_M1 = f0$rate$log_evidence)
}
