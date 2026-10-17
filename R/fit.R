#' Prior distributions
#'
#' Variances of the Normal prior distributions of a latent space DyNAM:
#' \eqn{\xi \sim N(0, \sigma^2_\xi I)}{xi ~ N(0, sigma2_xi I)},
#' \eqn{\zeta \sim N(0, \sigma^2_\zeta I)}{zeta ~ N(0, sigma2_zeta I)},
#' \eqn{u_i \sim N(0, \sigma^2_u I)}{u_i ~ N(0, sigma2_u I)},
#' \eqn{v_i \sim N(0, \sigma^2_v I)}{v_i ~ N(0, sigma2_v I)},
#' \eqn{\mu_a \sim N(0, \sigma^2_\mu)}{mu_a ~ N(0, sigma2_mu)} and
#' \eqn{\kappa \sim N(0, \sigma^2_\kappa I)}{kappa ~ N(0, sigma2_kappa I)}
#' for the coefficients of the link. Since the basis functions of the link are
#' scaled to unit standard deviation across actors, \eqn{\sigma_\kappa} is
#' the prior standard deviation of the contribution of one basis function to
#' the activity effect of a typical actor.
#'
#' @param sigma2_xi,sigma2_zeta prior variances of the coefficients of the
#'   rate and choice statistics.
#' @param sigma2_u,sigma2_v prior variances of the coordinates of the sender
#'   and receiver positions.
#' @param sigma2_mu prior variance of the average activity effect.
#' @param sigma2_kappa prior variance of the coefficients of the link.
#'
#' @return A list of prior variances.
#' @examples
#' lsdynam_priors(sigma2_kappa = 0.25)
#' @export
lsdynam_priors <- function(sigma2_xi = 4, sigma2_zeta = 4, sigma2_u = 1, sigma2_v = 1,
                           sigma2_mu = 100, sigma2_kappa = 1) {
  list(sigma2_xi = sigma2_xi, sigma2_zeta = sigma2_zeta, sigma2_u = sigma2_u,
       sigma2_v = sigma2_v, sigma2_mu = sigma2_mu, sigma2_kappa = sigma2_kappa)
}

#' Fit a latent space DyNAM with a position-rate link
#'
#' Fits the model in two stages (modular inference): the choice module, with
#' the latent positions, from the choice of receivers alone
#' ([lsdynam_choice_stage()]), and then the rate module given summaries of the
#' positions ([lsdynam_rate_stage()]).
#'
#' Three models of the paper are special cases:
#' * \eqn{\mathcal{M}_0}{M0}: `D = 0`, `link = "none"`, a DyNAM with random
#'   activity effects and no latent positions;
#' * \eqn{\mathcal{M}_1}{M1}: `D >= 1`, `link = "none"`, the latent space
#'   DyNAM without link;
#' * \eqn{\mathcal{M}_2}{M2}: `D >= 1` with a link, by default the additive
#'   link in latent popularity and excess self-affinity, `c("pop", "rho")`.
#'
#' With `inference = "mcmc"` both stages are simulated by MALA, with the
#' summaries plugged in at the posterior means of the aligned positions. With
#' `inference = "map"` both stages use the Laplace approximation at the
#' posterior mode: this is fast, but the summaries are plugged in at the mode,
#' which attenuates the link when actors have few events. The analysis of the
#' paper used `inference = "mcmc"`, `n_iter = 12000`, `burn_in = 2000`,
#' `thin = 5`, `n_starts = 30` and `n_cut = 300`.
#'
#' @param events_df data frame of events ordered by time, with columns
#'   `sender` and `receiver` (actor labels \eqn{1, \dots, n}) and `time`.
#' @param n number of actors.
#' @param D dimension of the latent space; `D = 0` gives the model without
#'   latent positions.
#' @param stats names of the statistics, see [compute_stats()]: names starting
#'   with `rate_` enter the rate sub-process, names starting with `choice_` the
#'   choice sub-process.
#' @param link `"none"`, or the summaries of the positions in the link
#'   (columns of [position_summaries()]), one additive component each:
#'   `c("pop", "rho")` (default) for latent popularity and excess
#'   self-affinity, `"r"` for the self-affinity alone.
#' @param link_type `"spline"` (natural cubic splines) or `"linear"`.
#' @param df number of basis functions per component of the link.
#' @param winsorize_q if smaller than 1, each link variable is capped at its
#'   quantiles `1 - winsorize_q` and `winsorize_q` before the basis is built.
#' @param priors prior variances, see [lsdynam_priors()].
#' @param sigma_a standard deviation of the activity effects about the link:
#'   `"eb"` to set it by maximising the Laplace approximation to the marginal
#'   likelihood of the rate module (empirical Bayes), or a number.
#' @param nu degrees of freedom of Student-t activity effects (`Inf`: Normal).
#' @param inference `"mcmc"` (posterior simulation) or `"map"` (Laplace
#'   approximation).
#' @param sampler `"mala"` (Metropolis-adjusted Langevin algorithm,
#'   recommended) or `"rwm"` (random-walk Metropolis).
#' @param n_iter,burn_in,thin number of iterations of the sampler (including
#'   the burn-in), length of the burn-in and thinning interval; in
#'   [lsdynam_fit()] they apply to the choice module.
#' @param n_iter_rate,burn_in_rate,thin_rate the same for the rate module
#'   (by default those of the choice module).
#' @param n_starts number of random starting configurations for the posterior
#'   mode of the choice module, which has local modes.
#' @param n_cut number of draws of the positions used for the cut posterior
#'   (0: none). See [lsdynam_rate_stage()].
#' @param cut_refit_sigma if `TRUE`, `sigma_a` is re-estimated for every draw
#'   of the cut posterior.
#' @param rate_event logical vector over the events, `TRUE` for the events
#'   that count as a sending decision in the rate module (default all). With
#'   messages to several recipients stored as several events, set it to
#'   `TRUE` for the first event of each message only (see
#'   [prepare_events()]).
#' @param plugin `"positions"` evaluates the summaries at the posterior means
#'   of the aligned positions (at the mode for `inference = "map"`);
#'   `"summary"` uses the posterior means of the summaries computed draw by
#'   draw.
#' @param scale_basis scale the basis functions of the link to unit standard
#'   deviation across actors (recommended).
#' @param sigma_kappa prior standard deviation of the link coefficients:
#'   `NULL` uses `priors$sigma2_kappa`, a number sets it, `"eb"` estimates it
#'   by maximum marginal likelihood.
#' @param spectral_start add a starting configuration from the singular value
#'   decomposition of the matrix of event counts.
#' @param maxit iteration limit of BFGS for each start.
#' @param choice a previous result of [lsdynam_choice_stage()] (or
#'   `fit$choice`), to fit the rate module on the same positions, for example
#'   models with and without link.
#' @param sl statistics computed by [compute_stats()] for `events_df`
#'   (computed if `NULL`).
#' @param keep_stats keep the statistics in the fitted object.
#' @param seed seed of the random number generator; the functions call
#'   [set.seed()].
#' @param verbose print progress.
#'
#' @return An object of class `lsdynam`: a list with
#'   * `choice`: the result of [lsdynam_choice_stage()];
#'   * `rate`: the result of [lsdynam_rate_stage()];
#'   * `summaries`: [position_summaries()] at the plug-in positions;
#'   * `n`, `D`, `M`, `stats`, `priors`, the numbers of events sent and
#'     received by each actor `N_out` and `N_in`, the statistics
#'     `stats_list` (if `keep_stats`), and the posterior modes of the two
#'     modules `map_fit`.
#' @seealso [summary.lsdynam()], [plot.lsdynam()], [link_contrast()],
#'   [link_evidence()], [predictive_scores()], [mcmc_diagnostics()].
#' @examples
#' f <- function(x) 1.5 * (plogis(3 * x) - 0.5)
#' tr <- make_truth(n = 20, D = 1, link_fun = list(f, NULL), link_var = c("pop", "rho"),
#'                  mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 20, M = 800, seed = 2)
#' fit <- lsdynam_fit(ev, n = 20, D = 1, inference = "map", verbose = FALSE)
#' fit
#' link_contrast(fit, component = "pop")$mean
#' # the model without link, on the same positions
#' fit1 <- lsdynam_fit(ev, n = 20, link = "none", inference = "map",
#'                     choice = fit$choice, verbose = FALSE)
#' \donttest{
#' # posterior simulation with the cut posterior (short chains)
#' fit_mcmc <- lsdynam_fit(ev, n = 20, D = 1, n_iter = 1500, burn_in = 500, thin = 2,
#'                         n_starts = 2, n_cut = 50, verbose = FALSE)
#' summary(fit_mcmc)
#' }
#' @export
lsdynam_fit <- function(events_df, n, D = 1,
                        stats = c("rate_indeg", "choice_recip"),
                        link = if (D > 0) c("pop", "rho") else "none",
                        link_type = c("spline", "linear"), df = 4, winsorize_q = 1,
                        priors = lsdynam_priors(), sigma_a = "eb", nu = Inf,
                        inference = c("mcmc", "map"), sampler = c("mala", "rwm"),
                        n_iter = 6000, burn_in = 1000, thin = 5,
                        n_iter_rate = n_iter, burn_in_rate = burn_in, thin_rate = thin,
                        n_starts = 1, n_cut = 0, cut_refit_sigma = FALSE, rate_event = NULL,
                        plugin = c("positions", "summary"), scale_basis = TRUE, sigma_kappa = NULL,
                        spectral_start = TRUE, maxit = 5000,
                        choice = NULL, sl = NULL, keep_stats = FALSE, seed = 1, verbose = TRUE) {
  inference <- match.arg(inference); sampler <- match.arg(sampler); link_type <- match.arg(link_type)
  plugin <- match.arg(plugin)
  known <- c("rate_outdeg", "rate_indeg", "choice_outdeg", "choice_indeg", "choice_inertia", "choice_recip")
  if (!all(stats %in% known)) {
    # compute_stats() silently drops names it does not know
    stop("unknown statistic(s): ", paste(setdiff(stats, known), collapse = ", "))
  }
  if (is.null(sl)) sl <- compute_stats(events_df, n, stats = stats)
  if (is.null(choice)) {
    if (verbose) cat(sprintf("Stage 1 (choice module, D = %d, %s)\n", D, if (inference == "map") "Laplace" else sampler))
    choice <- lsdynam_choice_stage(events_df, n, D, sl, priors = priors, inference = inference,
                                   sampler = sampler, n_iter = n_iter, burn_in = burn_in, thin = thin,
                                   n_starts = n_starts, spectral_start = spectral_start, maxit = maxit,
                                   seed = seed, verbose = verbose)
    if (verbose && inference == "mcmc") cat(sprintf("  acceptance %.3f, %.1f s\n", choice$accept_rate, choice$seconds))
  }
  if (verbose) cat(sprintf("Stage 2 (rate module, link = %s)\n", paste(link, collapse = " + ")))
  rate <- lsdynam_rate_stage(events_df, n, sl, choice, link = link, link_type = link_type, df = df,
                             winsorize_q = winsorize_q, priors = priors, sigma_a = sigma_a, nu = nu,
                             inference = inference, sampler = sampler, n_iter = n_iter_rate,
                             burn_in = burn_in_rate, thin = thin_rate, n_cut = n_cut,
                             cut_refit_sigma = cut_refit_sigma, rate_event = rate_event,
                             plugin = plugin, scale_basis = scale_basis, sigma_kappa = sigma_kappa, seed = seed, verbose = verbose)
  if (verbose && inference == "mcmc") cat(sprintf("  acceptance %.3f, %.1f s\n", rate$accept_rate, rate$seconds))
  ps <- if (choice$D > 0) position_summaries(choice$U_mean, choice$V_mean) else NULL
  fit <- list(choice = choice, rate = rate, n = n, D = choice$D, stats = stats,
              stats_list = if (keep_stats) sl else NULL,
              summaries = ps,
              N_out = tabulate(events_df$sender, nbins = n), N_in = tabulate(events_df$receiver, nbins = n),
              M = nrow(events_df), priors = priors,
              map_fit = list(choice = choice$map, rate = rate$map))
  class(fit) <- c("lsdynam", "list")
  fit
}

#' @rdname summary.lsdynam
#' @export
print.lsdynam <- function(x, ...) {
  q <- function(v) sprintf("%7.3f [%7.3f, %7.3f]", mean(v), quantile(v, 0.025), quantile(v, 0.975))
  lk <- if (identical(x$rate$link$vars, "none")) "none" else
    sprintf("%s (%s), sigma_kappa = %.3f (%s)", paste(x$rate$link$vars, collapse = " + "), x$rate$link$type,
            x$rate$sigma_kappa, x$rate$sigma_kappa_method)
  cat(sprintf("Latent space DyNAM: n = %d, M = %d, D = %d, link = %s, sigma_a = %.3f (%s)\n",
              x$n, x$M, x$D, lk, x$rate$sigma_a, x$rate$sigma_a_method))
  cat(sprintf("Inference: choice %s, rate %s\n", x$choice$sampler, x$rate$sampler))
  for (nm in colnames(x$rate$xi))     cat(sprintf("  xi   %-16s %s\n", nm, q(x$rate$xi[, nm])))
  for (nm in colnames(x$choice$zeta)) cat(sprintf("  zeta %-16s %s\n", nm, q(x$choice$zeta[, nm])))
  cat(sprintf("  mu_a %-16s %s\n", "", q(x$rate$mu_a)))
  for (nm in colnames(x$rate$kappa_vec)) cat(sprintf("  kappa %-15s %s\n", nm, q(x$rate$kappa_vec[, nm])))
  invisible(x)
}

#' Summarise a fitted latent space DyNAM
#'
#' Posterior summaries of the coefficients and, for a model with a link, the
#' contrasts of its components, the proportion of the between-actor variation
#' of the activity effects that it explains, and the evidence for the whole
#' link and for each component given the others.
#'
#' @param object,x an object of class `lsdynam` (or `summary.lsdynam` for the
#'   print method).
#' @param level probability of the credible intervals.
#' @param evidence compute the Bayes factors with [link_evidence()].
#' @param digits number of significant digits to print.
#' @param ... not used.
#'
#' @return `summary()` returns an object of class `summary.lsdynam`: a list
#'   with the data frame `coefficients` (posterior mean, standard deviation and
#'   credible interval of \eqn{\xi}, \eqn{\zeta}, \eqn{\mu_a} and the link
#'   coefficients), `sigma_a`, the data frame `contrasts` (see
#'   [link_contrast()]), `R2` (see [link_r2()]), the data frame `evidence`
#'   (log Bayes factors by the Savage-Dickey ratio and by Laplace
#'   approximations) and information on the model. The print methods return
#'   their argument invisibly.
#' @seealso [lsdynam_fit()], [link_evidence()].
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' summary(fit)
#' @export
summary.lsdynam <- function(object, level = 0.95, evidence = TRUE, ...) {
  pr <- c((1 - level) / 2, 1 - (1 - level) / 2)
  one <- function(nm, d) data.frame(parameter = nm, mean = mean(d), sd = sd(d),
                                    lower = quantile(d, pr[1], names = FALSE),
                                    upper = quantile(d, pr[2], names = FALSE), stringsAsFactors = FALSE)
  rows <- list()
  for (nm in colnames(object$rate$xi)) rows[[length(rows) + 1]] <- one(paste0("xi[", nm, "]"), object$rate$xi[, nm])
  for (nm in colnames(object$choice$zeta)) rows[[length(rows) + 1]] <- one(paste0("zeta[", nm, "]"), object$choice$zeta[, nm])
  rows[[length(rows) + 1]] <- one("mu_a", object$rate$mu_a)
  for (nm in colnames(object$rate$kappa_vec)) rows[[length(rows) + 1]] <- one(paste0("kappa[", nm, "]"), object$rate$kappa_vec[, nm])
  coefs <- do.call(rbind, rows)
  vars <- object$rate$link$vars
  has_link <- !identical(vars, "none")
  ct <- NULL; ev <- NULL; r2 <- NA_real_
  if (has_link) {
    r2 <- link_r2(object)
    ct <- do.call(rbind, lapply(seq_along(vars), function(j) {
      cj <- link_contrast(object, component = j, level = level)
      data.frame(component = vars[j], x_10 = cj$x[1], x_90 = cj$x[2], mean = cj$mean, lower = cj$ci[1],
                 upper = cj$ci[2], prob_positive = cj$prob_positive, stringsAsFactors = FALSE)
    }))
    if (evidence) {
      e <- link_evidence(object)
      ev <- data.frame(test = "link", log_bf_savage_dickey = e$log_bf10_sd, log_bf_laplace = e$log_bf10_laplace,
                       stringsAsFactors = FALSE)
      if (length(vars) > 1) for (j in seq_along(vars)) {
        e <- link_evidence(object, components = j)
        ev <- rbind(ev, data.frame(test = paste(vars[j], "|", paste(vars[-j], collapse = " + ")),
                                   log_bf_savage_dickey = e$log_bf10_sd, log_bf_laplace = e$log_bf10_laplace,
                                   stringsAsFactors = FALSE))
      }
    }
  }
  out <- list(coefficients = coefs, sigma_a = object$rate$sigma_a, sigma_a_method = object$rate$sigma_a_method,
              sigma_kappa = object$rate$sigma_kappa, sigma_kappa_method = object$rate$sigma_kappa_method,
              link = vars, link_type = object$rate$link$type, contrasts = ct, R2 = r2, evidence = ev,
              n = object$n, M = object$M, D = object$D,
              inference = c(choice = object$choice$sampler, rate = object$rate$sampler),
              accept = c(choice = object$choice$accept_rate, rate = object$rate$accept_rate), level = level)
  class(out) <- "summary.lsdynam"
  out
}

#' @rdname summary.lsdynam
#' @export
print.summary.lsdynam <- function(x, digits = 3, ...) {
  cat(sprintf("Latent space DyNAM: n = %d actors, M = %d events, D = %d\n", x$n, x$M, x$D))
  cat(sprintf("Inference: choice module %s, rate module %s", x$inference[1], x$inference[2]))
  if (!all(is.na(x$accept))) cat(sprintf(" (acceptance %s)", paste(sprintf("%.2f", x$accept), collapse = ", ")))
  cat("\n")
  lk <- if (identical(x$link, "none")) "none" else paste0(paste(x$link, collapse = " + "), " (", x$link_type, ")")
  cat(sprintf("Link: %s; sigma_a = %.3f (%s)", lk, x$sigma_a, x$sigma_a_method))
  if (!identical(x$link, "none")) cat(sprintf(", sigma_kappa = %.3f (%s)", x$sigma_kappa, x$sigma_kappa_method))
  cat("\n\nCoefficients (posterior mean, sd and ", 100 * x$level, "% credible interval):\n", sep = "")
  cf <- x$coefficients
  rownames(cf) <- cf$parameter
  print(signif(as.matrix(cf[, -1]), digits))
  if (!is.null(x$contrasts)) {
    cat("\nContrasts of the components (level at the 90% minus level at the 10% quantile):\n")
    ct <- x$contrasts; rownames(ct) <- ct$component
    print(signif(as.matrix(ct[, -1]), digits))
    cat(sprintf("R2 of the link: %.3f\n", x$R2))
  }
  if (!is.null(x$evidence)) {
    cat("\nLog Bayes factors:\n")
    ev <- x$evidence; rownames(ev) <- ev$test
    print(signif(as.matrix(ev[, -1]), digits))
  }
  invisible(x)
}

#' Plot a fitted latent space DyNAM
#'
#' Plots the components of the link (one panel per component, see
#' [plot_link_level()]), the latent positions for \eqn{D = 1}
#' ([plot_positions()]), the activity effects ([plot_activity_effects()]) or
#' trace plots ([plot_traces()]).
#'
#' @param x an object of class `lsdynam`.
#' @param which what to plot: `"link"`, `"positions"`, `"activity"` or
#'   `"traces"`.
#' @param ... further arguments for the plotting function.
#' @return `x`, invisibly.
#' @examples
#' tr <- make_truth(n = 15, D = 1, link_fun = list(function(x) x, NULL),
#'                  link_var = c("pop", "rho"), mean_u = 1, seed = 1)
#' ev <- simulate_lsdynam(tr$theta, n = 15, M = 500, seed = 2)
#' fit <- lsdynam_fit(ev, n = 15, inference = "map", verbose = FALSE)
#' plot(fit)
#' plot(fit, which = "positions")
#' @export
plot.lsdynam <- function(x, which = c("link", "positions", "activity", "traces"), ...) {
  which <- match.arg(which)
  if (which == "link") {
    vars <- x$rate$link$vars
    if (identical(vars, "none")) stop("this model has no link")
    nc <- length(vars)
    if (nc > 1) { op <- par(mfrow = c(1, nc)); on.exit(par(op)) }
    for (j in seq_len(nc)) {
      plot_link_level(x, component = j, xlab = summary_label(vars[j]),
                      main = paste("Component in", summary_label(vars[j])),
                      legend_pos = if (j == 1) "topleft" else NULL, ...)
    }
  } else if (which == "positions") {
    plot_positions(x, ...)
  } else if (which == "activity") {
    plot_activity_effects(x, ...)
  } else {
    plot_traces(x, ...)
  }
  invisible(x)
}
