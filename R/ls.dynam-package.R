#' ls.dynam: latent space dynamic network actor models with a position-rate link
#'
#' Bayesian inference for latent space dynamic network actor models (DyNAMs)
#' of relational event data, in which a position-rate link relates invariant
#' summaries of the latent positions to the activity of each actor.
#'
#' @section Model:
#' An actor \eqn{i} becomes active with intensity
#' \deqn{\lambda_i(t) = \exp\{\xi^\top s(i;t) + a_i\}}{lambda_i(t) = exp{xi' s(i;t) + a_i}}
#' and, when active, selects receiver \eqn{j \ne i}{j != i} with probability
#' proportional to \eqn{\exp\{\zeta^\top s(i,j;t) + u_i^\top v_j\}}{exp{zeta' s(i,j;t) + u_i' v_j}},
#' where \eqn{u_i} and \eqn{v_j} are latent sender and receiver positions in
#' \eqn{D} dimensions. The position-rate link is the hierarchical prior
#' \deqn{a_i \sim N(\mu_a + g_\pi(\pi_i) + g_\rho(\rho_i), \sigma_a^2),}{a_i ~ N(mu_a + g_pi(pi_i) + g_rho(rho_i), sigma_a^2),}
#' where \eqn{\pi_i} (latent popularity) and \eqn{\rho_i} (excess
#' self-affinity) are summaries of the positions that the choice model
#' identifies (see [position_summaries()]) and the components \eqn{g} are
#' natural cubic splines.
#'
#' @section Inference:
#' Inference is modular. Stage 1 estimates the choice module (positions and
#' choice coefficients) from the choice of receivers alone; stage 2 estimates
#' the rate module given summaries of the positions, either at plug-in values
#' (the plug-in posterior) or propagating their uncertainty (the cut
#' posterior). Both stages use a Metropolis-adjusted Langevin algorithm
#' preconditioned with the Laplace approximation at the posterior mode, or the
#' Laplace approximation alone (`inference = "map"`).
#'
#' @section Main functions:
#' * Data: [read_konect_events()], [prepare_events()], [select_time_window()],
#'   [filter_active_both()], [compute_stats()].
#' * Fitting: [lsdynam_fit()], [lsdynam_choice_stage()],
#'   [lsdynam_rate_stage()], [lsdynam_priors()].
#' * Summaries of the positions: [position_summaries()], [summaries_draws()].
#' * The link: [link_curve()], [link_contrast()], [link_evidence()],
#'   [link_r2()], [link_comparison_table()].
#' * Plots: [plot.lsdynam()], [plot_link_level()], [plot_link_slope()],
#'   [plot_positions()], [plot_activity_effects()], [plot_traces()].
#' * Prediction: [predictive_scores()], [score_difference()].
#' * Diagnostics and calibration: [mcmc_diagnostics()],
#'   [mode_stability_check()], [null_link_calibration()],
#'   [link_shrinkage_check()], [link_leverage_check()],
#'   [summary_agreement_check()].
#' * Simulation: [make_truth()], [simulate_lsdynam()].
#'
#' @references
#' Caimo, A. and Gollini, I. Latent space dynamic network actor models with
#' position-rate link. Manuscript.
#'
#' Stadtfeld, C. and Block, P. (2017). Interactions, actors, and time:
#' dynamic network actor models for relational events. *Sociological
#' Science*, 4, 318--352.
#'
#' Hoff, P. D. (2021). Additive and multiplicative effects network models.
#' *Statistical Science*, 36, 34--50.
#'
#' Plummer, M. (2015). Cuts in Bayesian graphical models. *Statistics and
#' Computing*, 25, 37--43.
#'
#' Roberts, G. O. and Tweedie, R. L. (1996). Exponential convergence of
#' Langevin distributions and their discrete approximations. *Bernoulli*, 2,
#' 341--363.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats acf cor cor.test cov median optim optimHess optimize rexp
#'   qchisq quantile rnorm runif sd setNames var
#' @importFrom graphics abline axis legend lines matplot par plot points
#'   polygon rug segments text
NULL
