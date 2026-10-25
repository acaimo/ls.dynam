# ls.dynam

Latent space dynamic network actor models (DyNAMs) with a position-rate link,
for relational event data such as e-mail.

In a DyNAM each event is the outcome of two decisions: when an actor becomes
active (the rate sub-process) and whom it selects (the choice sub-process).
**ls.dynam** adds latent sender and receiver positions to the choice
sub-process and a *position-rate link* that relates summaries of the
positions identified by the choice model, latent popularity and excess
self-affinity, to the activity of each actor:

* rate: `lambda_i(t) = exp{xi' s(i;t) + a_i}`
* choice: `p(j | i, t)` proportional to `exp{zeta' s(i,j;t) + u_i' v_j}`
* link: `a_i ~ N(mu_a + g_pi(pi_i) + g_rho(rho_i), sigma_a^2)`, with natural
  cubic splines `g`.

Inference is Bayesian and modular: the positions are estimated from the choice
of receivers alone, by a Metropolis-adjusted Langevin algorithm or a Laplace
approximation, and the link under the plug-in or the cut posterior.

## Installation

From the source file:

```r
install.packages("ls.dynam_0.1.0.tar.gz", repos = NULL, type = "source")
```

The package needs only base R (packages `stats`, `graphics`, `grDevices`,
`splines` and `utils`).

## Example

```r
library(ls.dynam)

# simulate events in which activity depends on latent popularity
f  <- function(x) 1.5 * (plogis(3 * x) - 0.5)
tr <- make_truth(n = 40, D = 1, link_fun = list(f, NULL), link_var = c("pop", "rho"),
                 mean_u = 1, seed = 1)
ev <- simulate_lsdynam(tr$theta, n = 40, M = 3000, seed = 2)

# fit (Laplace approximation; use inference = "mcmc", the default, for the
# full posterior)
fit <- lsdynam_fit(ev, n = 40, D = 1, inference = "map", n_starts = 3)
summary(fit)                                  # coefficients, contrasts, Bayes factors
plot(fit)                                     # components of the link
link_evidence(fit, components = "rho")        # excess self-affinity given popularity
```

See `vignette("ls.dynam")` for a longer introduction, with predictive
comparisons and the settings used in the paper. The scripts that reproduce
the analyses of the paper are in `system.file("scripts", package = "ls.dynam")`
(see the README there).

## Main functions

| Task | Functions |
|---|---|
| Data | `read_konect_events()`, `prepare_events()`, `select_time_window()`, `filter_active_both()`, `exclude_actors()`, `compute_stats()` |
| Fitting | `lsdynam_fit()`, `lsdynam_choice_stage()`, `lsdynam_rate_stage()`, `lsdynam_priors()` |
| Summaries of the positions | `position_summaries()`, `summaries_draws()`, `summary_label()` |
| The link | `summary()`, `link_curve()`, `link_contrast()`, `link_evidence()`, `link_r2()`, `link_comparison_table()` |
| Plots | `plot()`, `plot_link_level()`, `plot_link_slope()`, `plot_positions()`, `plot_activity_effects()`, `plot_traces()` |
| Prediction | `predictive_scores()`, `score_difference()` |
| Diagnostics and calibration | `mcmc_diagnostics()`, `ess()`, `mode_stability_check()`, `null_link_calibration()`, `link_shrinkage_check()`, `link_leverage_check()`, `summary_agreement_check()` |
| Simulation | `make_truth()`, `simulate_lsdynam()` |

Functions with a `seed` argument call `set.seed()`, so that results are
reproducible.

## Reference

Caimo, A. and Gollini, I. Latent space dynamic network actor models with
position-rate link. Manuscript.
