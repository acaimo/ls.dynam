## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 3.6,
                      out.width = "100%")

## ----simulate-----------------------------------------------------------------
library(ls.dynam)
f <- function(x) 1.5 * (plogis(3 * x) - 0.5)          # a smooth step
tr <- make_truth(n = 40, D = 1, link_fun = list(pop = f, rho = NULL),
                 link_var = c("pop", "rho"), mean_u = 1, sigma_a = 0.5,
                 xi = c(rate_indeg = 0.1), zeta = c(choice_recip = 0.5), seed = 1)
ev <- simulate_lsdynam(tr$theta, n = 40, M = 3000, seed = 2)
head(ev, 3)

## ----fit----------------------------------------------------------------------
fit <- lsdynam_fit(ev, n = 40, D = 1, inference = "map", n_starts = 3, verbose = FALSE)
summary(fit)

## ----link-plot----------------------------------------------------------------
plot(fit)

## ----link-summaries-----------------------------------------------------------
ct <- link_contrast(fit, "pop")
round(c(estimate = ct$mean, lower = ct$ci[1], upper = ct$ci[2]), 2)
link_evidence(fit, components = "rho")$log_bf10_laplace   # excess self-affinity given popularity

## ----comparison---------------------------------------------------------------
cmp <- link_comparison_table(fit, list("pop + rho" = c("pop", "rho"), pop = "pop",
                                       rho = "rho", r = "r"))
cmp$table[, c("link", "K", "sigma_a", "log_ml_vs_M1", "R2")]

## ----positions, fig.height = 4.5, fig.width = 5, out.width = "70%"------------
plot_positions(fit, ellipses = FALSE)

## ----prediction---------------------------------------------------------------
m_tr <- 2400
train <- ev[seq_len(m_tr), ]
m0 <- lsdynam_fit(train, n = 40, D = 0, inference = "map", verbose = FALSE)
m2 <- lsdynam_fit(train, n = 40, D = 1, inference = "map", verbose = FALSE)
m1 <- lsdynam_fit(train, n = 40, D = 1, link = "none", inference = "map",
                  choice = m2$choice, verbose = FALSE)
sc <- lapply(list(M0 = m0, M1 = m1, M2 = m2), predictive_scores,
             events_full = ev, n = 40, m_train = m_tr, n_draws = 100)
rbind("M1 - M0" = score_difference(sc$M1, sc$M0)[2, 2:4],
      "M2 - M1" = score_difference(sc$M2, sc$M1)[1, 2:4])

## ----mcmc, eval = FALSE-------------------------------------------------------
#  fit <- lsdynam_fit(ev, n = 40, D = 1, n_iter = 12000, burn_in = 2000, thin = 5,
#                     n_starts = 30, n_cut = 300)
#  mcmc_diagnostics(fit)                       # effective sample sizes, acceptance
#  plot_traces(fit)
#  mode_stability_check(fit)                   # the link at each mode of the choice module
#  plot(fit)                                   # plug-in (shaded) and cut (dashed) posteriors
#  link_evidence(fit, components = "rho", sigma_kappa_grid = c(0.1, 0.25, 0.5, 1, 2))
#  # reference distributions: no association, and popularity only
#  cal0 <- null_link_calibration(fit, n_rep = 20)
#  cal1 <- null_link_calibration(fit, n_rep = 20, keep = "pop")

## ----real, eval = FALSE-------------------------------------------------------
#  ev <- read_konect_events("out.radoslaw_email_email")
#  ev <- prepare_events(ev)                       # days; messages; ties separated
#  ev <- select_time_window(ev, start_frac = 0.2, end_frac = 0.4)
#  fa <- filter_active_both(ev)
#  fit <- lsdynam_fit(fa$events_df, n = fa$n, D = 1, n_iter = 12000, burn_in = 2000,
#                     thin = 5, n_starts = 30, n_cut = 300)

