test_that("lsdynam_fit with the Laplace approximation", {
  d <- sim_data(n = 15, M = 500)
  fit <- lsdynam_fit(d$events, d$n, D = 1, inference = "map", n_starts = 2, verbose = FALSE)
  expect_s3_class(fit, "lsdynam")
  expect_equal(fit$rate$link$vars, c("pop", "rho"))
  expect_equal(dim(fit$choice$U)[2:3], c(15, 1))
  expect_equal(ncol(fit$rate$a), 15)
  expect_output(print(fit), "Latent space DyNAM")
  s <- summary(fit)
  expect_s3_class(s, "summary.lsdynam")
  expect_output(print(s), "Log Bayes factors")
  expect_equal(nrow(s$evidence), 3)
  # components by name or by position
  expect_equal(link_contrast(fit, "rho")$mean, link_contrast(fit, 2)$mean)
  expect_error(link_contrast(fit, "phi"), "unknown component")
  e <- link_evidence(fit, components = "rho", sigma_kappa_grid = c(0.5, 1))
  expect_true(is.finite(e$log_bf10_laplace) && is.finite(e$log_bf10_sd))
  expect_length(e$grid, 2)
  expect_true(link_r2(fit) > -1 && link_r2(fit) < 1)
  # the same positions without link
  fit1 <- lsdynam_fit(d$events, d$n, link = "none", inference = "map", choice = fit$choice, verbose = FALSE)
  expect_identical(fit1$choice$U_mean, fit$choice$U_mean)
  expect_equal(link_r2(fit1), 0)
  tab <- link_comparison_table(fit, list(both = c("pop", "rho"), pop = "pop", r = "r"))$table
  expect_equal(tab$link, c("both", "pop", "r"))
  ms <- mode_stability_check(fit, plot = FALSE, verbose = FALSE)
  expect_true(nrow(ms$table) >= 2)
  expect_equal(min(ms$table$delta), 0)
})

test_that("the model without latent positions", {
  d <- sim_data(n = 10, M = 200)
  f0 <- lsdynam_fit(d$events, d$n, D = 0, inference = "map", verbose = FALSE)
  expect_equal(f0$D, 0)
  expect_identical(f0$rate$link$vars, "none")
  expect_null(f0$summaries)
  expect_error(lsdynam_fit(d$events, d$n, D = 0, link = "pop", inference = "map", verbose = FALSE),
               "needs latent positions")
})

test_that("lsdynam_fit with MALA and the cut posterior", {
  d <- sim_data(n = 10, M = 250)
  fit <- lsdynam_fit(d$events, d$n, D = 1, n_iter = 400, burn_in = 200, thin = 2, n_cut = 20,
                     verbose = FALSE)
  expect_equal(nrow(fit$choice$zeta), 100)
  expect_true(fit$choice$accept_rate > 0 && fit$choice$accept_rate < 1)
  expect_true(fit$rate$accept_rate > 0 && fit$rate$accept_rate < 1)
  expect_equal(nrow(fit$rate$cut$kappa_vec), 20)
  lc <- link_curve(fit, "pop", source = "cut", n_grid = 7)
  expect_equal(dim(lc$level_ci), c(7, 2))
  dg <- mcmc_diagnostics(fit, verbose = FALSE)
  expect_true(all(c("zeta", "pop", "rho", "a") %in% names(dg$ess)))
})

test_that("predictive scores and their differences", {
  d <- sim_data(n = 12, M = 400)
  m_tr <- 320
  f0 <- lsdynam_fit(d$events[1:m_tr, ], d$n, D = 0, inference = "map", verbose = FALSE)
  f1 <- lsdynam_fit(d$events[1:m_tr, ], d$n, D = 1, link = "none", inference = "map", verbose = FALSE)
  s0 <- predictive_scores(f0, d$events, d$n, m_train = m_tr, n_draws = 50)
  s1 <- predictive_scores(f1, d$events, d$n, m_train = m_tr, n_draws = 50)
  expect_equal(s1$n_test, 80)
  expect_length(s1$choice_by_event, 80)
  expect_length(s1$rate_by_actor, 12)
  expect_true(is.finite(s1$choice) && is.finite(s1$rate))
  # with the same rate module the rate scores agree
  expect_equal(s0$rate, s1$rate, tolerance = 0.5)
  sd <- score_difference(s1, s0, groups = factor(s1$N_out_train > median(s1$N_out_train)))
  expect_equal(sd$part[1:2], c("rate", "choice"))
  expect_equal(nrow(sd), 4)
})
