test_that("simulate_lsdynam returns valid events", {
  d <- sim_data(n = 10, M = 300)
  ev <- d$events
  expect_equal(nrow(ev), 300)
  expect_true(all(ev$sender != ev$receiver))
  expect_true(all(diff(ev$time) > 0))
  expect_true(all(ev$sender %in% 1:10) && all(ev$receiver %in% 1:10))
  tr <- d$truth
  expect_equal(dim(tr$theta$U), c(10, 1))
  expect_equal(length(tr$theta$a), 10)
  # the true component in popularity is centred over the actors
  expect_equal(mean(tr$component_fun$pop(tr$x[, "pop"])), 0, tolerance = 1e-12)
  expect_true(all(tr$component_fun$rho(tr$x[, "rho"]) == 0))
})

test_that("simulate_lsdynam stops at T_max and is reproducible", {
  tr <- make_truth(n = 8, D = 2, seed = 3)
  ev <- simulate_lsdynam(tr$theta, n = 8, T_max = 2, seed = 4)
  expect_true(all(ev$time <= 2))
  expect_identical(ev, simulate_lsdynam(tr$theta, n = 8, T_max = 2, seed = 4))
})
