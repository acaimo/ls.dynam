test_that("compute_stats counts the history before each event", {
  ev <- data.frame(sender = c(1, 2, 1, 3), receiver = c(2, 1, 2, 1), time = c(1, 2, 4, 7))
  sl <- compute_stats(ev, n = 3)
  expect_equal(sl$M, 4)
  expect_equal(sl$delta_t, c(1, 1, 2, 3))
  expect_equal(sl$R, 2)
  expect_equal(sl$S, 4)
  # before event 3 (sender 1): actor 2 has sent one event to actor 1
  expect_equal(sl$choice_recip[3, ], log1p(c(0, 1, 0)))
  # before event 3: actor 1 has sent one event to actor 2
  expect_equal(sl$choice_inertia[3, ], log1p(c(0, 1, 0)))
  # before event 4: in-degrees (1, 2, 0) and out-degrees (2, 1, 0)
  expect_equal(sl$rate_indeg[4, ], log1p(c(1, 2, 0)))
  expect_equal(sl$rate_outdeg[4, ], log1p(c(2, 1, 0)))
  expect_equal(sl$choice_indeg[4, ], log1p(c(1, 2, 0)))
  expect_true(all(sl$rate_indeg[1, ] == 0))
})

test_that("compute_stats returns only the statistics requested", {
  ev <- data.frame(sender = c(1, 2), receiver = c(2, 1), time = c(1, 2))
  sl <- compute_stats(ev, n = 2, stats = c("rate_indeg", "choice_recip"))
  expect_setequal(grep("^(rate|choice)_", names(sl), value = TRUE), c("rate_indeg", "choice_recip"))
  expect_equal(c(sl$R, sl$S), c(1, 1))
})
