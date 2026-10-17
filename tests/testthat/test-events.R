test_that("prepare_events removes self-loops, identifies messages and separates ties", {
  ev <- data.frame(sender = c(1, 1, 2, 3, 3, 2), receiver = c(2, 3, 2, 1, 2, 1),
                   time = c(0, 0, 3600, 7200, 7200, 86400))
  out <- prepare_events(ev, verbose = FALSE)
  expect_equal(nrow(out), 5)
  expect_true(all(out$sender != out$receiver))
  expect_equal(out$msg_id, c(1, 1, 2, 2, 3))
  expect_true(all(diff(out$time) > 0))
  expect_equal(max(out$time), 1, tolerance = 1e-5)
})

test_that("select_time_window and exclude_actors", {
  ev <- data.frame(sender = c(1, 2, 3, 4, 1), receiver = c(2, 3, 4, 1, 3), time = 0:4)
  w <- select_time_window(ev, start_frac = 0.25, end_frac = 0.75, verbose = FALSE)
  expect_equal(w$time, c(0, 1, 2))
  ex <- exclude_actors(ev, n = 4, actors_to_exclude = 2, verbose = FALSE)
  expect_equal(ex$n, 3)
  expect_equal(nrow(ex$events_df), 3)
  expect_true(all(unlist(ex$events_df[, c("sender", "receiver")]) %in% 1:3))
})

test_that("filter_active_both repeats until every actor sends and receives", {
  ev <- data.frame(sender = c(1, 2, 3, 1, 4), receiver = c(2, 1, 1, 4, 3), time = 1:5)
  fa <- filter_active_both(ev, verbose = FALSE)
  s <- fa$events_df
  expect_true(all(seq_len(fa$n) %in% s$sender) && all(seq_len(fa$n) %in% s$receiver))
  ev2 <- data.frame(sender = c(1, 2, 3, 1), receiver = c(2, 1, 1, 4), time = 1:4)
  fa2 <- filter_active_both(ev2, verbose = FALSE)
  expect_equal(fa2$n, 2)
  expect_equal(fa2$old_id, c(1, 2))
  expect_equal(fa2$n_pass, 2)
})

test_that("read_konect_events reads the KONECT format", {
  f <- tempfile()
  writeLines(c("% asym unweighted", "% 3 3 3", "1 2 1 100", "2 1 1 160", "1 3 1 130"), f)
  ev <- read_konect_events(f)
  expect_equal(ev$time, c(100, 130, 160))
  expect_equal(ev$receiver, c(2, 3, 1))
})
