inv_cols <- c("r", "pop", "rho", "phi", "delta", "align", "z")

test_that("summaries are invariant to the transformations of the positions", {
  set.seed(42)
  for (D in 1:3) {
    n <- 12
    U <- matrix(rnorm(n * D, 1), n, D); V <- matrix(rnorm(n * D), n, D)
    A <- matrix(rnorm(D * D), D, D) + diag(2, D); w <- rnorm(D)
    U2 <- U %*% t(A)                                 # u_i -> A u_i
    V2 <- sweep(V, 2, w, "+") %*% solve(A)           # v_j -> A^{-T}(v_j + w)
    # the latent affinities, and so the choice probabilities, are unchanged up to a sender term
    aff  <- U  %*% t(sweep(V,  2, colMeans(V)))
    aff2 <- U2 %*% t(sweep(V2, 2, colMeans(V2)))
    expect_equal(aff, aff2, tolerance = 1e-10)
    ps <- position_summaries(U, V); ps2 <- position_summaries(U2, V2)
    expect_equal(ps[, inv_cols], ps2[, inv_cols], tolerance = 1e-8)
    expect_equal(ps$r, ps$pop + ps$rho)
  }
})

test_that("a change of sign of all positions leaves the summaries unchanged for D = 1", {
  set.seed(1)
  U <- matrix(rnorm(8), 8, 1); V <- matrix(rnorm(8), 8, 1)
  expect_equal(position_summaries(U, V)[, inv_cols], position_summaries(-U, -V)[, inv_cols])
})

test_that("summaries_draws agrees with position_summaries", {
  set.seed(2)
  Ud <- array(rnorm(4 * 10 * 2), c(4, 10, 2)); Vd <- array(rnorm(4 * 10 * 2), c(4, 10, 2))
  sd <- summaries_draws(Ud, Vd, vars = c("pop", "rho"))
  expect_equal(dim(sd), c(4, 10, 2))
  expect_equal(sd[3, , "pop"], position_summaries(Ud[3, , ], Vd[3, , ])$pop)
  expect_equal(summary_label(c("pop", "rho", "unknown")), c("latent popularity", "excess self-affinity", "unknown"))
})

test_that("position_summaries returns NULL without latent dimensions", {
  expect_null(position_summaries(matrix(0, 5, 0), matrix(0, 5, 0)))
})
