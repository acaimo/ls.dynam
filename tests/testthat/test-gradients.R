num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(k) {
    e <- replace(numeric(length(x)), k, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(1))
}

test_that("the compressed rate likelihood equals the direct computation", {
  d <- sim_data(n = 8, M = 150)
  sl <- compute_stats(d$events, d$n, stats = c("rate_indeg", "rate_outdeg"))
  cr <- ls.dynam:::compress_rate_stats(sl, d$events$sender, d$n)
  set.seed(5)
  xi <- c(0.3, -0.2); a <- rnorm(d$n, -1, 0.5)
  lin <- xi[1] * sl[[cr$rnames[1]]] + xi[2] * sl[[cr$rnames[2]]] + matrix(a, sl$M, d$n, byrow = TRUE)
  direct <- sum(lin[cbind(seq_len(sl$M), d$events$sender)]) - sum(sl$delta_t * rowSums(exp(lin)))
  par <- c(xi, a, 0)
  ev <- ls.dynam:::rate_mod_eval(par, cr, matrix(0, d$n, 0), lsdynam_priors(), sigma2_a = 1, deriv = 0)
  expect_equal(ev$loglik, direct, tolerance = 1e-10)
})

test_that("gradient and Hessian of the rate module are correct", {
  d <- sim_data(n = 7, M = 120)
  sl <- compute_stats(d$events, d$n, stats = c("rate_indeg", "choice_recip"))
  cr <- ls.dynam:::compress_rate_stats(sl, d$events$sender, d$n)
  set.seed(6)
  B <- matrix(rnorm(d$n * 2), d$n, 2)
  hp <- lsdynam_priors()
  par <- c(0.1, rnorm(d$n, -1, 0.3), -0.8, 0.2, -0.1)
  for (nu in c(Inf, 4)) {
    f <- function(p) ls.dynam:::rate_mod_eval(p, cr, B, hp, sigma2_a = 0.5, nu = nu, deriv = 0)$value
    ev <- ls.dynam:::rate_mod_eval(par, cr, B, hp, sigma2_a = 0.5, nu = nu, deriv = 2)
    expect_equal(unname(ev$gradient), num_grad(f, par), tolerance = 1e-5)
    g <- function(p) ls.dynam:::rate_mod_eval(p, cr, B, hp, sigma2_a = 0.5, nu = nu, deriv = 1)$gradient
    Hn <- -sapply(seq_along(par), function(k) {
      e <- replace(numeric(length(par)), k, 1e-6); (g(par + e) - g(par - e)) / 2e-6
    })
    expect_equal(unname(ev$hessian), unname(Hn), tolerance = 1e-4)
  }
})

test_that("gradient of the choice module is correct", {
  d <- sim_data(n = 6, M = 100, D = 2)
  sl <- compute_stats(d$events, d$n, stats = c("choice_recip", "choice_inertia"))
  cn <- grep("^choice_", names(sl), value = TRUE)
  n <- d$n; D <- 2; S <- length(cn)
  set.seed(7)
  par <- rnorm(S + 2 * n * D, 0, 0.5)
  args <- list(n = n, D = D, S = S, sl = sl, senders = d$events$sender, receivers = d$events$receiver,
               cnames = cn, sigma2_zeta = 4, sigma2_u = 1, sigma2_v = 1)
  vg <- do.call(ls.dynam:::choice_val_grad, c(list(par = par), args))
  nll <- function(p) do.call(ls.dynam:::choice_negloglik, c(list(par = p), args))
  expect_equal(vg$value, -nll(par))
  expect_equal(vg$gradient, -num_grad(nll, par), tolerance = 1e-5)
  expect_equal(vg$gradient, -do.call(ls.dynam:::choice_neggrad, c(list(par = par), args)), tolerance = 1e-10)
})
