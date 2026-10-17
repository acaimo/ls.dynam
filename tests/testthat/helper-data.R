# Small simulated data set shared by the tests
sim_data <- function(n = 15, M = 400, D = 1, seed = 1) {
  f <- function(x) 1.5 * (plogis(3 * x) - 0.5)
  tr <- make_truth(n = n, D = D, link_fun = list(f, NULL), link_var = c("pop", "rho"),
                   mean_u = c(1, rep(0, D - 1)), seed = seed)
  ev <- simulate_lsdynam(tr$theta, n = n, M = M, seed = seed + 1)
  list(truth = tr, events = ev, n = n)
}
