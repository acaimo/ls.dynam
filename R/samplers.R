# Metropolis-adjusted Langevin algorithm, preconditioned with a fixed
# covariance (the Laplace covariance at the mode). fn(par) must return
# list(value = log-posterior, gradient = its gradient). The step size is
# adapted during burn-in only (target acceptance 0.57) and then frozen, so
# the retained draws come from a valid Metropolis-Hastings chain. With d
# parameters a random walk needs of the order of d iterations per independent
# draw, MALA of the order of d^(1/3).
# align (optional, used for the choice module): a function of par returning
# list(par = the configuration rotated to the reference orientation, rot = a
# function applying the same rotation to another configuration). The
# posterior is invariant to rotations of the latent space, but the
# preconditioner is only appropriate near the orientation of the mode; with
# align the kernel commutes with rotations, the Hastings ratio evaluates the
# reverse proposal in the frame of the proposed configuration and the state
# is kept in the reference orientation.
mala_sample <- function(fn, par0, cov0, n_iter = 5000, burn_in = 1000, thin = 1,
                        eps = NULL, target_accept = 0.57, align = NULL, seed = 1, verbose = FALSE) {
  set.seed(seed)
  d <- length(par0)
  L <- t(chol(make_pd(cov0)))                       # cov0 = L L'
  if (is.null(eps)) eps <- 1.2 * d^(-1 / 6)
  cur <- if (is.null(align)) par0 else align(par0)$par
  f <- fn(cur)
  if (!is.finite(f$value)) stop("mala_sample: non-finite log-posterior at the starting point")
  g <- as.numeric(crossprod(L, f$gradient))         # gradient in whitened coordinates
  keep_it <- seq(burn_in + 1, n_iter, by = thin)
  samples <- matrix(NA_real_, length(keep_it), d)
  lp <- numeric(length(keep_it))
  n_acc <- 0; k <- 0
  for (it in seq_len(n_iter)) {
    z <- rnorm(d)
    step <- 0.5 * eps^2 * g + eps * z
    prop <- cur + as.numeric(L %*% step)
    if (is.null(align)) {
      back_w <- -step                                # L^{-1} (cur - prop)
    } else {
      al <- align(prop)
      back_w <- forwardsolve(L, al$rot(cur) - al$par)
      prop <- al$par
    }
    fp <- fn(prop)
    acc_prob <- 0
    if (is.finite(fp$value)) {
      gp <- as.numeric(crossprod(L, fp$gradient))
      back <- back_w - 0.5 * eps^2 * gp
      log_alpha <- fp$value - f$value - sum(back^2) / (2 * eps^2) + sum(z^2) / 2
      if (is.finite(log_alpha)) acc_prob <- min(1, exp(log_alpha))
    }
    if (runif(1) < acc_prob) {
      cur <- prop; f <- fp; g <- gp
      if (it > burn_in) n_acc <- n_acc + 1
    }
    if (it <= burn_in) eps <- eps * exp((acc_prob - target_accept) * it^(-0.6))
    if (it > burn_in && (it - burn_in - 1) %% thin == 0) {
      k <- k + 1
      samples[k, ] <- cur; lp[k] <- f$value
    }
    if (verbose && it %% 1000 == 0) cat(sprintf("  MALA iteration %d / %d (eps = %.3f)\n", it, n_iter, eps))
  }
  list(samples = samples, accept_rate = n_acc / max(1, n_iter - burn_in), eps = eps, log_post = lp)
}

# Random-walk Metropolis with a fixed proposal covariance (sampler = "rwm").
rwm_sample <- function(log_post, par0, cov0, n_iter = 8000, burn_in = 2000, thin = 5,
                        scale = 2.38 / sqrt(length(par0)), seed = 1) {
  set.seed(seed)
  d <- length(par0)
  L <- chol(make_pd(cov0))
  cur <- par0
  cur_lp <- log_post(cur)
  keep_every <- thin
  n_keep <- length(seq(burn_in + 1, n_iter, by = keep_every))
  samples <- matrix(NA_real_, n_keep, d)
  n_accept <- 0
  k <- 0
  for (it in seq_len(n_iter)) {
    prop <- cur + scale * as.numeric(rnorm(d) %*% L)
    prop_lp <- log_post(prop)
    if (is.finite(prop_lp) && log(runif(1)) < (prop_lp - cur_lp)) {
      cur <- prop; cur_lp <- prop_lp; n_accept <- n_accept + 1
    }
    if (it > burn_in && (it - burn_in) %% keep_every == 0) {
      k <- k + 1
      samples[k, ] <- cur
    }
  }
  list(samples = samples[seq_len(k), , drop = FALSE], accept_rate = n_accept / n_iter)
}

# Symmetric positive definite version of M (eigenvalues floored at
# min_eig_ratio times the largest), as a guard for the Cholesky factorisation.
make_pd <- function(M, min_eig_ratio = 1e-10) {
  eig <- eigen((M + t(M)) / 2, symmetric = TRUE)
  floor_val <- min_eig_ratio * max(eig$values)
  vals <- pmax(eig$values, floor_val)
  eig$vectors %*% diag(vals, length(vals)) %*% t(eig$vectors)
}

# Covariance of the Laplace approximation from a (numerical) Hessian of the
# negative log-posterior. Eigenvalues below the floor are raised to the floor
# before inversion, so that flat directions get a large but finite variance.
hessian_to_cov <- function(H, abs_floor = 0, rel_floor = 1e-10) {
  H  <- (H + t(H)) / 2
  eg <- eigen(H, symmetric = TRUE)
  fl <- max(abs_floor, rel_floor * max(eg$values))
  vals <- pmax(eg$values, fl)
  list(cov = eg$vectors %*% (t(eg$vectors) / vals),
       n_floored = sum(eg$values < fl), eigenvalues = eg$values, floor = fl)
}
