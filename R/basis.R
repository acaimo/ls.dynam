# Basis of the position-rate link. X is the n x p matrix of link variables
# (one column per additive component). Each column is winsorised (optional,
# two-sided), centred, expanded in a natural cubic spline basis with df
# functions (interior knots at quantiles, boundary knots at the extremes) or
# kept linear; the basis columns are centred across actors, so that mu_a is
# the average level, and with scale = TRUE divided by their standard
# deviation, so that the prior kappa_k ~ N(0, sigma2_kappa) means the same
# thing for every link variable: each basis function contributes a term with
# prior standard deviation sigma_kappa to the activity effect of a typical
# actor. The returned spec evaluates the same basis at new values (grids,
# posterior draws of the summaries).
build_link_basis <- function(X, type = c("spline", "linear"), df = 4, winsorize_q = 1, scale = TRUE) {
  type <- match.arg(type)
  X <- as.matrix(X); p <- ncol(X); n <- nrow(X)
  vars <- colnames(X)
  if (is.null(vars)) vars <- paste0("x", seq_len(p))
  comps <- vector("list", p); Bl <- vector("list", p)
  for (j in seq_len(p)) {
    x <- X[, j]
    lim <- if (winsorize_q < 1) quantile(x, c(1 - winsorize_q, winsorize_q), names = FALSE) else c(-Inf, Inf)
    xw <- pmin(pmax(x, lim[1]), lim[2])
    center <- mean(xw)
    xc <- xw - center
    knots <- NULL; bk <- range(xc)
    if (type == "linear" || df <= 1) {
      Bj <- matrix(xc, ncol = 1)
      type_j <- "linear"
    } else {
      probs <- seq(0, 1, length.out = df + 1)[-c(1, df + 1)]
      knots <- unique(quantile(xc, probs, names = FALSE))
      knots <- knots[knots > bk[1] & knots < bk[2]]
      if (length(knots) < df - 1) {
        warning(sprintf("build_link_basis: '%s' has tied quantiles; %d interior knot(s) instead of %d. Many actors share the same value (typically actors with few sent events, shrunk to the prior).",
                        vars[j], length(knots), df - 1))
      }
      Bj <- matrix(splines::ns(xc, knots = knots, Boundary.knots = bk), nrow = n)
      type_j <- "spline"
    }
        cm <- colMeans(Bj)
    Bj <- sweep(Bj, 2, cm)
    cs <- if (scale) apply(Bj, 2, sd) else rep(1, ncol(Bj))
    cs[!is.finite(cs) | cs < 1e-12] <- 1
    Bj <- sweep(Bj, 2, cs, "/")
    colnames(Bj) <- paste0(vars[j], "_", seq_len(ncol(Bj)))
    comps[[j]] <- list(var = vars[j], type = type_j, center = center, limits = lim, knots = knots,
                       boundary = bk, col_means = cm, col_sds = cs, K = ncol(Bj))
    Bl[[j]] <- Bj
  }
  B <- do.call(cbind, Bl)
  Ks <- vapply(comps, function(cc) cc$K, numeric(1))
  cols <- split(seq_len(ncol(B)), rep(seq_len(p), Ks))
  list(B = B, K = ncol(B), spec = list(comps = comps, cols = cols, vars = vars))
}

# Basis of one component at raw (uncentred) values x.
eval_link_component <- function(comp, x) {
  xw <- pmin(pmax(x, comp$limits[1]), comp$limits[2])
  xc <- xw - comp$center
  Bj <- if (comp$type == "linear") matrix(xc, ncol = 1) else
        matrix(splines::ns(xc, knots = comp$knots, Boundary.knots = comp$boundary), nrow = length(xc))
  sweep(sweep(Bj, 2, comp$col_means), 2, comp$col_sds, "/")
}

eval_link_basis <- function(spec, X) {
  X <- as.matrix(X)
  B <- do.call(cbind, lapply(seq_along(spec$comps), function(j) eval_link_component(spec$comps[[j]], X[, j])))
  colnames(B) <- unlist(lapply(spec$comps, function(cc) paste0(cc$var, "_", seq_len(cc$K))))
  B
}
