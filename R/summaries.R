#' Invariant summaries of the latent positions
#'
#' The choice model identifies the latent positions only up to the
#' transformations \eqn{u_i \mapsto A u_i}{u_i -> A u_i},
#' \eqn{v_j \mapsto A^{-\top}(v_j + w)}{v_j -> A^{-T}(v_j + w)} for any
#' invertible \eqn{A} and any \eqn{w}. This function computes summaries of
#' each actor's positions that are invariant to these transformations, and
#' so are identified.
#'
#' With \eqn{d_i = v_i - \bar v}{d_i = v_i - vbar} the centred receiver
#' position, \eqn{\bar u}{ubar} the mean sender position and
#' \eqn{\Sigma_v = n^{-1}\sum_j d_j d_j^\top}{Sigma_v = sum_j d_j d_j' / n}:
#' * `r`: self-affinity \eqn{r_i = u_i^\top d_i}{r_i = u_i' d_i}, the latent
#'   log-odds that \eqn{i} chooses an alter at its own receiver position
#'   rather than one at the centroid;
#' * `pop`: latent popularity \eqn{\pi_i = \bar u^\top d_i}{pi_i = ubar' d_i},
#'   how much the average sender favours \eqn{i} as a receiver;
#' * `rho`: excess self-affinity \eqn{\rho_i = (u_i - \bar u)^\top d_i}{rho_i = (u_i - ubar)' d_i},
#'   so that \eqn{r_i = \pi_i + \rho_i};
#' * `phi`: selectivity \eqn{\phi_i = (u_i^\top \Sigma_v u_i)^{1/2}}{phi_i = (u_i' Sigma_v u_i)^(1/2)};
#' * `delta`: distinctiveness \eqn{\delta_i = (d_i^\top \Sigma_v^{-1} d_i)^{1/2}}{delta_i = (d_i' Sigma_v^(-1) d_i)^(1/2)};
#' * `align`: alignment \eqn{c_i = r_i / (\phi_i \delta_i)}{c_i = r_i / (phi_i delta_i)};
#' * `z`: standardised self-affinity \eqn{r_i / \phi_i}{r_i / phi_i};
#' * `log_phi`, `log_delta`: logarithms of `phi` and `delta`;
#' * `norm_u` (\eqn{\|u_i\|^2}{||u_i||^2}) and `uv` (\eqn{u_i^\top v_i}{u_i' v_i}),
#'   which are not invariant and are returned for comparison only.
#'
#' @param U,V \eqn{n \times D}{n x D} matrices of sender and receiver
#'   positions.
#'
#' @return A data frame with one row per actor and the columns described
#'   above (`NULL` when \eqn{D = 0}).
#' @seealso [summaries_draws()] for posterior draws, [summary_label()].
#' @examples
#' set.seed(1)
#' U <- matrix(rnorm(20, 1), 10, 2); V <- matrix(rnorm(20), 10, 2)
#' ps <- position_summaries(U, V)
#' all.equal(ps$r, ps$pop + ps$rho)
#' # invariance to the transformations that leave the likelihood unchanged
#' A <- matrix(c(2, 1, 0, 1), 2, 2); w <- c(3, -1)
#' ps2 <- position_summaries(U %*% t(A), sweep(V, 2, w, "+") %*% solve(A))
#' all.equal(ps[, c("r", "pop", "rho", "phi", "delta", "z")],
#'           ps2[, c("r", "pop", "rho", "phi", "delta", "z")])
#' @export
position_summaries <- function(U, V) {
  U <- as.matrix(U); V <- as.matrix(V)
  n <- nrow(U); D <- ncol(U)
  if (D == 0) return(NULL)
  d    <- sweep(V, 2, colMeans(V))            # receiver positions, centred
  ubar <- colMeans(U)
  r    <- rowSums(U * d)                      # self-affinity
  pop  <- as.numeric(d %*% ubar)              # latent popularity
  Sv   <- crossprod(d) / n
  phi  <- sqrt(pmax(rowSums((U %*% Sv) * U), 0))                # selectivity
  Sv_inv <- tryCatch(solve(Sv), error = function(e) NULL)
  delta <- if (is.null(Sv_inv)) rep(NA_real_, n) else
    sqrt(pmax(rowSums((d %*% Sv_inv) * d), 0))                  # distinctiveness
  z     <- ifelse(phi > 0, r / phi, 0)                          # standardised self-affinity
  align <- ifelse(phi * delta > 0, r / (phi * delta), 0)        # alignment
  data.frame(r = r, pop = pop, rho = r - pop, phi = phi, delta = delta, align = align, z = z,
             log_phi = log(pmax(phi, 1e-300)), log_delta = log(pmax(delta, 1e-300)),
             norm_u = rowSums(U^2), uv = rowSums(U * V))
}

#' Invariant summaries for posterior draws of the positions
#'
#' Applies [position_summaries()] to every draw (or to the draws `idx`) of
#' the latent positions.
#'
#' @param U_draws,V_draws arrays of dimension (draws, \eqn{n}, \eqn{D}), such
#'   as `fit$choice$U` and `fit$choice$V` of a fitted [lsdynam_fit()] object.
#' @param vars names of the summaries (columns of [position_summaries()]).
#' @param idx indices of the draws to use (default all).
#'
#' @return An array of dimension (draws, \eqn{n}, `length(vars)`).
#' @examples
#' set.seed(1)
#' Ud <- array(rnorm(3 * 10), c(3, 10, 1)); Vd <- array(rnorm(3 * 10), c(3, 10, 1))
#' dim(summaries_draws(Ud, Vd, vars = c("pop", "rho")))
#' @export
summaries_draws <- function(U_draws, V_draws, vars = "r", idx = NULL) {
  n <- dim(U_draws)[2]; D <- dim(U_draws)[3]
  if (is.null(idx)) idx <- seq_len(dim(U_draws)[1])
  out <- array(NA_real_, dim = c(length(idx), n, length(vars)), dimnames = list(NULL, NULL, vars))
  for (k in seq_along(idx)) {
    ps <- position_summaries(matrix(U_draws[idx[k], , ], n, D), matrix(V_draws[idx[k], , ], n, D))
    out[k, , ] <- as.matrix(ps[, vars, drop = FALSE])
  }
  out
}

#' Labels of the invariant summaries
#'
#' Full names of the summaries returned by [position_summaries()], for axis
#' labels and tables.
#'
#' @param v names of summaries, such as `"pop"` or `"rho"`.
#' @return A character vector of labels (unknown names are returned
#'   unchanged).
#' @examples
#' summary_label(c("pop", "rho", "r"))
#' @export
summary_label <- function(v) {
  lab <- c(pop = "latent popularity", rho = "excess self-affinity", r = "self-affinity",
           z = "standardised self-affinity", phi = "selectivity", delta = "distinctiveness",
           c = "alignment", align = "alignment", log_phi = "log selectivity",
           log_delta = "log distinctiveness", uv = "u'v (uncentred)", norm_u = "squared length of u")
  unname(ifelse(v %in% names(lab), lab[v], v))
}
