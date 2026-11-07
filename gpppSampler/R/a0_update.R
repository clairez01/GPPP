# Precompute all quantities needed for the O(n0) a0 log-posterior (Eq. 10).
# Must be called after each update to theta and sigma2.
#
# Returns a list with:
#   lambda  - eigenvalues of K_tilde = Kn0n0 - Kn0n (Knn + sigma2 I)^{-1} Knn0
#   d       - eigenvalues of Kn0n0
#   v       - Q^T (y0 - mu_{y0|y}),  Q = eigenvectors of K_tilde
#   w       - U^T y0,                 U = eigenvectors of Kn0n0
precompute_spectral <- function(y, y0, Km, sigma2, mu = 0) {
  yc  <- y  - mu
  y0c <- y0 - mu
  n  <- length(yc)
  ch_Knn <- safe_chol(Km$Knn + sigma2 * diag(n))
  if (is.null(ch_Knn)) return(NULL)

  # mu_{y0|y} = Kn0n (Knn + sigma2 I)^{-1} (y - mu)
  alpha_y <- backsolve(ch_Knn, forwardsolve(t(ch_Knn), yc))
  mu_y0_y <- as.vector(Km$Kn0n %*% alpha_y)

  # K_tilde = Kn0n0 - Kn0n (Knn + sigma2 I)^{-1} Knn0
  V       <- backsolve(ch_Knn, forwardsolve(t(ch_Knn), Km$Knn0))
  Kt      <- Km$Kn0n0 - Km$Kn0n %*% V
  Kt      <- 0.5 * (Kt + t(Kt))      # symmetrise against floating-point drift

  eig_Kt  <- eigen(Kt,          symmetric = TRUE)
  eig_K00 <- eigen(Km$Kn0n0,   symmetric = TRUE)

  list(
    lambda = pmax(eig_Kt$values,  0),
    d      = pmax(eig_K00$values, 0),
    v      = as.vector(t(eig_Kt$vectors)  %*% (y0c - mu_y0_y)),
    w      = as.vector(t(eig_K00$vectors) %*% y0c)
  )
}

# Eq. (10): O(n0) log-posterior for a0 (up to an additive constant in a0).
# nu = sigma2 / a0 is the historical pseudo-noise variance.
#
# log p(a0 | y, y0, sigma2) =
#   -1/2 * sum_i [ log(lambda_i + nu) + v_i^2 / (lambda_i + nu) ]   [numerator]
#   +1/2 * sum_i [ log(d_i     + nu) + w_i^2 / (d_i     + nu)   ]   [denominator]
log_post_a0 <- function(a0, sigma2, sp) {
  if (a0 <= 0 || a0 >= 1) return(-Inf)
  nu  <- sigma2 / a0
  -0.5 * sum(log(sp$lambda + nu) + sp$v^2 / (sp$lambda + nu)) +
   0.5 * sum(log(sp$d     + nu) + sp$w^2 / (sp$d     + nu))
}

# Griddy-Gibbs draw for a0 via discrete approximation on a uniform grid.
update_a0_griddy <- function(sigma2, sp, n_grid = 200) {
  grid  <- seq(1e-6, 1 - 1e-6, length.out = n_grid)
  log_w <- vapply(grid, log_post_a0, numeric(1L), sigma2 = sigma2, sp = sp)
  log_w <- log_w - max(log_w)          # stabilise before exp
  w     <- exp(log_w); w <- w / sum(w)
  sample(grid, 1L, prob = w)
}
