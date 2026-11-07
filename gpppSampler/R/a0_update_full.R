# Full GP power prior (Section 3.3, unknown tau): collapsed full conditional
# for a0 given (tau, theta), with f* marginalized out.
#
# Reuses precompute_spectral(y, y0, Km, sigma2 = 1/tau, mu) from
# a0_update.R for the numerator spectral quantities (lambda = eigenvalues of
# Ktilde = Kn0n0 - Kn0n (Knn + I/tau)^{-1} Knn0, v = Q^T(y0 - mu_{y0|y})) and
# the denominator eigendecomposition of Kn0n0 (d = eigenvalues, w = U^T y0c)
# that feeds the tau-integrated normalizer log_C_a0_theta() (Eq. 51-52).
#
# Eq. (43): log p(a0 | y, y0, tau, theta) [up to the additive constant
# log pi(a0), which is 0 for the Uniform(0,1) prior] equals
#   -1/2 sum_i [log(lambda_i + nu) + v_i^2/(lambda_i + nu)]     nu = 1/(a0*tau)
#   - log C(a0, theta)
#   - (n0/2) log(a0) - (n0(1-a0)/2) log(tau / 2*pi)
log_post_a0_full <- function(a0, tau, sp, alpha0, beta0, quad) {
  if (a0 <= 0 || a0 >= 1) return(-Inf)
  n0  <- length(sp$lambda)
  nu  <- 1 / (a0 * tau)
  num  <- -0.5 * sum(log(sp$lambda + nu) + sp$v^2 / (sp$lambda + nu))
  logC <- log_C_a0_theta(a0, sp$d, sp$w, alpha0, beta0, quad)
  num - logC - (n0 / 2) * log(a0) - (n0 * (1 - a0) / 2) * log(tau / (2 * pi))
}

# log pi(a0) for a0 ~ Beta(s1, s2), up to a constant. Used by the slice paths.
log_prior_a0 <- function(a0, s1, s2) {
  (s1 - 1) * log(a0) + (s2 - 1) * log(1 - a0)
}

# Log prior weight of each griddy-Gibbs grid point under a0 ~ Beta(s1, s2):
# the Beta mass of the k-th of n_grid equal-width cells of [0, 1]. Using cell
# masses rather than the density at the grid point keeps U-shaped priors
# (s1 or s2 < 1) from piling weight onto the end points (1e-6, 1 - 1e-6),
# where the density is unbounded. Uniform(0, 1) gives equal weights, i.e. the
# original unweighted grid.
log_prior_a0_grid <- function(n_grid, s1, s2) {
  log(diff(pbeta(seq(0, 1, length.out = n_grid + 1L), s1, s2)))
}

# Griddy-Gibbs draw for a0 under the full power prior.
update_a0_griddy_full <- function(tau, sp, alpha0, beta0, quad, n_grid = 200) {
  grid  <- seq(1e-6, 1 - 1e-6, length.out = n_grid)
  log_w <- vapply(grid, log_post_a0_full, numeric(1L),
                   tau = tau, sp = sp, alpha0 = alpha0, beta0 = beta0, quad = quad)
  log_w <- log_w - max(log_w)
  w     <- exp(log_w); w <- w / sum(w)
  sample(grid, 1L, prob = w)
}
