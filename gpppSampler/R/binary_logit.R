# Logit-link, Polya-Gamma-augmented power prior for binary outcomes
# (Section 2.4, revised). Discounting is applied to the true observed-data
# logistic likelihood L0(f0)^a0, not to the complete-data (augmented-latent)
# likelihood -- see the design note at the top of the binary branch in
# sampler.R for the full derivation and why the two are not equivalent.

# Sum of the observed-data Bernoulli(logit) log-likelihood at linear
# predictor eta = mu + f.
log_lik_logit_sum <- function(eta, y) {
  sum(ifelse(y == 1L, plogis(eta, log.p = TRUE), plogis(-eta, log.p = TRUE)))
}

# Laplace approximation to the historical normalizing constant
#   C(a0, theta) = integral  L0(mu + f0)^a0 * N(f0 | 0, K) d f0,
# which has no closed form for the logistic likelihood and is required for
# a0's full conditional (the a0-weighted Polya-Gamma augmentation gives exact
# draws of f0, but not this normalizer -- it does not cancel).
#
# Standard GP-classification Newton-Raphson Laplace approximation (Rasmussen
# & Williams, 2006, Algorithm 3.1), generalized so the discount a0 scales
# both the log-likelihood and its Hessian weight W. The objective
#   Psi(f0) = a0 * loglik(mu + f0, y0) - 0.5 f0' K^{-1} f0
# is concave (log-logistic-likelihood is concave in f0, GP prior is a
# concave quadratic), so Newton's method converges reliably from any start.
#
# K       - historical prior covariance Kn0n0(theta), n0 x n0
# y0      - historical binary outcomes
# mu      - current intercept (linear predictor is mu + f0)
# a0      - discounting weight, the Polya-Gamma / Laplace shape parameter
# f0_init - warm start for the mode (e.g. the previous grid point's mode);
#           NULL starts at 0
# Returns list(f0_hat, logC, converged), or NULL if a Cholesky fails.
laplace_logit_gp <- function(K, y0, mu, a0, f0_init = NULL,
                              tol = 1e-8, max_iter = 100L) {
  n0 <- length(y0)
  f0 <- if (is.null(f0_init)) numeric(n0) else f0_init
  I_n0 <- diag(n0)
  psi_prev  <- -Inf
  ch_B <- NULL; a_vec <- NULL; converged <- FALSE

  for (iter in seq_len(max_iter)) {
    eta   <- mu + f0
    p     <- plogis(eta)
    W     <- a0 * p * (1 - p)
    sqrtW <- sqrt(W)
    B     <- I_n0 + outer(sqrtW, sqrtW) * K
    ch_B  <- safe_chol(B)
    if (is.null(ch_B)) return(NULL)

    b     <- W * f0 + a0 * (y0 - p)
    v     <- backsolve(ch_B, forwardsolve(t(ch_B), sqrtW * as.vector(K %*% b)))
    a_vec <- b - sqrtW * v
    f0    <- as.vector(K %*% a_vec)

    psi <- -0.5 * sum(a_vec * f0) + a0 * log_lik_logit_sum(mu + f0, y0)
    if (abs(psi - psi_prev) < tol * (abs(psi) + 1)) { converged <- TRUE; break }
    psi_prev <- psi
  }

  logC <- -0.5 * sum(a_vec * f0) + a0 * log_lik_logit_sum(mu + f0, y0) -
    sum(log(diag(ch_B)))
  list(f0_hat = f0, logC = logC, converged = converged)
}

# a0's full conditional (up to the Uniform(0,1) base prior handled by the
# caller): log[ L0(f0_current)^a0 / C_Laplace(a0, theta) ]. Only the
# normalizing constant is approximated -- the numerator uses the actual,
# exactly-drawn current f0, not the Laplace mode.
#
# Returns list(lp, f0_hat) so the Laplace mode can warm-start the next call
# (griddy-Gibbs / slice sampling evaluate this many times per iteration).
log_post_a0_laplace <- function(a0, f0_current, y0, K0, mu, f0_init = NULL) {
  if (a0 <= 0 || a0 >= 1) return(list(lp = -Inf, f0_hat = f0_init))
  lap <- laplace_logit_gp(K0, y0, mu, a0, f0_init = f0_init)
  if (is.null(lap)) return(list(lp = -Inf, f0_hat = f0_init))
  num <- a0 * log_lik_logit_sum(mu + f0_current, y0)
  list(lp = num - lap$logC, f0_hat = lap$f0_hat)
}

# Griddy-Gibbs draw for a0 (logit branch). Each grid point requires one
# Newton solve for C_Laplace; warm-starting each solve at the previous grid
# point's mode keeps this from being 200 independent from-scratch solves.
update_a0_logit_griddy <- function(f0_current, y0, K0, mu, n_grid = 200L,
                                   log_prior = 0) {
  grid   <- seq(1e-6, 1 - 1e-6, length.out = n_grid)
  log_w  <- numeric(n_grid)
  f0_hat <- NULL
  for (i in seq_along(grid)) {
    res <- log_post_a0_laplace(grid[i], f0_current, y0, K0, mu, f0_init = f0_hat)
    log_w[i] <- res$lp
    if (!is.null(res$f0_hat)) f0_hat <- res$f0_hat
  }
  log_w <- log_w + log_prior
  log_w <- log_w - max(log_w)
  w <- exp(log_w); w <- w / sum(w)
  sample(grid, 1L, prob = w)
}

# Collapsed Gaussian marginal likelihood for theta = (eta2, ell), conditional
# on the current Polya-Gamma weights omega (both current and historical
# blocks) -- exact and closed-form given omega, since z_star | f, theta is
# then a standard Gaussian working-response model. No normalizing-constant
# approximation is needed here: C(a0, theta) only enters a0's own full
# conditional (log_post_a0_laplace), not theta's.
log_marginal_lik_hetero <- function(z_star, Km, var_vec) {
  K_full <- rbind(cbind(Km$Knn, Km$Knn0), cbind(Km$Kn0n, Km$Kn0n0))
  M <- K_full + diag(var_vec)
  ch_M <- safe_chol(M)
  if (is.null(ch_M)) return(-Inf)
  log_mvn(z_star, ch_M)
}

# Conjugate Gibbs update for the intercept mu given the Polya-Gamma pseudo
# responses z, z0, the current latent draws f, f0, and weights omega, omega0
# (z_i = mu + f_i + N(0, 1/omega_i)).
update_mu_gibbs_pg <- function(z, z0, f, f0, omega, omega0, m0, v0) {
  prior_prec <- if (is.finite(v0)) 1 / v0 else 0
  prior_wt   <- if (is.finite(v0)) m0 / v0 else 0
  prec     <- sum(omega) + sum(omega0) + prior_prec
  mean_num <- sum(omega * (z - f)) + sum(omega0 * (z0 - f0)) + prior_wt
  if (prec <= 0) return(m0)
  post_var  <- 1 / prec
  post_mean <- post_var * mean_num
  rnorm(1L, post_mean, sqrt(post_var))
}
