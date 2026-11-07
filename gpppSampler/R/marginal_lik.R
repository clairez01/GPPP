# log N(y | 0, Sigma) via its upper-triangular Cholesky factor.
log_mvn <- function(y, chol_S) {
  n <- length(y)
  alpha <- backsolve(chol_S, forwardsolve(t(chol_S), y))
  -0.5 * n * log(2 * pi) - sum(log(diag(chol_S))) - 0.5 * drop(crossprod(y, alpha))
}

# Eq. (20) extended to include a scalar intercept mu:
#   log p(y | y0, a0, sigma2, theta, mu)
#   = log N(y* - mu | 0, M) - log N(y0 - mu | 0, M0)
# where y* = c(y, y0).
log_marginal_lik <- function(y, y0, Km, sigma2, a0, mu = 0) {
  M  <- build_M(Km, sigma2, a0)
  M0 <- build_M0(Km, sigma2, a0)
  ch_M  <- safe_chol(M)
  ch_M0 <- safe_chol(M0)
  if (is.null(ch_M) || is.null(ch_M0)) return(-Inf)
  yc  <- y  - mu
  y0c <- y0 - mu
  log_mvn(c(yc, y0c), ch_M) - log_mvn(y0c, ch_M0)
}

# Conjugate Gibbs update for the scalar intercept mu.
# Model: y = mu + f + eps,  f ~ GP(0, K),  mu ~ N(m0, v0).
# Full conditional:  mu | rest ~ N(post_mean, post_var).
# Returns list(mu, lml) with the new draw and updated log marginal likelihood.
update_mu_gibbs <- function(y, y0, Km, sigma2, a0, m0, v0) {
  n  <- length(y)
  n0 <- length(y0)
  M  <- build_M(Km, sigma2, a0)
  M0 <- build_M0(Km, sigma2, a0)
  ch_M  <- safe_chol(M)
  ch_M0 <- safe_chol(M0)
  if (is.null(ch_M) || is.null(ch_M0)) return(NULL)

  y_star  <- c(y, y0)
  ones_np <- rep(1, n + n0)
  ones_n0 <- rep(1, n0)

  Mi1   <- backsolve(ch_M,  forwardsolve(t(ch_M),  ones_np))
  M0i1  <- backsolve(ch_M0, forwardsolve(t(ch_M0), ones_n0))
  Miy   <- backsolve(ch_M,  forwardsolve(t(ch_M),  y_star))
  M0iy0 <- backsolve(ch_M0, forwardsolve(t(ch_M0), y0))

  prior_prec <- if (is.finite(v0)) 1 / v0 else 0
  prior_wt   <- if (is.finite(v0)) m0 / v0 else 0

  prec     <- sum(Mi1) - sum(M0i1) + prior_prec
  mean_num <- sum(Miy) - sum(M0iy0) + prior_wt

  if (prec <= 0) return(NULL)

  post_var  <- 1 / prec
  post_mean <- post_var * mean_num

  mu_new  <- rnorm(1L, post_mean, sqrt(post_var))
  yc      <- y  - mu_new
  y0c     <- y0 - mu_new
  lml_new <- log_mvn(c(yc, y0c), ch_M) - log_mvn(y0c, ch_M0)
  list(mu = mu_new, lml = lml_new)
}

# Full GP power prior (Section 3.5, unknown tau): collapsed marginal
# likelihood for the kernel hyperparameters theta, at fixed (a0, tau).
#
# Eq. (49): log p(theta | y, y0, a0, tau) [up to log pi(theta), added by the
# caller] equals log N(y* | 0, M) - log C(a0, theta), where M = K(theta) + Sigma
# with Sigma = diag(I_n/tau, I_n0/(a0*tau)) (Eq. 29) -- identical to
# build_M(Km, sigma2 = 1/tau, a0) -- and C(a0, theta) is the tau-integrated
# normalizer (Eq. 51-52), recomputed at each proposed theta because Kn0n0(theta)
# changes with it.
log_marginal_lik_full <- function(y, y0, Km, tau, a0, mu, alpha0, beta0, quad) {
  sigma2 <- 1 / tau
  M  <- build_M(Km, sigma2, a0)
  ch_M <- safe_chol(M)
  if (is.null(ch_M)) return(-Inf)

  yc  <- y  - mu
  y0c <- y0 - mu
  lp_data <- log_mvn(c(yc, y0c), ch_M)

  eig_K00 <- eigen(Km$Kn0n0, symmetric = TRUE)
  D <- pmax(eig_K00$values, 0)
  w <- as.vector(t(eig_K00$vectors) %*% y0c)
  logC <- log_C_a0_theta(a0, D, w, alpha0, beta0, quad)

  lp_data - logC
}

# Full GP power prior (Section 3.4, unknown tau): conjugate Gibbs draw for the
# noise precision tau, given the retained latent vectors f, f0 (Eq. 48).
# a0 scales both the historical degrees of freedom and the historical sum of
# squares -- the gamma analogue of discounting.
update_tau_gibbs <- function(y, y0, f, f0, a0, mu, alpha0, beta0) {
  n  <- length(y)
  n0 <- length(y0)
  SS <- sum((y - mu - f)^2) + a0 * sum((y0 - mu - f0)^2)
  rgamma(1L, shape = alpha0 + (n + a0 * n0) / 2, rate = beta0 + SS / 2)
}

# Full GP power prior: conjugate Gibbs draw for the intercept mu, given the
# retained latent vectors f, f0 and the current tau, a0.  Because f, f0 are
# already drawn, this is a plain Bayesian-regression update on the residuals
# (y - f) and (y0 - f0), with the historical residuals discounted by a0 -- no
# M / M0 machinery required (contrast update_mu_gibbs(), used when f* is
# marginalized out).
update_mu_gibbs_full <- function(y, y0, f, f0, tau, a0, m0, v0) {
  n  <- length(y)
  n0 <- length(y0)
  prior_prec <- if (is.finite(v0)) 1 / v0 else 0
  prior_wt   <- if (is.finite(v0)) m0 / v0 else 0

  prec     <- tau * n + a0 * tau * n0 + prior_prec
  mean_num <- tau * sum(y - f) + a0 * tau * sum(y0 - f0) + prior_wt

  if (prec <= 0) return(m0)
  post_var  <- 1 / prec
  post_mean <- post_var * mean_num
  rnorm(1L, post_mean, sqrt(post_var))
}
