#' Sample from a Truncated Normal Distribution
#'
#' Uses the inverse-CDF method to draw from N(mu, sigma^2) truncated to
#' (lower, upper).  All arguments except sigma are recycled so that
#' element-wise draws with varying means and bounds can be obtained in one call.
#'
#' @param n     Number of draws when mu is scalar; overridden by max(length(mu),
#'   length(lower), length(upper)) when any of those are vectors.
#' @param mu    Mean(s) of the underlying normal.
#' @param sigma Standard deviation (positive scalar).
#' @param lower Lower truncation bound(s); use \code{-Inf} for no lower bound.
#' @param upper Upper truncation bound(s); use \code{Inf} for no upper bound.
#' @return Numeric vector of length \code{n} (or the recycled length).
#' @export
rtnorm <- function(n = 1L, mu, sigma = 1, lower = -Inf, upper = Inf) {
  # Infer n from vector arguments when any have length > 1
  m <- max(n, length(mu), length(lower), length(upper))
  mu    <- rep_len(as.double(mu),    m)
  lower <- rep_len(as.double(lower), m)
  upper <- rep_len(as.double(upper), m)
  sigma <- as.double(sigma)[1L]

  p_lo <- ifelse(is.finite(lower), pnorm((lower - mu) / sigma), 0)
  p_hi <- ifelse(is.finite(upper), pnorm((upper - mu) / sigma), 1)

  # Degenerate intervals (nearly no probability mass inside the window):
  # return the finite boundary closest to mu.
  deg <- (p_hi - p_lo) < .Machine$double.eps^0.5
  out <- numeric(m)

  if (any(!deg)) {
    u <- runif(sum(!deg), p_lo[!deg], p_hi[!deg])
    out[!deg] <- mu[!deg] + sigma * qnorm(u)
  }
  if (any(deg)) {
    lo_d <- lower[deg]; hi_d <- upper[deg]; mu_d <- mu[deg]
    bd <- ifelse(!is.finite(lo_d), hi_d,
           ifelse(!is.finite(hi_d), lo_d,
           ifelse(abs(lo_d - mu_d) <= abs(hi_d - mu_d), lo_d, hi_d)))
    out[deg] <- bd
  }
  out
}

# Internal: draw the joint latent vector f* = (f, f0) from its GP posterior
# under a general per-observation noise-variance vector (used by the
# Polya-Gamma logit binary model, where variance_i = 1/omega_i varies by
# observation instead of following the scalar-sigma2 + a0-scaled-block
# pattern below). Returns NULL if the required Cholesky fails; falls back to
# the posterior mean on a failed covariance Cholesky.
draw_gp_posterior_hetero <- function(z_star, Km, var_vec) {
  M    <- build_M_hetero(Km, var_vec)
  ch_M <- safe_chol(M)
  if (is.null(ch_M)) return(NULL)

  n_total <- length(z_star)
  K_full  <- rbind(cbind(Km$Knn, Km$Knn0), cbind(Km$Kn0n, Km$Kn0n0))

  # Posterior mean: K M^{-1} z*
  Mi_z <- backsolve(ch_M, forwardsolve(t(ch_M), z_star))
  mu_f <- as.vector(K_full %*% Mi_z)

  # Posterior covariance: K - K M^{-1} K
  V     <- backsolve(ch_M, forwardsolve(t(ch_M), K_full))
  Sig_f <- K_full - K_full %*% V
  Sig_f <- 0.5 * (Sig_f + t(Sig_f))

  ch_Sf <- safe_chol(Sig_f)
  if (is.null(ch_Sf)) return(mu_f)
  mu_f + as.vector(t(ch_Sf) %*% rnorm(n_total))
}

# Internal: draw the joint latent vector f* = (f, f0) from its GP posterior
# given augmented observations z_star = c(z, z0), kernel matrices Km,
# noise variance sigma2, and discounting parameter a0 (Sigma = diag(sigma2 I_n,
# sigma2/a0 I_n0), Eq. 19). Used by the continuous/TTE branches, where the
# variance-inflation identity is exact because the underlying likelihood is
# genuinely Gaussian (see design note in sampler.R for why this does NOT
# carry over to the binary/logit case).
draw_gp_posterior <- function(z_star, Km, sigma2, a0) {
  n  <- nrow(Km$Knn)
  n0 <- nrow(Km$Kn0n0)
  draw_gp_posterior_hetero(z_star, Km, c(rep(sigma2, n), rep(sigma2 / a0, n0)))
}

# Internal: draw the Polya-Gamma augmentation variables for the logit-link
# power prior binary model (Section 2.4, revised). Current data use weight
# b = 1 (the full, undiscounted logistic likelihood); historical data use
# b = a0, which is the exact Polya-Gamma identity for the **observed-data**
# logistic likelihood raised to the power a0 -- unlike the previous
# probit/Albert-Chib scheme, which discounted the complete-data
# augmented-normal likelihood, a different, non-equivalent quantity (see
# design note in sampler.R).
#
# eta, eta0 - current linear predictors mu + f, mu + f0
# y, y0     - binary outcomes in {0, 1}
# a0        - current discounting parameter
# Returns list(omega, omega0, z, z0), where z_i = kappa_i / omega_i is the
# Gaussian pseudo-response (kappa = y - 0.5 for current, a0*(y0 - 0.5) for
# historical) with conditional variance 1/omega_i.
augment_binary_pg <- function(eta, eta0, y, y0, a0) {
  omega  <- rpg(length(eta),  b = 1,  c = eta)
  omega0 <- rpg(length(eta0), b = a0, c = eta0)
  list(
    omega  = omega,
    omega0 = omega0,
    z      = (y - 0.5)        / omega,
    z0     = (a0 * (y0 - 0.5)) / omega0
  )
}

# Internal: one data-augmentation sweep for time-to-event (AFT) outcomes.
#
# For censored observations the true latent log-event-time T_i is unknown
# beyond the censoring bound C_i.  We draw T_i from a truncated normal
# bounded below by the observed censoring time.  Observed events are unchanged.
#
# T_aug, T0_aug   - current augmented (log) event times
# y_obs, y0_obs   - observed min(T, C): event times for events, censoring times for censored
# delta, delta0   - event indicators (1 = event observed, 0 = censored)
# Km              - kernel matrix list
# sigma2, a0      - current noise variance and discounting parameter
# mu              - current intercept
#
# The observation model is T_i = mu + f(x_i) + eps_i, so the GP draw must
# condition on the mu-centred pseudo-observations (T_aug - mu, T0_aug - mu);
# the truncated-normal draw for censored observations is then centred back
# at mu + f.
augment_tte <- function(T_aug, T0_aug, y_obs, y0_obs, delta, delta0,
                        Km, sigma2, a0, mu) {
  n  <- length(T_aug)
  n0 <- length(T0_aug)

  f_star <- draw_gp_posterior(c(T_aug - mu, T0_aug - mu), Km, sigma2, a0)
  if (is.null(f_star)) return(list(T_aug = T_aug, T0_aug = T0_aug))

  f  <- f_star[seq_len(n)]
  f0 <- f_star[seq_len(n0) + n]

  # T_i | f_i, mu, C_i, delta_i=0 ~ TN(mu + f_i, sigma2, C_i, Inf)
  cens_cur <- which(delta  == 0L)
  cens_his <- which(delta0 == 0L)

  if (length(cens_cur) > 0L)
    T_aug[cens_cur]  <- rtnorm(length(cens_cur),
                               mu    = mu + f[cens_cur],
                               sigma = sqrt(sigma2),
                               lower = y_obs[cens_cur])
  if (length(cens_his) > 0L)
    T0_aug[cens_his] <- rtnorm(length(cens_his),
                               mu    = mu + f0[cens_his],
                               sigma = sqrt(sigma2),
                               lower = y0_obs[cens_his])

  list(T_aug = T_aug, T0_aug = T0_aug)
}
