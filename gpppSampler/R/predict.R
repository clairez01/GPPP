#' Posterior Predictive Draws of the Mean Response
#'
#' Compiled hot path: for a supported kernel (se, matern52, linear), the
#' per-draw loop (kernel matrix build, Cholesky, conditional mean/covariance,
#' joint draw) runs in one compiled call (predict_gppp_draws_cpp) instead of
#' an R for-loop over n_save draws. Falls back to the R loop for kernels not
#' covered by the compiled path (there are none currently) or when the
#' compiled call fails for some other reason.
#'
#' Binary outcomes (logit link, Section 2.4) bypass the compiled path
#' entirely: fit$samples$f/f0 already hold an exact joint GP draw at each
#' retained iteration (no noise term, since f* are exact function values,
#' not noisy pseudo-observations), so prediction is the plain GP conditional
#' on f* -- apply plogis() to the returned draws to get probabilities.
#' @export
predict_gppp <- function(fit, y, X, y0, X0, Xp = NULL) {
  y  <- as.vector(y);  X  <- as.matrix(X)
  y0 <- as.vector(y0); X0 <- as.matrix(X0)
  if (is.null(Xp)) Xp <- X
  Xp <- as.matrix(Xp)

  kernel <- fit$kernel
  samp   <- fit$samples
  n_save <- length(samp$eta2)

  if (identical(fit$outcome_type, "binary")) {
    kfn <- switch(kernel, se = se_kernel, matern52 = matern52_kernel, linear = linear_kernel)
    np    <- nrow(Xp)
    F_out <- matrix(NA_real_, nrow = n_save, ncol = np)
    for (i in seq_len(n_save)) {
      eta2 <- samp$eta2[i]; ell <- as.vector(samp$ell[i, ])
      f_i  <- c(samp$f[i, ], samp$f0[i, ])

      K_full <- rbind(cbind(kfn(X,  X,  eta2, ell), kfn(X,  X0, eta2, ell)),
                       cbind(kfn(X0, X,  eta2, ell), kfn(X0, X0, eta2, ell)))
      ch_K <- safe_chol(K_full)
      if (is.null(ch_K)) next

      Kp_star <- cbind(kfn(Xp, X, eta2, ell), kfn(Xp, X0, eta2, ell))
      alpha   <- backsolve(ch_K, forwardsolve(t(ch_K), f_i))
      mu_p    <- as.vector(Kp_star %*% alpha) + samp$mu[i]

      V       <- backsolve(ch_K, forwardsolve(t(ch_K), t(Kp_star)))
      Sigma_p <- kfn(Xp, Xp, eta2, ell) - Kp_star %*% V
      Sigma_p <- 0.5 * (Sigma_p + t(Sigma_p))

      ch_Sp <- safe_chol(Sigma_p)
      if (is.null(ch_Sp)) { F_out[i, ] <- mu_p; next }
      F_out[i, ] <- mu_p + as.vector(t(ch_Sp) %*% rnorm(np))
    }
    return(F_out)
  }

  has_aug <- !is.null(samp$y_aug)
  y_eff  <- if (has_aug) samp$y_aug  else matrix(y,  nrow = 1)
  y0_eff <- if (has_aug) samp$y0_aug else matrix(y0, nrow = 1)

  if (kernel %in% c("se", "matern52", "linear")) {
    return(predict_gppp_draws_cpp(
      X, X0, Xp, kernel,
      eta2_draws = samp$eta2,
      ell_draws  = samp$ell,
      tau_draws  = samp$tau,
      a0_draws   = samp$a0,
      mu_draws   = samp$mu,
      y_eff      = y_eff,
      y0_eff     = y0_eff,
      rng_seeds_unused = 0
    ))
  }

  # Fallback: R reference loop (unreachable given match.arg() in
  # gppp_sampler(), kept for completeness).
  kfn <- switch(kernel,
    se       = se_kernel,
    matern52 = matern52_kernel,
    linear   = linear_kernel
  )
  np <- nrow(Xp)
  F_out <- matrix(NA_real_, nrow = n_save, ncol = np)
  for (i in seq_len(n_save)) {
    eta2   <- samp$eta2[i]
    ell    <- as.vector(samp$ell[i, ])
    sigma2 <- 1 / samp$tau[i]
    a0     <- samp$a0[i]
    mu_i   <- samp$mu[i]
    y_i    <- if (has_aug) samp$y_aug[ i, ] else y
    y0_i   <- if (has_aug) samp$y0_aug[i, ] else y0
    y_star_c <- c(y_i - mu_i, y0_i - mu_i)

    Km   <- build_kernel_mats(X, X0, eta2, ell, kernel)
    ch_M <- safe_chol(build_M(Km, sigma2, a0))
    if (is.null(ch_M)) next

    Kp_star <- cbind(kfn(Xp, X, eta2, ell), kfn(Xp, X0, eta2, ell))
    alpha <- backsolve(ch_M, forwardsolve(t(ch_M), y_star_c))
    mu_p  <- as.vector(Kp_star %*% alpha) + mu_i

    V       <- backsolve(ch_M, forwardsolve(t(ch_M), t(Kp_star)))
    Sigma_p <- kfn(Xp, Xp, eta2, ell) - Kp_star %*% V
    Sigma_p <- 0.5 * (Sigma_p + t(Sigma_p))

    ch_Sp <- safe_chol(Sigma_p)
    if (is.null(ch_Sp)) { F_out[i, ] <- mu_p; next }
    F_out[i, ] <- mu_p + as.vector(t(ch_Sp) %*% rnorm(np))
  }
  F_out
}
