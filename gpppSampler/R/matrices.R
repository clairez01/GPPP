# Build the four kernel sub-blocks for current (X, n) and historical (X0, n0).
# kernel: one of "se", "matern52", "linear".
build_kernel_mats <- function(X, X0, eta2, ell, kernel = "se") {
  kfn <- switch(kernel,
    se       = se_kernel,
    matern52 = matern52_kernel,
    linear   = linear_kernel,
    stop("Unknown kernel: '", kernel, "'. Choose 'se', 'matern52', or 'linear'.")
  )
  Knn   <- kfn(X,  X,  eta2, ell)
  Knn0  <- kfn(X,  X0, eta2, ell)
  Kn0n0 <- kfn(X0, X0, eta2, ell)
  list(Knn = Knn, Knn0 = Knn0, Kn0n = t(Knn0), Kn0n0 = Kn0n0)
}

# Cholesky with automatic diagonal jitter for numerical stability.
safe_chol <- function(M) {
  n <- nrow(M)
  for (jit in c(0, 1e-8, 1e-6, 1e-4, 1e-2)) {
    ch <- tryCatch(chol(M + jit * diag(n)), error = function(e) NULL)
    if (!is.null(ch)) return(ch)
  }
  NULL
}

# Joint covariance M = K + Sigma  (Eq. 19), size (n + n0) x (n + n0).
# Sigma = diag(sigma2 * I_n,  sigma2/a0 * I_n0)
build_M <- function(Km, sigma2, a0) {
  n  <- nrow(Km$Knn)
  n0 <- nrow(Km$Kn0n0)
  rbind(
    cbind(Km$Knn + sigma2 * diag(n),     Km$Knn0),
    cbind(Km$Kn0n, Km$Kn0n0 + (sigma2 / a0) * diag(n0))
  )
}

# Historical normaliser covariance M0 = Kn0n0 + sigma2/a0 * I_n0  (Eq. 19).
build_M0 <- function(Km, sigma2, a0) {
  n0 <- nrow(Km$Kn0n0)
  Km$Kn0n0 + (sigma2 / a0) * diag(n0)
}

# General joint covariance M = K + diag(var_vec), for heteroskedastic
# per-observation noise variance (used by the Polya-Gamma-augmented logit
# binary model, where each observation's working-likelihood variance
# 1/omega_i differs, rather than a single sigma2 with an a0-scaled
# historical block).
build_M_hetero <- function(Km, var_vec) {
  K_full <- rbind(cbind(Km$Knn, Km$Knn0), cbind(Km$Kn0n, Km$Kn0n0))
  K_full + diag(var_vec)
}
