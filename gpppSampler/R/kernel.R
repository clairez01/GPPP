# Isotropic squared Euclidean distance matrix (n1 x n2), numerically safe.
sq_dist <- function(X1, X2) {
  r1 <- rowSums(X1^2)
  r2 <- rowSums(X2^2)
  D  <- outer(r1, r2, "+") - 2 * tcrossprod(X1, X2)
  pmax(D, 0)
}

# ARD-scaled squared distance: column j is divided by ell[j] before computing
# the isotropic squared distance.  ell may be a positive scalar (isotropic,
# recycled across all d columns) or a positive d-vector (ARD).
sq_dist_ard <- function(X1, X2, ell) {
  X1s <- sweep(X1, 2, ell, "/")
  X2s <- sweep(X2, 2, ell, "/")
  sq_dist(X1s, X2s)
}

#' Squared Exponential (RBF) Kernel Matrix with Optional ARD
#'
#' k(x, x') = eta2 * exp(-sum_j (x_j - x'_j)^2 / (2 * ell_j^2))
#'
#' When \code{ell} is a positive scalar the kernel is isotropic.  When
#' \code{ell} is a d-vector each input dimension gets its own length-scale
#' (Automatic Relevance Determination), which is recommended when covariates
#' have heterogeneous scales (e.g. age and weight on their natural scales).
#'
#' @param X1 n1 x d numeric matrix of inputs.
#' @param X2 n2 x d numeric matrix of inputs.
#' @param eta2 Signal variance (positive scalar).
#' @param ell  Length-scale: positive scalar (isotropic) or d-vector (ARD).
#' @return n1 x n2 covariance matrix.
#' @export
se_kernel <- function(X1, X2, eta2, ell) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  eta2 * exp(-sq_dist_ard(X1, X2, ell) / 2)
}

#' Matern 5/2 Kernel Matrix with Optional ARD
#'
#' k(x, x') = eta2 * (1 + sqrt(5)*r + 5*r^2/3) * exp(-sqrt(5)*r)
#' where r = sqrt(sum_j (x_j - x'_j)^2 / ell_j^2).
#'
#' Twice differentiable; often preferred over the squared exponential for
#' biological or clinical data where infinite smoothness is unrealistic.
#' Supports the same isotropic (scalar \code{ell}) or ARD (d-vector \code{ell})
#' parameterisation as \code{se_kernel}.
#'
#' @param X1 n1 x d numeric matrix of inputs.
#' @param X2 n2 x d numeric matrix of inputs.
#' @param eta2 Signal variance (positive scalar).
#' @param ell  Length-scale: positive scalar (isotropic) or d-vector (ARD).
#' @return n1 x n2 covariance matrix.
#' @export
matern52_kernel <- function(X1, X2, eta2, ell) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  r  <- sqrt(sq_dist_ard(X1, X2, ell))
  r5 <- sqrt(5) * r
  eta2 * (1 + r5 + r5^2 / 3) * exp(-r5)
}

#' Linear Kernel Matrix with Optional ARD
#'
#' k(x, x') = eta2 * sum_j x_j * x'_j / ell_j^2
#'
#' Equivalent to Bayesian linear regression with independent
#' N(0, eta2 / ell_j^2) priors on each regression coefficient.  ARD is
#' recommended when covariates have heterogeneous scales: without it the dot
#' product is dominated by large-scale features.  Larger ell_j shrinks the
#' j-th coefficient toward zero; smaller ell_j allows it to vary more freely.
#' The model intercept mu handles the overall mean; no separate bias term is
#' needed.
#'
#' @param X1 n1 x d numeric matrix of inputs.
#' @param X2 n2 x d numeric matrix of inputs.
#' @param eta2 Signal variance (positive scalar).
#' @param ell  Length-scale: positive scalar (isotropic) or d-vector (ARD).
#' @return n1 x n2 covariance matrix.
#' @export
linear_kernel <- function(X1, X2, eta2, ell) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  X1s <- sweep(X1, 2, ell, "/")
  X2s <- sweep(X2, 2, ell, "/")
  eta2 * tcrossprod(X1s, X2s)
}
