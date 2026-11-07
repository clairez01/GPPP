# Gauss-Legendre nodes/weights on [-1, 1] via the Golub-Welsch algorithm:
# the nodes are eigenvalues of the symmetric tridiagonal Jacobi matrix with
# zero diagonal and off-diagonal beta_k = k / sqrt(4k^2 - 1); weights are
# 2 * (first eigenvector component)^2.  Computed once and reused throughout
# the sampler run (n_nodes is fixed).
gauss_legendre_nodes <- function(n) {
  if (n < 2L) stop("gauss_legendre_nodes: n must be >= 2")
  J <- matrix(0, n, n)
  for (k in seq_len(n - 1L)) {
    b <- k / sqrt(4 * k^2 - 1)
    J[k, k + 1L] <- b
    J[k + 1L, k] <- b
  }
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  list(nodes = eig$values[ord], weights = 2 * eig$vectors[1, ord]^2)
}

# Fixed quadrature nodes/log-weights for integrating a function of tau > 0
# against the shape-rate Gamma(alpha0, beta0) prior on tau, via the rational
# map tau(x) = c*(1+x)/(1-x) for x in (-1, 1), with c = alpha0/beta0 (the
# prior mean of tau) chosen to centre the nodes near where the prior mass is.
#
# Returns list(tau = quadrature nodes, log_w = log(Gauss-Legendre weight *
# change-of-variables Jacobian)), so that for any integrand g(tau),
#   log(integral_0^Inf g(tau) dtau) ~= logsumexp(log_w + log(g(tau))).
make_tau_quadrature <- function(n_nodes, alpha0, beta0) {
  gl <- gauss_legendre_nodes(n_nodes)
  c0 <- alpha0 / beta0
  if (!is.finite(c0) || c0 <= 0) c0 <- 1
  x   <- gl$nodes
  tau <- c0 * (1 + x) / (1 - x)
  jac <- 2 * c0 / (1 - x)^2
  list(tau = tau, log_w = log(gl$weights) + log(jac))
}

# log density of the shape-rate Gamma(alpha0, beta0) prior on tau, i.e.
# pi0(tau) \propto tau^(alpha0 - 1) exp(-beta0 * tau).
log_pi0_tau <- function(tau, alpha0, beta0) {
  alpha0 * log(beta0) - lgamma(alpha0) + (alpha0 - 1) * log(tau) - beta0 * tau
}

# Numerically stable log(sum(exp(x))).
logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Eq. (51)-(52): the tau-integrated normalizer of the full GP power prior,
#   C(a0, theta) = integral_0^Inf C(a0, tau, theta) * pi0(tau) dtau,
# evaluated on the fixed quadrature nodes `quad` (see make_tau_quadrature),
# reusing the eigendecomposition Kn0n0 = U D U^T (D = eigenvalues of Kn0n0,
# w = U^T (y0 - mu), for the current kernel hyperparameters theta).
# Cost: O(n0) per quadrature node, O(n0 * length(quad$tau)) total.
log_C_a0_theta <- function(a0, D, w, alpha0, beta0, quad) {
  n0 <- length(D)
  vals <- vapply(seq_along(quad$tau), function(k) {
    tau <- quad$tau[k]
    nu  <- 1 / (a0 * tau)
    logC_tau <- -(n0 * a0 / 2) * log(2 * pi) - (n0 / 2) * log(a0) -
                (n0 * (1 - a0) / 2) * log(tau) -
                0.5 * sum(log(D + nu) + w^2 / (D + nu))
    quad$log_w[k] + logC_tau + log_pi0_tau(tau, alpha0, beta0)
  }, numeric(1))
  logsumexp(vals)
}
