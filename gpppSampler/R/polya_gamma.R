#' Draw from a Polya-Gamma Distribution
#'
#' Approximate sampler for PG(b, c) via the truncated sum-of-gammas series
#' (Polson, Scott & Windle, 2013, Sec. 2.2):
#'   J ~ (1 / 2*pi^2) * sum_{k=1}^Inf  g_k / ((k - 1/2)^2 + c^2/(4*pi^2)),
#'   g_k ~ Gamma(b, 1) iid.
#' This representation holds for any real b > 0, not just integers, which is
#' what makes it usable here: the historical Polya-Gamma draws in the GPPP
#' logit power-prior model (Section 2.4) need shape b = a0, a continuous
#' discounting weight in (0, 1), not the integer shape assumed by the
#' Devroye-style exact sampler most implementations use for b = 1.
#'
#' The series is truncated at `trunc` terms; the (deterministic) mean of the
#' dropped tail is added back as a bias correction, since simply truncating
#' would systematically under-draw. With `trunc = 200` the correction itself
#' is already below 1e-3 in relative size for the |c| ranges seen here, so
#' the remaining approximation error is far smaller than ordinary MCMC noise.
#'
#' @param n     Number of draws when b and c are scalars; overridden by
#'   max(n, length(b), length(c)) when either is a vector.
#' @param b     Shape parameter(s), b > 0 (may be a single continuous
#'   discounting weight such as a0, or an integer such as 1).
#' @param c     Tilting parameter(s) (the current linear predictor).
#' @param trunc Number of retained series terms.
#' @return Numeric vector of PG(b, c) draws, length max(n, length(b), length(c)).
#' @export
rpg <- function(n = 1L, b = 1, c = 0, trunc = 200L) {
  m <- max(n, length(b), length(c))
  b <- rep_len(as.double(b), m)
  c <- rep_len(as.double(c), m)
  stopifnot(all(b > 0))

  k     <- seq_len(trunc) - 0.5
  denom <- outer(k^2, (c^2) / (4 * pi^2), `+`)                # trunc x m
  g     <- matrix(rgamma(trunc * m, shape = rep(b, each = trunc), rate = 1),
                   nrow = trunc, ncol = m)

  tail_mean <- b / (2 * pi^2 * (trunc - 0.5))                 # dropped-tail correction
  as.vector(colSums(g / denom)) / (2 * pi^2) + tail_mean
}
