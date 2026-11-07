// gpppSampler: compiled hot-path replacements for the sampler's R functions.
//
// Ported (per Rprof evidence on a representative fit):
//   1. log_post_a0_grid_cpp    -- replaces update_a0_griddy_full()'s
//      vapply(grid, log_post_a0_full, ...) [200-point grid x 40-node tau
//      quadrature per grid point = 8000 R-level calls/iteration; 54% of
//      sampled runtime].
//   2. se_kernel_cpp / matern52_kernel_cpp / linear_kernel_cpp
//      + log_marginal_lik_full_cpp -- replaces build_kernel_mats() +
//      build_M() + safe_chol() + log_mvn() + eigen(Kn0n0) +
//      log_C_a0_theta(), combined into one call per likelihood evaluation
//      (used repeatedly inside the eta2/ell slice-sampling steps; ~35% of
//      sampled runtime, much of which is R-level call-chaining/allocation
//      overhead around small (n+n0)-sized matrices).
//   3. precompute_spectral_cpp -- same Armadillo building blocks, ported for
//      completeness/consistency (called once per iteration; not itself a
//      profiled hot spot, but shares all its machinery with #2).
//   4. predict_gppp_draws_cpp  -- the predict_gppp() per-draw loop (build
//      kernel mats + chol(M) + conditional mean/covariance + Cholesky
//      draw), run once per retained posterior sample; not covered by the
//      gppp_sampler() profile above but structurally identical cost per
//      call and invoked heavily (n_save x n_prediction_targets) by
//      gppp_sims.R's Bayesian-bootstrap ATE step.
//
// All R-level control flow (the MCMC loop, slice-sampling stepping-out /
// shrinkage, Gibbs updates for tau/mu) stays in R -- only the repeated small
// numerical evaluations move to C++.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

static const double LOG2PI = std::log(2.0 * M_PI);

// ---------------------------------------------------------------------
// Kernel matrices (ARD squared exponential / Matern 5/2 / linear)
// ---------------------------------------------------------------------

// Squared Euclidean distance matrix between rows of X1s and X2s (already
// ARD-scaled, i.e. divided by ell), matching R's sq_dist().
static arma::mat sq_dist_scaled(const arma::mat& X1s, const arma::mat& X2s) {
  arma::vec r1 = arma::sum(arma::square(X1s), 1);
  arma::vec r2 = arma::sum(arma::square(X2s), 1);
  arma::mat D = arma::repmat(r1, 1, X2s.n_rows) + arma::repmat(r2.t(), X1s.n_rows, 1)
                - 2.0 * (X1s * X2s.t());
  D.elem(arma::find(D < 0)).zeros();
  return D;
}

static arma::mat ard_scale(const arma::mat& X, const arma::vec& ell) {
  return X.each_row() / ell.t();
}

// [[Rcpp::export]]
arma::mat se_kernel_cpp(const arma::mat& X1, const arma::mat& X2,
                         double eta2, const arma::vec& ell) {
  arma::mat D = sq_dist_scaled(ard_scale(X1, ell), ard_scale(X2, ell));
  return eta2 * arma::exp(-D / 2.0);
}

// [[Rcpp::export]]
arma::mat matern52_kernel_cpp(const arma::mat& X1, const arma::mat& X2,
                               double eta2, const arma::vec& ell) {
  arma::mat D = sq_dist_scaled(ard_scale(X1, ell), ard_scale(X2, ell));
  arma::mat r = arma::sqrt(D);
  arma::mat r5 = std::sqrt(5.0) * r;
  return eta2 * ((1.0 + r5 + arma::square(r5) / 3.0) % arma::exp(-r5));
}

// [[Rcpp::export]]
arma::mat linear_kernel_cpp(const arma::mat& X1, const arma::mat& X2,
                             double eta2, const arma::vec& ell) {
  arma::mat X1s = ard_scale(X1, ell);
  arma::mat X2s = ard_scale(X2, ell);
  return eta2 * (X1s * X2s.t());
}

static arma::mat kernel_dispatch(const std::string& kernel,
                                  const arma::mat& X1, const arma::mat& X2,
                                  double eta2, const arma::vec& ell) {
  if (kernel == "se")       return se_kernel_cpp(X1, X2, eta2, ell);
  if (kernel == "matern52") return matern52_kernel_cpp(X1, X2, eta2, ell);
  if (kernel == "linear")   return linear_kernel_cpp(X1, X2, eta2, ell);
  Rcpp::stop("Unknown kernel: " + kernel);
}

// Cholesky with automatic diagonal jitter for numerical stability, matching
// R's safe_chol(). Returns success flag; fills `ch` (upper-triangular, R
// chol() convention) on success.
static bool safe_chol_cpp(const arma::mat& M, arma::mat& ch) {
  static const double jitters[5] = {0.0, 1e-8, 1e-6, 1e-4, 1e-2};
  const arma::uword n = M.n_rows;
  for (double jit : jitters) {
    arma::mat Mj = (jit > 0.0) ? (M + jit * arma::eye(n, n)) : M;
    if (arma::chol(ch, Mj)) return true;  // arma::chol default = upper triangular
  }
  return false;
}

// log N(y | 0, Sigma) given the upper-triangular Cholesky factor of Sigma.
static double log_mvn_cpp(const arma::vec& y, const arma::mat& ch_S) {
  arma::vec alpha = arma::solve(arma::trimatu(ch_S),
                     arma::solve(arma::trimatl(ch_S.t()), y));
  double log_det = arma::sum(arma::log(ch_S.diag()));
  return -0.5 * y.n_elem * LOG2PI - log_det - 0.5 * arma::dot(y, alpha);
}

// ---------------------------------------------------------------------
// tau-integrated normalizer C(a0, theta) -- Gauss-Legendre quadrature over
// the fixed nodes precomputed once in R (make_tau_quadrature()).
// ---------------------------------------------------------------------

static double log_pi0_tau_cpp(double tau, double alpha0, double beta0) {
  return alpha0 * std::log(beta0) - std::lgamma(alpha0) +
         (alpha0 - 1.0) * std::log(tau) - beta0 * tau;
}

// D, w: eigenvalues/rotated-y0 of Kn0n0 (see precompute_spectral_cpp / R's
// eigen(Kn0n0)). Returns log C(a0, theta).
static double log_C_a0_theta_cpp(double a0, const arma::vec& D, const arma::vec& w,
                                  double alpha0, double beta0,
                                  const arma::vec& quad_tau, const arma::vec& quad_logw) {
  const int n0 = D.n_elem;
  const int nq = quad_tau.n_elem;
  arma::vec vals(nq);
  for (int k = 0; k < nq; ++k) {
    double tau_k = quad_tau[k];
    double nu    = 1.0 / (a0 * tau_k);
    double logC_tau = -(n0 * a0 / 2.0) * LOG2PI - (n0 / 2.0) * std::log(a0) -
                       (n0 * (1.0 - a0) / 2.0) * std::log(tau_k) -
                       0.5 * arma::sum(arma::log(D + nu) + arma::square(w) / (D + nu));
    vals[k] = quad_logw[k] + logC_tau + log_pi0_tau_cpp(tau_k, alpha0, beta0);
  }
  double m = vals.max();
  if (!std::isfinite(m)) return m;
  return m + std::log(arma::sum(arma::exp(vals - m)));
}

// ---------------------------------------------------------------------
// #1 target: the a0 griddy-Gibbs grid, evaluated entirely in C++.
// Mirrors a0_update_full.R's log_post_a0_full(), vectorized over `grid`.
// sp$lambda, sp$v, sp$d, sp$w are precompute_spectral()'s outputs (still
// computed once per iteration -- via precompute_spectral_cpp() below, or
// the existing R version -- and passed in here).
// ---------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector log_post_a0_grid_cpp(NumericVector grid, double tau,
                                    arma::vec lambda, arma::vec v,
                                    arma::vec d, arma::vec w,
                                    double alpha0, double beta0,
                                    arma::vec quad_tau, arma::vec quad_logw) {
  const int n0 = d.n_elem;
  const int ng = grid.size();
  NumericVector out(ng);
  for (int i = 0; i < ng; ++i) {
    double a0 = grid[i];
    if (a0 <= 0.0 || a0 >= 1.0) { out[i] = R_NegInf; continue; }
    double nu  = 1.0 / (a0 * tau);
    double num = -0.5 * arma::sum(arma::log(lambda + nu) + arma::square(v) / (lambda + nu));
    double logC = log_C_a0_theta_cpp(a0, d, w, alpha0, beta0, quad_tau, quad_logw);
    out[i] = num - logC - (n0 / 2.0) * std::log(a0) -
             (n0 * (1.0 - a0) / 2.0) * std::log(tau / (2.0 * M_PI));
  }
  return out;
}

// Same log-posterior at a single a0 (used by the a0 slice-sampling path,
// i.e. slice = c(..., "a0")), so that path also avoids the R-level
// build/eigen chain per proposal.
// [[Rcpp::export]]
double log_post_a0_single_cpp(double a0, double tau,
                               arma::vec lambda, arma::vec v,
                               arma::vec d, arma::vec w,
                               double alpha0, double beta0,
                               arma::vec quad_tau, arma::vec quad_logw) {
  if (a0 <= 0.0 || a0 >= 1.0) return R_NegInf;
  const int n0 = d.n_elem;
  double nu  = 1.0 / (a0 * tau);
  double num = -0.5 * arma::sum(arma::log(lambda + nu) + arma::square(v) / (lambda + nu));
  double logC = log_C_a0_theta_cpp(a0, d, w, alpha0, beta0, quad_tau, quad_logw);
  return num - logC - (n0 / 2.0) * std::log(a0) -
         (n0 * (1.0 - a0) / 2.0) * std::log(tau / (2.0 * M_PI));
}

// ---------------------------------------------------------------------
// #2 target: combined kernel-build + M-assembly + Cholesky + log_mvn +
// eigen(Kn0n0) + quadrature, in one call. Mirrors log_marginal_lik_full()
// in marginal_lik.R exactly (continuous/tte "full prior" branch).
// ---------------------------------------------------------------------

// [[Rcpp::export]]
List log_marginal_lik_full_cpp(const arma::mat& X, const arma::mat& X0,
                                double eta2, const arma::vec& ell,
                                std::string kernel,
                                const arma::vec& y, const arma::vec& y0,
                                double tau, double a0, double mu,
                                double alpha0, double beta0,
                                const arma::vec& quad_tau, const arma::vec& quad_logw) {
  const arma::uword n  = X.n_rows;
  const arma::uword n0 = X0.n_rows;
  const double sigma2 = 1.0 / tau;

  arma::mat Knn   = kernel_dispatch(kernel, X,  X,  eta2, ell);
  arma::mat Knn0  = kernel_dispatch(kernel, X,  X0, eta2, ell);
  arma::mat Kn0n0 = kernel_dispatch(kernel, X0, X0, eta2, ell);

  arma::mat M(n + n0, n + n0);
  M.submat(0, 0, n - 1, n - 1)         = Knn + sigma2 * arma::eye(n, n);
  M.submat(0, n, n - 1, n + n0 - 1)    = Knn0;
  M.submat(n, 0, n + n0 - 1, n - 1)    = Knn0.t();
  M.submat(n, n, n + n0 - 1, n + n0 - 1) = Kn0n0 + (sigma2 / a0) * arma::eye(n0, n0);

  arma::mat ch_M;
  if (!safe_chol_cpp(M, ch_M)) return List::create(_["lml"] = R_NegInf, _["ok"] = false);

  arma::vec yc  = y  - mu;
  arma::vec y0c = y0 - mu;
  arma::vec ystar = arma::join_cols(yc, y0c);
  double lp_data = log_mvn_cpp(ystar, ch_M);

  arma::mat Kn0n0_sym = 0.5 * (Kn0n0 + Kn0n0.t());
  arma::vec D;
  arma::mat U;
  arma::eig_sym(D, U, Kn0n0_sym);
  D.elem(arma::find(D < 0)).zeros();
  arma::vec w = U.t() * y0c;

  double logC = log_C_a0_theta_cpp(a0, D, w, alpha0, beta0, quad_tau, quad_logw);

  return List::create(_["lml"] = lp_data - logC, _["ok"] = true,
                       _["D"] = D, _["w"] = w);
}

// ---------------------------------------------------------------------
// #3: precompute_spectral(), for completeness (shares the same machinery
// as #2; called once per iteration ahead of the a0 update).
// ---------------------------------------------------------------------

// [[Rcpp::export]]
List precompute_spectral_cpp(const arma::mat& X, const arma::mat& X0,
                              double eta2, const arma::vec& ell, std::string kernel,
                              const arma::vec& y, const arma::vec& y0,
                              double sigma2, double mu) {
  arma::mat Knn   = kernel_dispatch(kernel, X,  X,  eta2, ell);
  arma::mat Knn0  = kernel_dispatch(kernel, X,  X0, eta2, ell);
  arma::mat Kn0n0 = kernel_dispatch(kernel, X0, X0, eta2, ell);

  arma::vec yc  = y  - mu;
  arma::vec y0c = y0 - mu;
  const arma::uword n = X.n_rows;

  arma::mat Mnn = Knn + sigma2 * arma::eye(n, n);
  arma::mat ch_Knn;
  if (!safe_chol_cpp(Mnn, ch_Knn))
    return List::create(_["ok"] = false);

  arma::vec alpha_y = arma::solve(arma::trimatu(ch_Knn),
                       arma::solve(arma::trimatl(ch_Knn.t()), yc));
  arma::vec mu_y0_y = Knn0.t() * alpha_y;

  arma::mat V  = arma::solve(arma::trimatu(ch_Knn),
                  arma::solve(arma::trimatl(ch_Knn.t()), Knn0));
  arma::mat Kt = Kn0n0 - Knn0.t() * V;
  Kt = 0.5 * (Kt + Kt.t());

  arma::vec lambda, d;
  arma::mat Q, U;
  arma::eig_sym(lambda, Q, Kt);
  arma::eig_sym(d, U, 0.5 * (Kn0n0 + Kn0n0.t()));
  lambda.elem(arma::find(lambda < 0)).zeros();
  d.elem(arma::find(d < 0)).zeros();

  arma::vec v = Q.t() * (y0c - mu_y0_y);
  arma::vec w = U.t() * y0c;

  return List::create(_["ok"] = true, _["lambda"] = lambda, _["d"] = d,
                       _["v"] = v, _["w"] = w);
}

// ---------------------------------------------------------------------
// #4: predict_gppp()'s per-draw loop -- build kernel mats, conditional
// mean/covariance at Xp, and a joint draw of (mu + f_p), for ALL retained
// posterior draws in one C++ call (avoids one R-level chol/solve chain
// per draw x prediction-target combination).
// ---------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat predict_gppp_draws_cpp(const arma::mat& X, const arma::mat& X0,
                                  const arma::mat& Xp, std::string kernel,
                                  const arma::vec& eta2_draws,
                                  const arma::mat& ell_draws,   // n_save x d
                                  const arma::vec& tau_draws,
                                  const arma::vec& a0_draws,
                                  const arma::vec& mu_draws,
                                  const arma::mat& y_eff,   // n_save x n  (or 1 x n, recycled)
                                  const arma::mat& y0_eff,  // n_save x n0 (or 1 x n0, recycled)
                                  arma::vec rng_seeds_unused) {
  const int n_save = eta2_draws.n_elem;
  const arma::uword np = Xp.n_rows;
  arma::mat F_out(n_save, np, arma::fill::value(NA_REAL));

  const bool y_recycled  = (y_eff.n_rows  == 1);
  const bool y0_recycled = (y0_eff.n_rows == 1);

  for (int i = 0; i < n_save; ++i) {
    double eta2 = eta2_draws[i];
    arma::vec ell = ell_draws.row(i).t();
    double sigma2 = 1.0 / tau_draws[i];
    double a0 = a0_draws[i];
    double mu_i = mu_draws[i];

    arma::vec y_i  = (y_recycled  ? y_eff.row(0).t()  : y_eff.row(i).t());
    arma::vec y0_i = (y0_recycled ? y0_eff.row(0).t() : y0_eff.row(i).t());

    const arma::uword n  = X.n_rows;
    const arma::uword n0 = X0.n_rows;

    arma::mat Knn   = kernel_dispatch(kernel, X,  X,  eta2, ell);
    arma::mat Knn0  = kernel_dispatch(kernel, X,  X0, eta2, ell);
    arma::mat Kn0n0 = kernel_dispatch(kernel, X0, X0, eta2, ell);

    arma::mat M(n + n0, n + n0);
    M.submat(0, 0, n - 1, n - 1)            = Knn + sigma2 * arma::eye(n, n);
    M.submat(0, n, n - 1, n + n0 - 1)       = Knn0;
    M.submat(n, 0, n + n0 - 1, n - 1)       = Knn0.t();
    M.submat(n, n, n + n0 - 1, n + n0 - 1)  = Kn0n0 + (sigma2 / a0) * arma::eye(n0, n0);

    arma::mat ch_M;
    if (!safe_chol_cpp(M, ch_M)) continue;

    arma::vec ystar_c = arma::join_cols(y_i - mu_i, y0_i - mu_i);

    arma::mat Kp_n  = kernel_dispatch(kernel, Xp, X,  eta2, ell);
    arma::mat Kp_n0 = kernel_dispatch(kernel, Xp, X0, eta2, ell);
    arma::mat Kp_star = arma::join_rows(Kp_n, Kp_n0);

    arma::vec alpha = arma::solve(arma::trimatu(ch_M),
                       arma::solve(arma::trimatl(ch_M.t()), ystar_c));
    arma::vec mu_p = Kp_star * alpha + mu_i;

    arma::mat V = arma::solve(arma::trimatu(ch_M),
                   arma::solve(arma::trimatl(ch_M.t()), Kp_star.t()));
    arma::mat Kpp = kernel_dispatch(kernel, Xp, Xp, eta2, ell);
    arma::mat Sigma_p = Kpp - Kp_star * V;
    Sigma_p = 0.5 * (Sigma_p + Sigma_p.t());

    arma::mat ch_Sp;
    if (!safe_chol_cpp(Sigma_p, ch_Sp)) {
      F_out.row(i) = mu_p.t();
      continue;
    }

    arma::vec z = arma::randn(np);
    arma::vec draw = mu_p + ch_Sp.t() * z;
    F_out.row(i) = draw.t();
  }

  return F_out;
}
