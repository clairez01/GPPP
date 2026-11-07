#' Metropolis-within-Gibbs Sampler for the GP Power Prior
#'
#' Compiled hot paths: for the continuous and time-to-event ("full prior",
#' unknown tau) case with a supported kernel (se, matern52, linear), the a0
#' griddy-Gibbs grid (200 points x 40-node tau quadrature) and the collapsed
#' kernel-hyperparameter likelihood evaluation used inside the eta2/ell
#' slice-sampling steps are computed via compiled Rcpp/RcppArmadillo
#' (log_post_a0_grid_cpp / log_post_a0_single_cpp / log_marginal_lik_full_cpp
#' / precompute_spectral_cpp) -- identified as the dominant runtime cost by
#' profiling. This is the path gppp_sims.R exercises. The R-level MCMC
#' control flow (main loop, slice-sampling stepping-out / shrinkage, Gibbs
#' updates for tau/mu, data augmentation) stays in R.
#'
#' Binary outcomes (sigma2_fixed = TRUE) use a logit link with Polya-Gamma
#' augmentation (Section 2.4, revised), discounting the true observed-data
#' logistic likelihood rather than the complete-data augmented-latent
#' likelihood -- see the design note above the "Logit-link Polya-Gamma power
#' prior" block below for why those are not the same prior, and
#' binary_logit.R for the Laplace-approximated normalizing constant that
#' a0's full conditional requires. This path is pure R, not profiled, and out
#' of scope for the compiled path.
#'
#' @inheritParams gpppSampler::gppp_sampler
#' @export
gppp_sampler <- function(
    y, X, y0, X0,
    outcome_type = "continuous",
    delta  = NULL,
    delta0 = NULL,
    kernel    = "se",
    n_iter    = 2000L,
    n_warmup  = 1000L,
    a0_fixed  = NULL,
    init      = list(),
    prior     = list(),
    mh_sd     = list(),
    slice     = FALSE,
    slice_w   = list(),
    n_grid_a0 = 200L,
    n_quad_tau = 40L,
    thin      = 1L,
    verbose   = TRUE
) {
  # ---- Input coercion ----
  y  <- as.vector(y);  X  <- as.matrix(X)
  y0 <- as.vector(y0); X0 <- as.matrix(X0)
  outcome_type <- match.arg(outcome_type, c("continuous", "binary", "tte"))
  stopifnot(length(y) == nrow(X), length(y0) == nrow(X0), ncol(X) == ncol(X0))
  stopifnot(n_warmup < n_iter)
  kernel <- match.arg(kernel, c("se", "matern52", "linear"))

  # ---- Outcome-type validation ----
  if (outcome_type == "binary") {
    if (!all(y  %in% c(0, 1))) stop("binary outcome_type requires y in {0, 1}.")
    if (!all(y0 %in% c(0, 1))) stop("binary outcome_type requires y0 in {0, 1}.")
    sigma2_fixed <- TRUE
  } else if (outcome_type == "tte") {
    if (is.null(delta) || is.null(delta0))
      stop("tte outcome_type requires delta and delta0.")
    delta  <- as.integer(delta)
    delta0 <- as.integer(delta0)
    if (length(delta)  != length(y))  stop("delta must have the same length as y.")
    if (length(delta0) != length(y0)) stop("delta0 must have the same length as y0.")
    if (!all(delta  %in% 0L:1L)) stop("delta must contain only 0s and 1s.")
    if (!all(delta0 %in% 0L:1L)) stop("delta0 must contain only 0s and 1s.")
    sigma2_fixed <- FALSE
  } else {
    sigma2_fixed <- FALSE
  }

  # Compiled fast path only covers the "full prior" (unknown tau) branch.
  use_cpp <- !sigma2_fixed

  d         <- ncol(X)
  n         <- length(y)
  n0        <- length(y0)
  random_a0 <- is.null(a0_fixed)
  if (!random_a0) stopifnot(is.numeric(a0_fixed), a0_fixed > 0, a0_fixed < 1)

  # ---- Defaults ----
  pr <- list(
    log_eta2_mean = 0, log_eta2_sd = 1,
    log_ell_mean  = 0, log_ell_sd  = 1,
    tau_alpha     = 2, tau_beta     = 1,
    mu_mean       = 0, mu_sd        = Inf,
    a0_shape1     = 1, a0_shape2    = 1      # a0 ~ Beta(a0_shape1, a0_shape2)
  )
  pr[names(prior)] <- prior
  alpha0 <- pr$tau_alpha
  beta0  <- pr$tau_beta
  stopifnot(pr$a0_shape1 > 0, pr$a0_shape2 > 0)

  # Beta prior on a0: grid weights for griddy Gibbs (0 when uniform, so the
  # default grid draw is unchanged), log density for the slice paths
  a0_unif    <- pr$a0_shape1 == 1 && pr$a0_shape2 == 1
  a0_grid_lp <- if (a0_unif) 0 else log_prior_a0_grid(n_grid_a0, pr$a0_shape1, pr$a0_shape2)
  lp_a0_prior <- function(a0_) log_prior_a0(a0_, pr$a0_shape1, pr$a0_shape2)

  mh <- list(log_eta2 = 0.3, log_ell = 0.3)
  mh[names(mh_sd)] <- mh_sd

  # ---- Slice sampling setup ----
  valid_sp  <- c("eta2", "ell", "a0")
  use_slice <- if (isTRUE(slice))           c("eta2", "ell")
               else if (is.character(slice)) intersect(slice, valid_sp)
               else                          character(0)

  sw <- list(log_eta2 = 1.0, log_ell = 1.0, logit_a0 = 1.0)
  sw[names(slice_w)] <- slice_w

  iv <- list(eta2 = 1, ell = 1, sigma2 = 1, a0 = 0.5, mu = NULL)
  iv[names(init)] <- init
  if (!random_a0) iv$a0 <- a0_fixed
  if (is.null(iv$mu)) iv$mu <- mean(y)

  if (length(iv$ell) == 1L) iv$ell <- rep(iv$ell, d)
  stopifnot(length(iv$ell) == d)

  pr$log_ell_mean <- rep_len(pr$log_ell_mean, d)
  pr$log_ell_sd   <- rep_len(pr$log_ell_sd,   d)

  quad <- if (!sigma2_fixed) make_tau_quadrature(n_quad_tau, alpha0, beta0) else NULL

  # ---- Log-prior helpers ----
  lp_log_eta2 <- function(le2) {
    dnorm(le2, pr$log_eta2_mean, pr$log_eta2_sd, log = TRUE)
  }
  lp_log_ell_j <- function(ll_j, j) {
    dnorm(ll_j, pr$log_ell_mean[j], pr$log_ell_sd[j], log = TRUE)
  }
  lp_theta <- function(le2, ll_vec) {
    lp_log_eta2(le2) +
    sum(dnorm(ll_vec, pr$log_ell_mean, pr$log_ell_sd, log = TRUE))
  }

  # ---- Storage ----
  n_save <- floor((n_iter - n_warmup) / thin)
  samp <- list(
    eta2   = numeric(n_save),
    ell    = matrix(NA_real_, nrow = n_save, ncol = d),
    tau    = numeric(n_save),
    a0     = numeric(n_save),
    mu     = numeric(n_save),
    log_ml = numeric(n_save),
    f      = if (outcome_type == "binary") matrix(NA_real_, n_save, n)  else NULL,
    f0     = if (outcome_type == "binary") matrix(NA_real_, n_save, n0) else NULL,
    y_aug  = if (outcome_type == "tte") matrix(NA_real_, n_save, length(y))  else NULL,
    y0_aug = if (outcome_type == "tte") matrix(NA_real_, n_save, length(y0)) else NULL
  )
  acc <- c(eta2 = 0L, ell = 0L)

  # ---- Initialize augmented observations (continuous/TTE only) ----
  y_aug  <- y
  y0_aug <- y0

  # ---- Initialise parameters ----
  eta2   <- iv$eta2
  ell    <- iv$ell
  sigma2 <- 1
  tau    <- if (sigma2_fixed) 1 else 1 / iv$sigma2
  a0     <- iv$a0
  mu     <- iv$mu
  f      <- if (outcome_type == "binary") numeric(n)  else NULL
  f0     <- if (outcome_type == "binary") numeric(n0) else NULL
  Km     <- build_kernel_mats(X, X0, eta2, ell, kernel)

  # Collapsed marginal log-likelihood at given (eta2, ell, a0, mu). For the
  # binary/logit branch, this is the Gaussian working-response marginal
  # collapsed over f given the CURRENT Polya-Gamma weights (aug$omega,
  # aug$omega0) -- exact given those weights, no normalizing-constant
  # approximation needed here (that only enters a0's own full conditional
  # below, via log_post_a0_laplace). a0_ is ignored: a0 already entered
  # through aug$z0/aug$omega0 when Step 0 drew the augmentation this sweep.
  # For the "full prior" (sigma2_fixed = FALSE) branch with a supported
  # kernel, this is the compiled fast path -- it rebuilds the kernel
  # matrices internally in C++ rather than taking a precomputed Km, so it
  # can be called directly from slice-sampling proposal closures without an
  # R-level build_kernel_mats() allocation per evaluation.
  lml_fn <- function(eta2_, ell_, a0_, mu_) {
    if (outcome_type == "binary") {
      Km_     <- build_kernel_mats(X, X0, eta2_, ell_, kernel)
      z_star  <- c(aug$z - mu_, aug$z0 - mu_)
      var_vec <- c(1 / aug$omega, 1 / aug$omega0)
      log_marginal_lik_hetero(z_star, Km_, var_vec)
    } else if (use_cpp) {
      res <- log_marginal_lik_full_cpp(X, X0, eta2_, ell_, kernel, y_aug, y0_aug,
                                        tau, a0_, mu_, alpha0, beta0, quad$tau, quad$log_w)
      if (!isTRUE(res$ok)) -Inf else res$lml
    } else {
      Km_ <- build_kernel_mats(X, X0, eta2_, ell_, kernel)
      log_marginal_lik_full(y_aug, y0_aug, Km_, tau, a0_, mu_, alpha0, beta0, quad)
    }
  }
  # Placeholder for binary: lml_fn() needs aug$z/omega, drawn by Step 0 of
  # iteration 1 before lml is used again.
  lml <- if (outcome_type == "binary") 0 else lml_fn(eta2, ell, a0, mu)

  si <- 0L

  # ---- Main loop ----
  for (iter in seq_len(n_iter)) {

    # === Step 0: Data augmentation (non-continuous outcomes only) ===
    if (outcome_type == "binary") {
      aug <- augment_binary_pg(mu + f, mu + f0, y, y0, a0)
      lml <- lml_fn(eta2, ell, a0, mu)
    } else if (outcome_type == "tte") {
      aug    <- augment_tte(y_aug, y0_aug, y, y0, delta, delta0, Km, 1 / tau, a0, mu)
      y_aug  <- aug$T_aug
      y0_aug <- aug$T0_aug
      lml    <- lml_fn(eta2, ell, a0, mu)
    }

    # === Step 1: Update eta2 ===
    if ("eta2" %in% use_slice) {
      le2_c     <- log(eta2)
      lp_cur_e2 <- lml + lp_theta(le2_c, log(ell))
      lp_e2 <- function(le2) {
        lml_p <- lml_fn(exp(le2), ell, a0, mu)
        if (!is.finite(lml_p)) return(-Inf)
        lml_p + lp_theta(le2, log(ell))
      }
      res  <- slice_sample_1d(le2_c, lp_cur_e2, lp_e2, w = sw$log_eta2)
      eta2 <- exp(res$x)
      lml  <- res$lp - lp_theta(res$x, log(ell))
    } else {
      le2_c <- log(eta2)
      le2_p <- le2_c + rnorm(1L, 0, mh$log_eta2)
      lml_p <- lml_fn(exp(le2_p), ell, a0, mu)
      log_r <- (lml_p - lml) + lp_log_eta2(le2_p) - lp_log_eta2(le2_c)
      if (is.finite(log_r) && log(runif(1L)) < log_r) {
        eta2 <- exp(le2_p); lml <- lml_p
        acc["eta2"] <- acc["eta2"] + 1L
      }
    }

    # === Step 2: Update each ell_j (ARD — one dimension at a time) ===
    if ("ell" %in% use_slice) {
      for (j in seq_len(d)) {
        ll_j_c   <- log(ell[j])
        lp_cur_j <- lml + lp_log_ell_j(ll_j_c, j)
        lp_ll_j  <- local({
          jj <- j
          function(ll_j) {
            ell_p    <- ell; ell_p[jj] <- exp(ll_j)
            lml_p    <- lml_fn(eta2, ell_p, a0, mu)
            if (!is.finite(lml_p)) return(-Inf)
            lml_p + lp_log_ell_j(ll_j, jj)
          }
        })
        res    <- slice_sample_1d(ll_j_c, lp_cur_j, lp_ll_j, w = sw$log_ell)
        ell[j] <- exp(res$x)
        lml    <- res$lp - lp_log_ell_j(res$x, j)
      }
    } else {
      for (j in seq_len(d)) {
        ll_j_c <- log(ell[j])
        ll_j_p <- ll_j_c + rnorm(1L, 0, mh$log_ell)
        ell_p  <- ell; ell_p[j] <- exp(ll_j_p)
        lml_p  <- lml_fn(eta2, ell_p, a0, mu)
        log_r  <- (lml_p - lml) + lp_log_ell_j(ll_j_p, j) - lp_log_ell_j(ll_j_c, j)
        if (is.finite(log_r) && log(runif(1L)) < log_r) {
          ell   <- ell_p
          lml   <- lml_p
          acc["ell"] <- acc["ell"] + 1L
        }
      }
    }

    # Kernel matrices at the now-finalised eta2/ell, needed below by
    # draw_gp_posterior() / update_mu_gibbs() / next iteration's
    # augmentation step (not yet ported to the compiled path).
    Km <- build_kernel_mats(X, X0, eta2, ell, kernel)

    if (sigma2_fixed) {
      # ===== Logit-link Polya-Gamma power prior (Section 2.4, revised) =====
      # a0 discounts the true observed-data logistic likelihood, not the
      # complete-data augmented-latent likelihood. Theta was already updated
      # above, collapsed over f given the current Polya-Gamma weights
      # (exact, no approximation). f* is drawn here at the finalised theta,
      # then mu, then a0. a0's full conditional needs the Laplace-
      # approximated normalizing constant C(a0, theta): the Polya-Gamma
      # augmentation gives exact f0 draws but not this normalizer -- it does
      # not cancel because C depends on theta as well as a0 (see
      # binary_logit.R).

      z_star  <- c(aug$z - mu, aug$z0 - mu)
      var_vec <- c(1 / aug$omega, 1 / aug$omega0)
      f_star  <- draw_gp_posterior_hetero(z_star, Km, var_vec)
      if (!is.null(f_star)) {
        f  <- f_star[seq_len(n)]
        f0 <- f_star[seq_len(n0) + n]
      }

      mu <- update_mu_gibbs_pg(aug$z, aug$z0, f, f0, aug$omega, aug$omega0,
                                pr$mu_mean, pr$mu_sd^2)

      if (random_a0) {
        K0 <- Km$Kn0n0
        if ("a0" %in% use_slice) {
          f0_hat_warm <- f0
          la0_c <- log(a0 / (1 - a0))
          cur   <- log_post_a0_laplace(a0, f0, y0, K0, mu, f0_init = f0_hat_warm)
          f0_hat_warm <- cur$f0_hat
          lp_cur_a0 <- cur$lp + lp_a0_prior(a0) + log(a0) + log(1 - a0)
          lp_a0 <- function(la0) {
            a0_ <- 1 / (1 + exp(-la0))
            res <- log_post_a0_laplace(a0_, f0, y0, K0, mu, f0_init = f0_hat_warm)
            f0_hat_warm <<- res$f0_hat
            if (!is.finite(res$lp)) return(-Inf)
            res$lp + lp_a0_prior(a0_) + log(a0_) + log(1 - a0_)
          }
          res <- slice_sample_1d(la0_c, lp_cur_a0, lp_a0, w = sw$logit_a0)
          a0  <- 1 / (1 + exp(-res$x))
        } else {
          a0 <- update_a0_logit_griddy(f0, y0, K0, mu, n_grid_a0, log_prior = a0_grid_lp)
        }
      }

      lml <- a0 * log_lik_logit_sum(mu + f0, y0) + log_lik_logit_sum(mu + f, y)

    } else {
      # ===== Full power prior (unknown tau): partially-collapsed sweep =====
      # a0 update and the kernel-hyperparameter likelihoods above route
      # through the compiled C++ path (log_marginal_lik_full_cpp /
      # precompute_spectral_cpp / log_post_a0_grid_cpp /
      # log_post_a0_single_cpp) for a supported kernel.

      if (random_a0) {
        sp <- if (use_cpp) {
          precompute_spectral_cpp(X, X0, eta2, ell, kernel, y_aug, y0_aug, 1 / tau, mu)
        } else {
          precompute_spectral(y_aug, y0_aug, Km, 1 / tau, mu)
        }
        if (isTRUE(sp$ok) || (!use_cpp && !is.null(sp))) {
          if ("a0" %in% use_slice) {
            lp_a0_val <- function(a0_) {
              if (use_cpp) {
                log_post_a0_single_cpp(a0_, tau, sp$lambda, sp$v, sp$d, sp$w,
                                        alpha0, beta0, quad$tau, quad$log_w)
              } else {
                log_post_a0_full(a0_, tau, sp, alpha0, beta0, quad)
              }
            }
            la0_c     <- log(a0 / (1 - a0))
            lp_cur_a0 <- lp_a0_val(a0) + lp_a0_prior(a0) + log(a0) + log(1 - a0)
            lp_a0 <- function(la0) {
              a0_ <- 1 / (1 + exp(-la0))
              lp  <- lp_a0_val(a0_)
              if (!is.finite(lp)) return(-Inf)
              lp + lp_a0_prior(a0_) + log(a0_) + log(1 - a0_)
            }
            res <- slice_sample_1d(la0_c, lp_cur_a0, lp_a0, w = sw$logit_a0)
            a0  <- 1 / (1 + exp(-res$x))
          } else {
            a0_grid <- seq(1e-6, 1 - 1e-6, length.out = n_grid_a0)
            log_w <- if (use_cpp) {
              log_post_a0_grid_cpp(a0_grid, tau, sp$lambda, sp$v, sp$d, sp$w,
                                    alpha0, beta0, quad$tau, quad$log_w)
            } else {
              vapply(a0_grid, log_post_a0_full, numeric(1L),
                     tau = tau, sp = sp, alpha0 = alpha0, beta0 = beta0, quad = quad)
            }
            log_w <- log_w + a0_grid_lp
            log_w <- log_w - max(log_w)
            wts   <- exp(log_w); wts <- wts / sum(wts)
            a0    <- sample(a0_grid, 1L, prob = wts)
          }
          lml <- lml_fn(eta2, ell, a0, mu)
        }
      }

      # === Draw f* = (f, f0) jointly, then tau (Gibbs) ===
      f_star <- draw_gp_posterior(c(y_aug - mu, y0_aug - mu), Km, 1 / tau, a0)
      if (!is.null(f_star)) {
        f   <- f_star[seq_len(n)]
        f0  <- f_star[seq_len(n0) + n]
        tau <- update_tau_gibbs(y_aug, y0_aug, f, f0, a0, mu, alpha0, beta0)
      }

      # === Update mu (conjugate Gibbs off the retained f*) ===
      if (!is.null(f)) {
        mu  <- update_mu_gibbs_full(y_aug, y0_aug, f, f0, tau, a0, pr$mu_mean, pr$mu_sd^2)
        lml <- lml_fn(eta2, ell, a0, mu)
      }
    }

    # === Save post-warmup ===
    if (iter > n_warmup && (iter - n_warmup) %% thin == 0L) {
      si <- si + 1L
      samp$eta2[si]     <- eta2
      samp$ell[si, ]    <- ell
      samp$tau[si]      <- tau
      samp$a0[si]       <- a0
      samp$mu[si]       <- mu
      samp$log_ml[si]   <- lml
      if (!is.null(samp$f))     samp$f[si, ]      <- f
      if (!is.null(samp$f0))    samp$f0[si, ]     <- f0
      if (!is.null(samp$y_aug))  samp$y_aug[ si, ] <- y_aug
      if (!is.null(samp$y0_aug)) samp$y0_aug[si, ] <- y0_aug
    }

    if (verbose && iter %% 500L == 0L) {
      ell_str <- if (d == 1L) sprintf("%.3f", ell) else
                 sprintf("[%s]", paste(sprintf("%.3f", ell), collapse = ","))
      message(sprintf(
        "[%d/%d]  eta2=%.3f  ell=%s  tau=%.3f  a0=%.3f  mu=%.3f  lml=%.2f",
        iter, n_iter, eta2, ell_str, tau, a0, mu, lml
      ))
    }
  }

  vec_names <- setdiff(names(samp), c("ell", "f", "f0", "y_aug", "y0_aug"))
  for (nm in vec_names) if (!is.null(samp[[nm]])) samp[[nm]] <- samp[[nm]][seq_len(si)]
  samp$ell <- samp$ell[seq_len(si), , drop = FALSE]
  if (!is.null(samp$f))      samp$f      <- samp$f[     seq_len(si), , drop = FALSE]
  if (!is.null(samp$f0))     samp$f0     <- samp$f0[    seq_len(si), , drop = FALSE]
  if (!is.null(samp$y_aug))  samp$y_aug  <- samp$y_aug[ seq_len(si), , drop = FALSE]
  if (!is.null(samp$y0_aug)) samp$y0_aug <- samp$y0_aug[seq_len(si), , drop = FALSE]

  acc_rates <- c(
    eta2 = as.numeric(acc["eta2"]) / n_iter,
    ell  = as.numeric(acc["ell"])  / (n_iter * d),
    tau  = NA_real_
  )
  for (p in intersect(use_slice, c("eta2", "ell")))
    acc_rates[p] <- NA_real_

  structure(
    list(
      samples          = samp,
      acceptance_rates = acc_rates,
      kernel           = kernel,
      outcome_type     = outcome_type,
      sigma2_fixed     = sigma2_fixed,
      d                = d,
      n_iter           = n_iter,
      n_warmup         = n_warmup,
      thin             = thin,
      random_a0        = random_a0,
      a0_fixed         = a0_fixed,
      prior            = pr,
      mh_sd            = mh,
      use_slice        = use_slice,
      slice_w          = sw
    ),
    class = "gppp_fit"
  )
}
