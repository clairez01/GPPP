# gpppSampler

R implementation of a collapsed Metropolis-within-Gibbs sampler for the **Gaussian Process Power Prior (GPPP)**. The model borrows strength from a historical dataset by discounting its likelihood by a parameter a0 ∈ (0, 1) before updating a GP prior with current data. All derivations follow Zhu (2026), *Gaussian Process Power Prior*.

---

## Quick start

```r
# Source all package files
pkg_dir <- "path/to/gpppSampler/R"
for (f in list.files(pkg_dir, full.names = TRUE)) source(f)

fit      <- gppp_sampler(y, X, y0, X0)
F_draws  <- predict_gppp(fit, y, X, y0, X0, Xp = Xp_new)
```

---

## `gppp_sampler()`

Runs the sampler and returns posterior draws of the kernel hyperparameters, noise variance, and (optionally) the borrowing parameter a0. The latent function f\* is integrated out analytically — only hyperparameters are stored.

### Sampler algorithm (Section 3.8)

Each iteration (default — random-walk MH / griddy-Gibbs):
1. Update η² via random-walk Metropolis on log η²
2. Update ℓ via random-walk Metropolis on log ℓ
3. Update σ² via random-walk Metropolis on log σ²
4. Update a0 via griddy-Gibbs on a discrete grid *(random a0 only)*

Steps 1–3 can optionally use **univariate slice sampling** (Neal 2003) instead of RW-MH; a0 can optionally use slice sampling on logit a0. See the **Slice sampling** section below.

The common target for steps 1–3 is the marginal likelihood (Eq. 20):

```
log p(y | y0, a0, σ², θ) = log N(y* | 0, M) − log N(y0 | 0, M0)
```

where M and M0 are the joint and historical covariance matrices (Eq. 19), and y\* = c(y, y0).

---

### Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `y` | numeric vector | — | Outcomes for the **current** dataset. Length n. |
| `X` | numeric matrix | — | Inputs for the current dataset. Dimensions n × d. Each row is one observation; columns are covariates (e.g., `cbind(dose, z1, z2)`). Must have the same number of columns as `X0`. |
| `y0` | numeric vector | — | Outcomes for the **historical** dataset. Length n0. |
| `X0` | numeric matrix | — | Inputs for the historical dataset. Dimensions n0 × d. Columns must match `X`. |
| `n_iter` | integer | `2000` | Total number of MCMC iterations, including warm-up. |
| `n_warmup` | integer | `1000` | Number of iterations discarded as warm-up (burn-in). Must be less than `n_iter`. Saved samples = floor((n_iter − n_warmup) / thin). |
| `a0_fixed` | numeric or `NULL` | `NULL` | Controls the borrowing parameter a0. If `NULL`, a0 is treated as unknown and updated at every iteration under a Uniform(0, 1) prior (recommended). If a number in (0, 1), a0 is held fixed at that value for the entire run — e.g., `a0_fixed = 1` borrows the historical data fully, `a0_fixed = 0.5` applies 50% discounting. |
| `init` | named list | see below | Starting values for the sampler. Recognised names: `eta2`, `ell`, `sigma2`, `a0`. Unspecified entries default to `eta2 = 1`, `ell = 1`, `sigma2 = 1`, `a0 = 0.5`. |
| `prior` | named list | see below | Prior hyperparameters. See **Priors** section below. |
| `mh_sd` | named list | see below | Standard deviations of the random-walk Metropolis proposals, on the log scale. See **Tuning** section below. Ignored for parameters that use slice sampling. |
| `slice` | logical or character | `FALSE` | Selects univariate slice sampling instead of RW-MH for one or more parameters. `TRUE` enables slice for `eta2`, `ell`, and `sigma2`; a character vector such as `c("eta2","ell","sigma2","a0")` gives per-parameter control. See **Slice sampling** section below. |
| `slice_w` | named list | see below | Initial stepping-out widths for the slice sampler, on the log or logit scale. Recognised names: `log_eta2`, `log_ell`, `log_sigma2`, `logit_a0`. All default to `1.0`. |
| `n_grid_a0` | integer | `200` | Number of equally-spaced grid points in (0, 1) used for the griddy-Gibbs a0 update. Ignored when `"a0"` is in `slice`. |
| `thin` | integer | `1` | Thinning interval applied after warm-up. Every `thin`-th post-warmup iterate is saved. `thin = 1` keeps all; `thin = 10` keeps every 10th. |
| `verbose` | logical | `TRUE` | If `TRUE`, prints a one-line progress update every 500 iterations showing the current parameter values and log marginal likelihood. |

---

### Priors

Passed as a named list to the `prior` argument. Any subset of names can be supplied; the rest use their defaults.

#### Kernel hyperparameters θ = (η², ℓ)

The squared exponential kernel is k(r) = η² exp(−r²/2ℓ²). Both parameters are given **log-normal** priors, meaning their logarithms follow normal distributions. Sampling is done on the log scale.

| Name | Default | Description |
|---|---|---|
| `log_eta2_mean` | `0` | Mean of the Normal prior on log(η²). Equivalently, the log of the prior median of η². |
| `log_eta2_sd` | `1` | Standard deviation of the Normal prior on log(η²). Larger values = more diffuse prior on the signal variance. |
| `log_ell_mean` | `0` | Mean of the Normal prior on log(ℓ). Setting this to log(your expected length-scale) is a reasonable starting point. |
| `log_ell_sd` | `1` | Standard deviation of the Normal prior on log(ℓ). |

**Guidance:** If your inputs are scaled to [0, 1], a prior centered near log(0.3) for ℓ (i.e., `log_ell_mean = log(0.3) ≈ −1.2`) is often reasonable. If inputs span a wider range, increase `log_ell_mean` accordingly.

#### Noise precision τ

τ (= 1/σ²) follows a **gamma** prior in shape-rate form: τ ~ Gamma(α, β), with mean α/β. This is numerically the same specification as the equivalent inverse-gamma prior σ² ~ IG(α, β), with mean β/(α − 1) for α > 1.

| Name | Default | Description |
|---|---|---|
| `tau_alpha` | `2` | Shape parameter α of the Gamma(τ) prior. Must be > 0. Values near 2 give a weakly informative prior. |
| `tau_beta` | `1` | Rate parameter β of the Gamma(τ) prior. The prior mean of τ is α/β (equivalently, prior mean of σ² is β/(α − 1)). Set `tau_beta` close to your expected noise variance to centre the prior. |

**Example:** If you expect σ² ≈ 0.25, set `tau_alpha = 3, tau_beta = 0.5` (implied prior mean of σ² = 0.5/2 = 0.25).

---

### Tuning (`mh_sd`)

Standard deviations of the Gaussian random-walk proposals on the log scale. Target MH acceptance rates of 20–40%. These are ignored for any parameter that uses slice sampling.

| Name | Default | Description |
|---|---|---|
| `log_eta2` | `0.3` | Proposal SD for log(η²). Increase if acceptance rate for η² is too high; decrease if too low. |
| `log_ell` | `0.3` | Proposal SD for log(ℓ). |
| `log_sigma2` | `0.3` | Proposal SD for log(σ²). |

---

### Slice sampling

Slice sampling (Neal 2003) replaces the RW-MH accept/reject step with a bracketing procedure that always accepts. The benefit is no step-size tuning and good mixing across multimodal or highly curved conditionals; the trade-off is more log-density evaluations per iteration (typically 3–10 per parameter vs. 1 for MH).

**How it works (stepping-out + shrinkage):**
1. Draw a uniform level below the current log-density: `log u ~ Uniform(-∞, log p(x_cur))`.
2. Expand a random interval `[L, R]` around `x_cur` in steps of width `w` until both endpoints are below the level.
3. Sample `x_new ~ Uniform(L, R)`, shrinking the interval on failures until `log p(x_new) > log u`.

Each univariate slice update is run on the log / logit scale matching the existing MH updates:

| Parameter | Scale | `slice_w` key |
|---|---|---|
| η² | log η² | `log_eta2` |
| ℓ  | log ℓ  | `log_ell`  |
| σ² | log σ² | `log_sigma2` |
| a0 | logit a0 | `logit_a0` |

**Slice widths (`slice_w`):** The initial stepping-out width `w` affects efficiency but not correctness. If the slice width is too narrow relative to the posterior, more stepping-out expansions are needed; if too wide, shrinkage takes longer. All widths default to `1.0`; values between `0.5` and `2.0` are usually reasonable on the log/logit scale.

**Acceptance rates:** Slice-sampled parameters appear as `NA` in `fit$acceptance_rates` — the algorithm always accepts, so the rate is not meaningful.

**a0 and griddy-Gibbs:** The griddy-Gibbs update for a0 is O(n0) and produces exact draws from the discrete approximation to the conditional, making it highly efficient. Slice sampling on logit a0 is a valid alternative but is generally slower and only worth choosing if you want a fully gradient-free, tuning-free setup.

#### Usage examples

```r
# Default: RW-MH for all parameters
fit_mh <- gppp_sampler(y, X, y0, X0)

# Slice sampling for all three MH parameters; a0 keeps griddy-Gibbs
fit_slice <- gppp_sampler(y, X, y0, X0, slice = TRUE)

# Slice only for ell (useful if ell mixes poorly under MH)
fit_sl_ell <- gppp_sampler(y, X, y0, X0, slice = "ell")

# Slice for everything including a0
fit_full_slice <- gppp_sampler(y, X, y0, X0,
                               slice = c("eta2", "ell", "sigma2", "a0"))

# Slice with wider stepping-out for ell (useful for longer length-scales)
fit_wide <- gppp_sampler(y, X, y0, X0,
                         slice   = TRUE,
                         slice_w = list(log_ell = 2.0))

# Mix: slice for eta2 and ell; tuned MH for sigma2; griddy-Gibbs for a0
fit_mixed <- gppp_sampler(y, X, y0, X0,
                          slice   = c("eta2", "ell"),
                          mh_sd   = list(log_sigma2 = 0.5))
```

---

### Return value

A list of class `"gppp_fit"` with the following elements.

| Element | Description |
|---|---|
| `samples$eta2` | Numeric vector of posterior draws of η² (signal variance). Length = floor((n_iter − n_warmup) / thin). |
| `samples$ell` | Posterior draws of ℓ (length-scale). |
| `samples$sigma2` | Posterior draws of σ² (noise variance). |
| `samples$a0` | Posterior draws of a0. Fixed at `a0_fixed` for every iteration if `a0_fixed` was set. |
| `samples$log_ml` | Log marginal likelihood log p(y \| y0, a0, σ², θ) at each saved iterate. Useful for diagnosing convergence. |
| `acceptance_rates` | Named vector of MH acceptance rates for `eta2`, `ell`, and `sigma2` over all `n_iter` iterations. Parameters updated via slice sampling appear as `NA`. |
| `random_a0` | Logical; `TRUE` if a0 was treated as random. |
| `a0_fixed` | The value of `a0_fixed` passed in, or `NULL`. |
| `prior` | The resolved prior list (defaults filled in). |
| `mh_sd` | The resolved proposal SD list. |
| `use_slice` | Character vector naming the parameters that used slice sampling (empty if `slice = FALSE`). |
| `slice_w` | The resolved slice-width list (all four keys, with defaults filled in). |
| `n_iter`, `n_warmup`, `thin` | Sampler settings as passed. |

Use `print(fit)` for a summary of posterior means and acceptance rates, and `summary(fit)` for a quantile table.

---

## `predict_gppp()`

For each saved MCMC sample of (η², ℓ, σ², a0), draws the latent function f jointly at a set of prediction inputs Xp from its posterior predictive distribution (Section 3.6, Eq. 26):

```
f_p | y*, a0, σ², θ  ~  N(μ_p, Σ_p)

μ_p    = K_{p*} (K + Σ)⁻¹ y*
Σ_p    = K_{pp} − K_{p*} (K + Σ)⁻¹ K_{*p}
```

where K_{p\*} = [K_{pn}, K_{pn0}] is the cross-covariance between the prediction points and **all** data (current and historical), and y\* = c(y, y0). The draw is **joint** across all rows of Xp, which preserves posterior correlations needed when computing contrasts such as f(dose=1, z) − f(dose=0, z) for the Bayesian bootstrap.

---

### Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `fit` | `gppp_fit` | — | The object returned by `gppp_sampler()`. |
| `y` | numeric vector | — | Current outcomes. Must be identical to what was passed to `gppp_sampler()`. |
| `X` | numeric matrix | — | Current inputs. Must be identical to what was passed to `gppp_sampler()`. |
| `y0` | numeric vector | — | Historical outcomes. Must be identical to what was passed to `gppp_sampler()`. |
| `X0` | numeric matrix | — | Historical inputs. Must be identical to what was passed to `gppp_sampler()`. |
| `Xp` | numeric matrix or `NULL` | `X` | Prediction inputs, n_p × d. These do not need to be the observed data points — they can be any inputs in the same space as X and X0. If `NULL`, defaults to the current inputs X. See the **Prediction inputs** section below for common patterns. |

---

### Return value

An **n_save × n_p** numeric matrix where n_save = `length(fit$samples$eta2)` and n_p = `nrow(Xp)`.

- **Row i** is a single joint draw of f at all n_p prediction points, under the i-th posterior hyperparameter sample.
- **Column j** corresponds to f evaluated at `Xp[j, ]`.

If the posterior covariance Σ_p is numerically non-positive-definite for a given iteration (rare), that row falls back to the posterior mean μ_p rather than a stochastic draw.

---

### Prediction inputs

The choice of Xp depends on the estimand.

**Posterior mean function** (fitted values at observed points):

```r
F_draws <- predict_gppp(fit, y, X, y0, X0)   # Xp defaults to X
f_mean  <- colMeans(F_draws)                  # posterior mean of f at each xi
```

**Average treatment effect via Bayesian bootstrap** (dose=1 vs dose=0 for current subjects):

```r
z_obs      <- X[, -dose_col, drop = FALSE]        # covariate columns only
Xp_treated <- cbind(dose = 1, z_obs)
Xp_control <- cbind(dose = 0, z_obs)
Xp_all     <- rbind(Xp_treated, Xp_control)       # 2n x d

F_draws <- predict_gppp(fit, y, X, y0, X0, Xp = Xp_all)
# cols 1:n     = f(dose=1, z_i)
# cols (n+1):2n = f(dose=0, z_i)

n <- nrow(z_obs)
ate_draws <- apply(F_draws, 1, function(fp) {
  w <- diff(c(0, sort(runif(n - 1)), 1))     # Dirichlet(1,...,1) weights
  sum(w * (fp[seq_len(n)] - fp[n + seq_len(n)]))
})
quantile(ate_draws, c(0.025, 0.5, 0.975))
```

**Dose-response curve** (fine grid of doses):

```r
doses   <- seq(0, 1, length.out = 50)
z_fixed <- rep(mean(X[, "z"]), 50)          # fix covariate at its mean
Xp_grid <- cbind(dose = doses, z = z_fixed)

F_draws  <- predict_gppp(fit, y, X, y0, X0, Xp = Xp_grid)
f_mean   <- colMeans(F_draws)
f_lower  <- apply(F_draws, 2, quantile, 0.025)
f_upper  <- apply(F_draws, 2, quantile, 0.975)
```

---

## Notes on computational cost

`predict_gppp()` loops over all n_save retained samples. Each iteration requires one (n + n0) × (n + n0) Cholesky factorisation (O((n + n0)³)) and one n_p × n_p Cholesky factorisation (O(n_p³)). For large datasets or many prediction points, consider subsampling the posterior before calling:

```r
idx     <- sample(length(fit$samples$eta2), 200)
fit_sub <- fit
fit_sub$samples <- lapply(fit$samples, `[`, idx)
F_draws <- predict_gppp(fit_sub, y, X, y0, X0, Xp = Xp_all)
```
