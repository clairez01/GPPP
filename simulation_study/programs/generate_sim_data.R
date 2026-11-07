# ------------------------------------------------------------------------------
# Program name: generate_sim_data.R
# Description.: Generate current + historical datasets for the GPPP simulation
#               study, for every scenario in grid.rds
# ------------------------------------------------------------------------------

library(haven)

## Directory to store files
paper.dir <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3"
save.dir  <- file.path(paper.dir, "Simulation_Study", "Inputs")

# Load grid of simulation scenarios
grid <- readRDS(file.path(save.dir, "grid.rds"))

# Load in CSL data
csl.dir    <- file.path(paper.dir, "Datasets")
current_real    <- read_sas(file.path(csl.dir, "csl830_current.sas7bdat"))
historical_real <- read_sas(file.path(csl.dir, "csl830_external.sas7bdat"))

current_real$DOSE_scaled    <- current_real$DOSE / 1000
historical_real$DOSE_scaled <- historical_real$DOSE / 1000

n_reps <- 10000   # Monte Carlo repetitions per scenario (matches n_reps_total in create_sim_grid.R)
n_curr <- 80      # Current trial sample size (fixed, matches grid$n)

#-------------------------------------------------------------------------------
# --- True Data-Generating Process ---
#-------------------------------------------------------------------------------
# True dose-response: mu0 + f(dose) + beta_z * z, f(dose) is linear, emax, or
# (Scenario B) a steep sigmoid. mu0/beta_z/TE_target/sd_truth are fit once from
# the real current-trial data and shared across all three; ED50 and h_B (curve
# location/steepness for emax and B) are hand-set since they aren't identifiable
# from this data (81 subjects, dose range 0-9.4).
fit_lin <- lm(CHG ~ DOSE_scaled + BASE, data = current_real)
print(summary(fit_lin))

mu0       <- unname(coef(fit_lin)["(Intercept)"])    # Baseline (placebo) mean response
beta_z    <- unname(coef(fit_lin)["BASE"])            # Covariate effect
TE_target <- unname(coef(fit_lin)["DOSE_scaled"]) * 5 # Treatment effect at dose = 5 vs dose = 0
sd_truth  <- sigma(fit_lin)                            # Residual SD
ED50      <- 5      # Half-maximal effective dose (emax truth)
h_B       <- 4      # Hill coefficient (Scenario B steepness)

# --- True mean functions (vectorized over dose/z; beta is a named list) ---
mean_fun_linear <- function(dose, z, beta) {
  beta$mu0 + beta$slope * dose + beta$beta_z * z
}
mean_fun_emax <- function(dose, z, beta) {
  beta$mu0 + (beta$Emax * dose) / (beta$ED50 + dose) + beta$beta_z * z
}

# Scenario B: Hill model with h = h_B, vs. emax's h = 1 - same mu0/ED50/Emax,
# different curvature.
mean_fun_B <- function(dose, z, beta) {
  beta$mu0 + (beta$Emax * dose^beta$h) / (beta$ED50^beta$h + dose^beta$h) + beta$beta_z * z
}

# --- True parameters per dose-response family ---
get_true_beta <- function(dose.resp) {
  if (dose.resp == "linear") {
    list(mu0 = mu0, beta_z = beta_z, slope = TE_target / 5)
  } else if (dose.resp == "emax") {
    list(mu0 = mu0, beta_z = beta_z, Emax = 2 * TE_target, ED50 = ED50)
  } else {
    stop("Unknown dose.resp: ", dose.resp)
  }
}

get_mean_fun <- function(dose.resp) {
  if (dose.resp == "linear") mean_fun_linear
  else if (dose.resp == "emax") mean_fun_emax
  else stop("Unknown dose.resp: ", dose.resp)
}

beta_lin  <- get_true_beta("linear")
beta_emax <- get_true_beta("emax")
beta_B    <- list(mu0 = mu0, beta_z = beta_z, Emax = 2 * TE_target, ED50 = ED50, h = h_B)

# --- Plot true dose-response shapes (at z = 0) ---
dose_seq  <- seq(0, 10, length.out = 200)
ylim_all  <- range(mean_fun_linear(dose_seq, 0, beta_lin),
                    mean_fun_emax(dose_seq, 0, beta_emax),
                    mean_fun_B(dose_seq, 0, beta_B))

curve(mean_fun_linear(x, 0, beta_lin), from = 0, to = 10, col = "steelblue", lwd = 2,
      ylim = ylim_all, xlab = "Dose", ylab = "Mean response",
      main = "True Dose-Response Curves")
curve(mean_fun_emax(x, 0, beta_emax), from = 0, to = 10, col = "darkorange", lwd = 2, add = TRUE)
curve(mean_fun_B(x, 0, beta_B), from = 0, to = 10, col = "forestgreen", lwd = 2, add = TRUE)
legend("bottomright", legend = c("linear", "emax", "B (sigmoid)"),
       col = c("steelblue", "darkorange", "forestgreen"), lwd = 2, bty = "n")

#-------------------------------------------------------------------------------
# Generate current data
# Independent of scenario (n0/prob.exch/nonexch.type/q) - only depends on
# dose.resp - so it's generated once and reused across all matching scenarios.
# Dose/trtgroup/weight/BASE are bootstrap-resampled (with replacement, by
# subject) from the real current-trial covariates each replicate.
#-------------------------------------------------------------------------------

set.seed(4001)

curr_idx     <- sample(nrow(current_real), n_curr * n_reps, replace = TRUE)
df_curr      <- current_real[curr_idx, ]

dose_mat     <- matrix(df_curr$DOSE_scaled, nrow = n_curr)
trtgroup_mat <- matrix(df_curr$trtgroup,    nrow = n_curr)
weight_mat   <- matrix(df_curr$WEIGHT,      nrow = n_curr)
base_mat     <- matrix(df_curr$BASE,        nrow = n_curr)  # real covariate, shared across truths

y_linear_mat <- mean_fun_linear(dose_mat, base_mat, beta_lin) +
  matrix(rnorm(n_curr * n_reps, 0, sd_truth), nrow = n_curr)
y_emax_mat <- mean_fun_emax(dose_mat, base_mat, beta_emax) +
  matrix(rnorm(n_curr * n_reps, 0, sd_truth), nrow = n_curr)

sim_data <- data.frame(
  iteration = rep(seq_len(n_reps), each = n_curr),
  trtpn     = as.integer(as.vector(dose_mat) > 0),
  trtgroup  = as.vector(trtgroup_mat),
  weight    = as.vector(weight_mat),
  dose      = as.vector(dose_mat),
  base      = as.vector(base_mat),
  y_linear  = as.vector(y_linear_mat),
  y_emax    = as.vector(y_emax_mat)
)

str(sim_data)
head(sim_data)

saveRDS(sim_data, file = file.path(save.dir, "simdata_current.rds"))

#-------------------------------------------------------------------------------
# Generate historical data
# One representative dataset per unique scenario (n0, prob.exch, nonexch.type,
# q, dose.resp) - not one per current-data iteration; the same historical
# dataset is reused across all of that scenario's current-data replicates.
# niter_search candidates are generated per scenario (vectorized), and the one
# whose residuals (y0 - true conditional mean) look most typical of iid
# N(0, sd_truth) noise is kept.
#-------------------------------------------------------------------------------

scenarios <- unique(grid[, c("q", "prob.exch", "nonexch.type", "dose.resp", "n0")])
rownames(scenarios) <- NULL

niter_search <- 10000

set.seed(4022)
hist_list <- vector("list", nrow(scenarios))

start_time <- Sys.time()
cat("Starting historical data generation at:", format(start_time), "\n\n")

for (i in seq_len(nrow(scenarios))) {

  sc          <- scenarios[i, ]
  n0          <- as.numeric(sc$n0)
  p           <- as.numeric(sc$prob.exch)
  q           <- as.numeric(sc$q)
  nonexch.type <- sc$nonexch.type
  dose.resp   <- as.character(sc$dose.resp)

  mean_fun <- get_mean_fun(dose.resp)
  beta     <- get_true_beta(dose.resp)

  n_exch   <- round(n0 * p)
  n_unexch <- n0 - n_exch

  idx0        <- sample(nrow(historical_real), n0 * niter_search, replace = TRUE)
  df_hist     <- historical_real[idx0, ]
  dose0_mat     <- matrix(df_hist$DOSE_scaled, nrow = n0)
  trtgroup0_mat <- matrix(df_hist$trtgroup,    nrow = n0)
  weight0_mat   <- matrix(df_hist$WEIGHT,      nrow = n0)
  base0_mat     <- matrix(df_hist$BASE,        nrow = n0)  # real covariate, shared across truths

  mu0_mat <- matrix(NA_real_, nrow = n0, ncol = niter_search)

  if (n_exch > 0) {
    idx_exch <- seq_len(n_exch)
    mu0_mat[idx_exch, ] <- mean_fun(dose0_mat[idx_exch, , drop = FALSE],
                                     base0_mat[idx_exch, , drop = FALSE], beta)
  }

  if (n_unexch > 0) {
    idx_unexch <- (n_exch + 1):n0
    if (nonexch.type == "A") {
      # Type A: mean shift, scaled by the residual SD (effect-size units)
      beta_shift     <- beta
      beta_shift$mu0 <- beta$mu0 + q * sd_truth
      mu0_mat[idx_unexch, ] <- mean_fun(dose0_mat[idx_unexch, , drop = FALSE],
                                         base0_mat[idx_unexch, , drop = FALSE], beta_shift)
    } else if (nonexch.type == "B") {
      # Type B: structural difference (steep sigmoid truth instead of the
      # current dose.resp's mean function)
      mu0_mat[idx_unexch, ] <- mean_fun_B(dose0_mat[idx_unexch, , drop = FALSE],
                                           base0_mat[idx_unexch, , drop = FALSE], beta_B)
    }
  }

  y0_mat <- mu0_mat + matrix(rnorm(n0 * niter_search, 0, sd_truth), nrow = n0)

  # --- Pick the candidate whose noise looks most "typical" ---
  resid_mat  <- y0_mat - mu0_mat
  mean_resid <- colMeans(resid_mat)
  sd_resid   <- sqrt((colMeans(resid_mat^2) - mean_resid^2) * n0 / (n0 - 1))
  se_mean    <- sd_truth / sqrt(n0)
  se_sd      <- sd_truth / sqrt(2 * (n0 - 1))
  score      <- (mean_resid / se_mean)^2 + ((sd_resid - sd_truth) / se_sd)^2
  best       <- which.min(score)

  hist_list[[i]] <- data.frame(
    trtpn0       = as.integer(dose0_mat[, best] > 0),
    trtgroup0    = trtgroup0_mat[, best],
    weight0      = weight0_mat[, best],
    dose0        = dose0_mat[, best],
    base0        = base0_mat[, best],
    y0           = y0_mat[, best],
    is_exch      = c(rep(1L, n_exch), rep(0L, n_unexch)),
    q            = q,
    prob.exch    = p,
    nonexch.type = nonexch.type,
    dose.resp    = dose.resp,
    n0           = n0
  )

  cat(sprintf("[%s] Selected best of %d: dose.resp=%-6s | n0=%2d | prob.exch=%.1f | q=%s | type=%s | score=%.2f\n",
              format(Sys.time(), "%H:%M:%S"), niter_search, dose.resp, n0, p,
              ifelse(is.na(q), "NA", q), ifelse(is.na(nonexch.type), "NA", nonexch.type), score[best]))
  flush.console()
}

cat("\nFinished! Total run time:",
    round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes.\n")

hist_data <- do.call(rbind, hist_list)
hist_data$dose.resp <- factor(hist_data$dose.resp, levels = c("linear", "emax"))

str(hist_data)
table(hist_data$dose.resp, hist_data$nonexch.type, useNA = "ifany")
table(hist_data$n0, hist_data$prob.exch)

saveRDS(hist_data, file = file.path(save.dir, "simdata_historical.rds"))
