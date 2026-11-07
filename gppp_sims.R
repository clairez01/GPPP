# ------------------------------------------------------------------------------
# Program name: gppp_sims.R
# Description.: Simulation study on data generated from the CSL830-based DGP
#               (generate_sim_data.R). Fits the GP power prior (GPPP; random
#               and fixed a0) for each scenario/replicate. Companion program
#               comparator_sims.R runs the hdbayes comparator methods
#               separately (GPPP is ~99% of per-dataset runtime, so it's
#               split into its own finer-grained job array). Designed to run
#               as one SLURM array task per scenario/sub-job, using
#               grid_gppp.rds.
# ------------------------------------------------------------------------------

# Load Libraries
#devtools::install("gpppSampler", quiet = TRUE)
library(gpppSampler)
library(posterior)

## Set directory paths
sim.dir   <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3/Simulation_Study"
input.dir <- file.path(sim.dir, "Inputs")
save.dir  <- file.path(sim.dir, "Outputs", "GPPP")
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

# Load in simulation data
current    <- readRDS(file.path(input.dir, "simdata_current.rds"))
historical <- readRDS(file.path(input.dir, "simdata_historical.rds"))

# Load in grid of simulation parameters
grid <- readRDS(file.path(input.dir, "grid_gppp.rds"))

## get task ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(id))
  id <- 1

## get this job's sim parameters
grid.id          <- grid[id, ]
seed.id          <- as.integer(grid.id$seed)
n0.id            <- as.numeric(grid.id$n0)
q.id             <- as.numeric(grid.id$q)
prob.exch.id     <- as.numeric(grid.id$prob.exch)
nonexch.type.id  <- as.character(grid.id$nonexch.type)
dose.resp.id     <- as.character(grid.id$dose.resp)
iter.start.id    <- as.integer(grid.id$iter.start)
iter.end.id      <- as.integer(grid.id$iter.end)

## Obtain file name based on id
filename <- file.path(save.dir, paste0('id_', id, '_n0_', n0.id, '_probexch_', prob.exch.id,
                                       '_nonexch_', nonexch.type.id, '_q_', q.id,
                                       '_doseresp_', dose.resp.id, '.rds'))

set.seed(seed.id)


################################################################################
# ------------------------------------------------------------------------------
# Setup code for GPPP model
# ------------------------------------------------------------------------------

# Match this job's scenario to its historical dataset (one per scenario; see
# generate_sim_data.R), NA-safe on q/nonexch.type
match_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)

hist_i <- historical[
  as.character(historical$dose.resp) == dose.resp.id &
    historical$n0        == n0.id &
    historical$prob.exch == prob.exch.id &
    match_na(historical$q, q.id) &
    match_na(as.character(historical$nonexch.type), nonexch.type.id),
]

y_col <- if (dose.resp.id == "linear") "y_linear" else "y_emax"

X0_raw <- cbind(dose = hist_i$dose0, base = hist_i$base0)
y0_raw <- hist_i$y0

# Standardized covariates/outcome. GPPP's eta2 prior (log-normal(0,1))
# assumes data on a comparable scale. Constants are computed once from the
# full current-trial sample (stable across all MC replicates for a given
# dose.resp) and reused for the historical data and every prediction grid
# below.
dose_mean <- mean(current$dose);     dose_sd <- sd(current$dose)
base_mean <- mean(current$base);     base_sd <- sd(current$base)
y_mean    <- mean(current[[y_col]]); y_sd    <- sd(current[[y_col]])
std_dose  <- function(d) (d - dose_mean) / dose_sd
std_base  <- function(b) (b - base_mean) / base_sd
std_y     <- function(v) (v - y_mean)    / y_sd

X0 <- cbind(dose = std_dose(X0_raw[, "dose"]), base = std_base(X0_raw[, "base"]))
y0 <- std_y(y0_raw)

# Prior scale for the noise variance, from the (standardized) historical outcome
sig2_prior <- 1

gppp_prior <- list(log_eta2_mean = 0, log_eta2_sd = 1,
                    log_ell_mean  = 0, log_ell_sd  = 1,
                    tau_alpha     = 2, tau_beta     = sig2_prior)
gppp_mhsd  <- list(log_eta2 = 0.4, log_ell = 0.4)

# n_iter = post-warmup draws retained; gppp_sampler()'s n_iter argument is
# the total (including warmup), so n_iter_total is what's actually passed in.
n_warmup <- 2000L
n_iter   <- 5000L
n_iter_total <- n_warmup + n_iter

# ATE contrast doses (IU/kg); per-subject dose = dose.per.kg * weight / 1000
dose_per_kg_lo <- 40
dose_per_kg_hi <- 60

gppp_a0_grid <- c(0, 0.5, 1)  # GPPP with fixed a0
# gppp_sampler() requires 0 < a0_fixed < 1, so the boundary values are
# approximated with this epsilon
gppp_a0_eps  <- 1e-3

################################################################################

# ------------------------------------------------------------------------------
# Begin simulation code

# Bayesian-bootstrap ATE: Dirichlet-weighted average of the per-subject
# predicted response difference, for each posterior draw
bb_ate <- function(Fd_hi, Fd_lo) {
  n_save <- nrow(Fd_hi)
  n_x    <- ncol(Fd_hi)
  vapply(seq_len(n_save), function(s) {
    w <- diff(c(0, sort(runif(n_x - 1L)), 1))
    sum(w * (Fd_hi[s, ] - Fd_lo[s, ]))
  }, numeric(1L))
}

summarize_ate <- function(draws_mat, method, iteration) {
  s <- data.frame(summarize_draws(draws_mat, "mean", "median", "sd", "mad",
                ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean", "ess_bulk", "ess_tail"))
  s$method    <- method
  s$iteration <- iteration
  s
}

results.all <- data.frame()

start <- Sys.time()
for (dfn in iter.start.id:iter.end.id) {

  current_i <- subset(current, iteration == dfn)

  X_raw <- cbind(dose = current_i$dose, base = current_i$base)
  y_raw <- current_i[[y_col]]

  X <- cbind(dose = std_dose(X_raw[, "dose"]), base = std_base(X_raw[, "base"]))
  y <- std_y(y_raw)

  # Prediction grids (40/60 IU/kg vs 0) at each subject's own weight and base;
  # ATE is rescaled to CHG via *y_sd after prediction
  base_pred   <- X[, "base"]
  weight_pred <- current_i$weight
  Xp_dlo <- cbind(dose = std_dose(dose_per_kg_lo * weight_pred / 1000), base = base_pred)
  Xp_dhi <- cbind(dose = std_dose(dose_per_kg_hi * weight_pred / 1000), base = base_pred)
  Xp_d0  <- cbind(dose = std_dose(0),                                  base = base_pred)

  ## ---- GPPP, random a0 ~ Uniform(0,1) ----
  fit <- gppp_sampler(
    y        = y,
    X        = X,
    y0       = y0,
    X0       = X0,
    n_iter   = n_iter_total,
    n_warmup = n_warmup,
    slice    = TRUE,
    kernel   = "matern52",
    a0_fixed = NULL,      # random a0 ~ Uniform(0,1)
    prior    = gppp_prior,
    mh_sd    = gppp_mhsd,
    verbose  = FALSE
  )

  Fdlo <- predict_gppp(fit, y, X, y0, X0, Xp = Xp_dlo)
  Fdhi <- predict_gppp(fit, y, X, y0, X0, Xp = Xp_dhi)
  Fd0  <- predict_gppp(fit, y, X, y0, X0, Xp = Xp_d0)

  draws <- cbind(a0 = fit$samples$a0,
                 ATE_4 = bb_ate(Fdlo, Fd0) * y_sd,
                 ATE_6 = bb_ate(Fdhi, Fd0) * y_sd)
  results.all <- rbind(results.all, summarize_ate(draws, "gppp", dfn))

  ## ---- GPPP, fixed a0 ----
  for (a0 in gppp_a0_grid) {
    a0_use <- if (a0 == 0) gppp_a0_eps else if (a0 == 1) 1 - gppp_a0_eps else a0

    fit_fixed <- gppp_sampler(
      y        = y,
      X        = X,
      y0       = y0,
      X0       = X0,
      n_iter   = n_iter_total,
      n_warmup = n_warmup,
      slice    = TRUE,
      kernel   = "matern52",
      a0_fixed = a0_use,
      prior    = gppp_prior,
      mh_sd    = gppp_mhsd,
      verbose  = FALSE
    )

    Fdlo_f <- predict_gppp(fit_fixed, y, X, y0, X0, Xp = Xp_dlo)
    Fdhi_f <- predict_gppp(fit_fixed, y, X, y0, X0, Xp = Xp_dhi)
    Fd0_f  <- predict_gppp(fit_fixed, y, X, y0, X0, Xp = Xp_d0)

    draws_fixed <- cbind(ATE_4 = bb_ate(Fdlo_f, Fd0_f) * y_sd,
                          ATE_6 = bb_ate(Fdhi_f, Fd0_f) * y_sd)
    results.all <- rbind(results.all, summarize_ate(draws_fixed, paste0("gppp_a0_", a0), dfn))
  }

  cat("Iteration:", dfn, "\n")
}

end <- Sys.time()
end - start

# End sim code
# ------------------------------------------

## SAVE THE RESULTS

lst <- list(
  'simscen' = grid.id
  , 'id'      = id
  , 'simres'  = results.all
)

saveRDS(lst, filename)
