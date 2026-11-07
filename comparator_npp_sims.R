# ------------------------------------------------------------------------------
# Program name: comparator_npp_sims.R
# Description.: Fits the normalized power prior (NPP) comparator for each
#               scenario/replicate, a0 ~ Beta(1, 1) on [0, 1], via hdbayes'
#               lm.npp (closed-form normalizing constant for the normal linear
#               model). Estimates the ATE at 40 and 60 IU/kg vs placebo, and
#               keeps the posterior a0. Run separately from comparator_sims.R
#               (same grid, data, and standardization); one SLURM array task
#               per grid_comparator.rds row.
# ------------------------------------------------------------------------------

# Load Libraries
library(hdbayes)
library(posterior)

## Set directory paths
sim.dir   <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3/Simulation_Study"
input.dir <- file.path(sim.dir, "Inputs")
save.dir  <- file.path(sim.dir, "Outputs", "Comparators_NPP")
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

# Load in simulation data
current    <- readRDS(file.path(input.dir, "simdata_current.rds"))
historical <- readRDS(file.path(input.dir, "simdata_historical.rds"))

# Load in grid of simulation parameters
grid <- readRDS(file.path(input.dir, "grid_comparator.rds"))

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

## Obtain file name based on id (same naming as comparator_sims.R)
filename <- file.path(save.dir, paste0('id_', id, '_n0_', n0.id, '_probexch_', prob.exch.id,
                                       '_nonexch_', nonexch.type.id, '_q_', q.id,
                                       '_doseresp_', dose.resp.id, '.rds'))

set.seed(seed.id)


################################################################################
# ------------------------------------------------------------------------------
# Setup code for NPP (mirrors comparator_sims.R)
# ------------------------------------------------------------------------------

# Historical dataset for this scenario (NA-safe match on q/nonexch.type)
match_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)

hist_i <- historical[
  as.character(historical$dose.resp) == dose.resp.id &
    historical$n0        == n0.id &
    historical$prob.exch == prob.exch.id &
    match_na(historical$q, q.id) &
    match_na(as.character(historical$nonexch.type), nonexch.type.id),
]

y_col <- if (dose.resp.id == "linear") "y_linear" else "y_emax"

hist_df <- data.frame(y = hist_i$y0, dose = hist_i$dose0, base = hist_i$base0)

# ATE doses (IU/kg); per-subject dose = dose.per.kg * weight / 1000
dose_per_kg_lo <- 40
dose_per_kg_hi <- 60

# Standardize y/dose/base with the current-trial mean/SD (same as
# comparator_sims.R); dose slope back-transformed via beta_std * y_sd / dose_sd.
dose_mean <- mean(current$dose);     dose_sd <- sd(current$dose)
base_mean <- mean(current$base);     base_sd <- sd(current$base)
y_mean    <- mean(current[[y_col]]); y_sd    <- sd(current[[y_col]])
std_y    <- function(v) (v - y_mean)    / y_sd
std_dose <- function(d) (d - dose_mean) / dose_sd
std_base <- function(b) (b - base_mean) / base_sd
std_df   <- function(df) data.frame(y = std_y(df$y), dose = std_dose(df$dose), base = std_base(df$base))
hist_df_std <- std_df(hist_df)

# Priors (standardized scale). lm.npp's initial prior is normal-inverse-gamma:
# beta | sigma^2 ~ N(0, 10^2 sigma^2), sigma^2 ~ IG(2.1, 1.1) (hdbayes default).
# a0 ~ Beta(1, 1) on [0, 1].
hd_beta_sd        <- 10
npp_sigmasq_shape <- 2.1
npp_sigmasq_scale <- 1.1
npp_a0_shape1     <- 1
npp_a0_shape2     <- 1

hd_chains         <- 1L
hd_iter_warmup    <- 1000L
hd_iter_sampling  <- 2000L
hd_refresh        <- 0L

################################################################################

# ------------------------------------------------------------------------------
# Begin simulation code

summarize_ate <- function(draws_mat, method, iteration) {
  s <- data.frame(summarize_draws(draws_mat, "mean", "median", "sd", "mad",
                ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean", "ess_bulk", "ess_tail"))
  s$method    <- method
  s$iteration <- iteration
  s
}

results.all  <- data.frame()
failed.iters <- data.frame()   # replicates skipped due to a fitting error

# Fit the NPP to one replicate
run_one_iteration <- function(dfn) {

  current_i <- subset(current, iteration == dfn)

  curr_df     <- data.frame(y = current_i[[y_col]], dose = current_i$dose, base = current_i$base)
  curr_df_std <- std_df(curr_df)

  # Mean per-subject dose at 40/60 IU/kg (1000 IU)
  dose_lo <- mean(dose_per_kg_lo * current_i$weight / 1000)
  dose_hi <- mean(dose_per_kg_hi * current_i$weight / 1000)

  ## Normalized power prior, y ~ dose + base
  fit_npp <- lm.npp(
    formula       = y ~ dose + base,
    data.list     = list(curr_df_std, hist_df_std),
    beta.sd       = hd_beta_sd,
    sigmasq.shape = npp_sigmasq_shape,
    sigmasq.scale = npp_sigmasq_scale,
    a0.shape1     = npp_a0_shape1,
    a0.shape2     = npp_a0_shape2,
    chains        = hd_chains,
    iter_warmup   = hd_iter_warmup,
    iter_sampling = hd_iter_sampling,
    refresh       = hd_refresh
  )
  beta_dose <- as.numeric(fit_npp[["dose"]]) * y_sd / dose_sd
  draws_npp <- cbind(ATE_4 = beta_dose * dose_lo, ATE_6 = beta_dose * dose_hi,
                     a0    = as.numeric(fit_npp[["a0_hist_1"]]))
  summarize_ate(draws_npp, "npp", dfn)
}

# Run every replicate; a fitting error skips only that replicate
start <- Sys.time()
for (dfn in iter.start.id:iter.end.id) {

  iter.out <- tryCatch(
    run_one_iteration(dfn),
    error = function(e) {
      cat("Iteration:", dfn, "FAILED -", conditionMessage(e), "\n")
      failed.iters <<- rbind(failed.iters,
                              data.frame(iteration = dfn, error = conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(iter.out)) results.all <- rbind(results.all, iter.out)

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
  , 'failed'  = failed.iters
)

saveRDS(lst, filename)
