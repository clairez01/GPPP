# ------------------------------------------------------------------------------
# Program name: comparator_sims.R
# Description.: Fits the comparator methods for each scenario/replicate:
#               power prior (fixed a0 = 0, 0.5, 1; a0 = 0 = no borrowing),
#               robust MAP, and frequentist ANCOVA/ANOVA. Estimates the ATE
#               at 40 and 60 IU/kg vs placebo. GPPP runs separately in
#               gppp_sims.R. One SLURM array task per grid_comparator.rds row.
# ------------------------------------------------------------------------------

# Load Libraries
library(hdbayes)
library(posterior)

## Set directory paths
sim.dir   <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3/Simulation_Study"
input.dir <- file.path(sim.dir, "Inputs")
save.dir  <- file.path(sim.dir, "Outputs", "Comparators")
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

## Obtain file name based on id
filename <- file.path(save.dir, paste0('id_', id, '_n0_', n0.id, '_probexch_', prob.exch.id,
                                       '_nonexch_', nonexch.type.id, '_q_', q.id,
                                       '_doseresp_', dose.resp.id, '.rds'))

set.seed(seed.id)


################################################################################
# ------------------------------------------------------------------------------
# Setup code for comparator methods
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

X0_raw <- cbind(dose = hist_i$dose0, base = hist_i$base0)
y0_raw <- hist_i$y0

hist_df <- data.frame(y = y0_raw, dose = X0_raw[, "dose"], base = X0_raw[, "base"])

# ATE doses (IU/kg); per-subject dose = dose.per.kg * weight / 1000
dose_per_kg_lo <- 40
dose_per_kg_hi <- 60

pp_a0_grid <- c(0, 0.5, 1)    # power prior a0 values

# Standardize y/dose/base with the current-trial mean/SD (same as gppp_sims.R).
# Used for the Bayesian methods; the dose slope is back-transformed to the CHG
# scale via beta_raw = beta_std * y_sd / dose_sd.
dose_mean <- mean(current$dose);     dose_sd <- sd(current$dose)
base_mean <- mean(current$base);     base_sd <- sd(current$base)
y_mean    <- mean(current[[y_col]]); y_sd    <- sd(current[[y_col]])
std_y    <- function(v) (v - y_mean)    / y_sd
std_dose <- function(d) (d - dose_mean) / dose_sd
std_base <- function(b) (b - base_mean) / base_sd
std_df   <- function(df) data.frame(y = std_y(df$y), dose = std_dose(df$dose), base = std_base(df$base))
hist_df_std <- std_df(hist_df)

# Priors (standardized scale)
hd_beta_sd <- 10     # N(0, 10^2) on each coefficient
hd_disp_sd <- 10     # half-normal(0, 10^2) on sigma^2 (hdbayes' "dispersion")

# Robust MAP. w = prior weight on the MAP component (0.2 on the vague one).
# Vague component N(0, 1): too diffuse and it never wins the mixture weight.
# Between-trial SD prior: half-normal, scale = half each coefficient's
# unit-information SD, sigma / sd(x_j).
rmap_w          <- 0.8
rmap_vague_sd   <- 1
sig0_hist       <- sigma(lm(y ~ dose + base, data = hist_df_std))
rmap_meta_sd_sd <- 0.5 * sig0_hist / c(1, sd(hist_df_std$dose), sd(hist_df_std$base))  # (Intercept), dose, base

hd_chains        <- 1L
hd_iter_warmup   <- 1000L
hd_iter_sampling <- 2000L
hd_refresh       <- 0L

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

# lm() version of summarize_ate(): estimate, SE (in sd), t-based 95% CI
summarize_lm <- function(fit, method, iteration) {
  est <- coef(fit)[c("dosegrp40", "dosegrp60")]
  se  <- sqrt(diag(vcov(fit)))[c("dosegrp40", "dosegrp60")]
  ci  <- confint(fit, parm = c("dosegrp40", "dosegrp60"), level = 0.95)
  data.frame(variable = c("ATE_4", "ATE_6"), mean = unname(est), median = unname(est),
             sd = unname(se), mad = NA_real_, X2.5. = unname(ci[, 1]), X97.5. = unname(ci[, 2]),
             rhat = NA_real_, mcse_mean = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
             method = method, iteration = iteration)
}

results.all  <- data.frame()
failed.iters <- data.frame()   # replicates skipped due to a fitting error

# Fit all methods to one replicate
run_one_iteration <- function(dfn) {

  current_i <- subset(current, iteration == dfn)

  X_raw <- cbind(dose = current_i$dose, base = current_i$base)
  y_raw <- current_i[[y_col]]

  curr_df     <- data.frame(y = y_raw, dose = X_raw[, "dose"], base = X_raw[, "base"])
  curr_df_std <- std_df(curr_df)

  # Mean per-subject dose at 40/60 IU/kg (1000 IU)
  dose_lo <- mean(dose_per_kg_lo * current_i$weight / 1000)
  dose_hi <- mean(dose_per_kg_hi * current_i$weight / 1000)

  iter.results <- data.frame()

  ## Power prior, y ~ dose + base
  for (a0 in pp_a0_grid) {
    fit_pp <- glm.pp(
      formula       = y ~ dose + base,
      family        = gaussian(),
      data.list     = list(curr_df_std, hist_df_std),
      a0.vals       = a0,
      beta.sd       = hd_beta_sd,
      disp.sd       = hd_disp_sd,
      chains        = hd_chains,
      iter_warmup   = hd_iter_warmup,
      iter_sampling = hd_iter_sampling,
      refresh       = hd_refresh
    )
    beta_dose <- as.numeric(fit_pp[["dose"]]) * y_sd / dose_sd
    draws_pp  <- cbind(ATE_4 = beta_dose * dose_lo, ATE_6 = beta_dose * dose_hi)
    iter.results <- rbind(iter.results, summarize_ate(draws_pp, paste0("pp_a0_", a0), dfn))
  }

  ## Robust MAP, y ~ dose + base
  fit_rmap <- glm.rmap(
    formula         = y ~ dose + base,
    family          = gaussian(),
    data.list       = list(curr_df_std, hist_df_std),
    w               = rmap_w,
    meta.mean.sd    = hd_beta_sd,
    meta.sd.sd      = rmap_meta_sd_sd,
    norm.vague.mean = 0,
    norm.vague.sd   = rmap_vague_sd,
    disp.sd         = hd_disp_sd,
    chains          = hd_chains,
    iter_warmup     = hd_iter_warmup,
    iter_sampling   = hd_iter_sampling,
    refresh         = hd_refresh
  )
  beta_dose_rmap <- as.numeric(fit_rmap$post.samples[["dose"]]) * y_sd / dose_sd
  draws_rmap <- cbind(ATE_4 = beta_dose_rmap * dose_lo, ATE_6 = beta_dose_rmap * dose_hi)
  iter.results <- rbind(iter.results, summarize_ate(draws_rmap, "rmap", dfn))

  ## ANCOVA / ANOVA: OLS on the raw scale, current data only, dose group as
  ## a factor (placebo, 40, 60 IU/kg)
  dosegrp  <- factor(round(X_raw[, "dose"] * 1000 / current_i$weight), levels = c(0, 40, 60))

  ## ANCOVA: dose group + baseline covariate
  fit_ancova <- lm(y ~ dosegrp + base, data = data.frame(y = y_raw, dosegrp = dosegrp, base = X_raw[, "base"]))
  iter.results <- rbind(iter.results, summarize_lm(fit_ancova, "ancova", dfn))

  ## ANOVA: dose group only
  fit_anova <- lm(y ~ dosegrp, data = data.frame(y = y_raw, dosegrp = dosegrp))
  iter.results <- rbind(iter.results, summarize_lm(fit_anova, "anova", dfn))

  iter.results
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
