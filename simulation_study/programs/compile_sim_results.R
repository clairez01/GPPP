# ------------------------------------------------------------------------------
# Program name: compile_sim_results.R
# Description.: Combine gppp_sims.R (GPPP/), comparator_sims.R
#               (Comparators/), and comparator_npp_sims.R (Comparators_NPP/)
#               output across all sub-jobs into one summary
#               table per scenario x method x variable, with bias/MSE/
#               relative MSE (vs ANOVA)/coverage/CI-width/interval score
#               against the oracle truth from the DGP.
# ------------------------------------------------------------------------------

library(haven)
library(dplyr)

## Set directory paths
paper.dir   <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3"
sim.dir     <- file.path(paper.dir, "Simulation_Study")
input.dir   <- file.path(sim.dir, "Inputs")
gppp.dir    <- file.path(sim.dir, "Outputs", "GPPP")
comp.dir    <- file.path(sim.dir, "Outputs", "Comparators")
npp.dir     <- file.path(sim.dir, "Outputs", "Comparators_NPP")
results.dir <- file.path(sim.dir, "Results")

# ------------------------------------------------------------------------------
# Oracle truth: true ATE_4 / ATE_6 (40 and 60 IU/kg vs placebo)
# = mean over real subjects of f(dose.per.kg * WEIGHT / 1000) - f(0)
# (exact, since the sims bootstrap these subjects uniformly). Used for every
# method except ancova/anova (see group-specific truth below).
# ------------------------------------------------------------------------------

csl.dir <- file.path(paper.dir, "Datasets")
current_real <- read_sas(file.path(csl.dir, "csl830_current.sas7bdat"))
current_real$DOSE_scaled <- current_real$DOSE / 1000

fit_lin_oracle    <- lm(CHG ~ DOSE_scaled + BASE, data = current_real)
TE_target_oracle  <- unname(coef(fit_lin_oracle)["DOSE_scaled"]) * 5
ED50_oracle       <- 5

# True effect vs placebo, f(dose) - f(0)
oracle_effect <- list(
  linear = function(dose) (TE_target_oracle / 5) * dose,
  emax   = function(dose) (2 * TE_target_oracle) * dose / (ED50_oracle + dose)
)

dose_per_kg_lo <- 40
dose_per_kg_hi <- 60

true_vals <- lapply(oracle_effect, function(eff) {
  c(ATE_4 = mean(eff(dose_per_kg_lo * current_real$WEIGHT / 1000)),
    ATE_6 = mean(eff(dose_per_kg_hi * current_real$WEIGHT / 1000)))
})
print(true_vals)
cat("Mean weight (kg):", mean(current_real$WEIGHT), "\n")

# ------------------------------------------------------------------------------
# Group-specific truth for ancova/anova, per replicate. These treat dose as a
# 3-level factor, so their 40/60 IU/kg contrast estimates the mean true
# effect among the subjects actually in that dose group,
#   mean_{i in group} f(dose_i) - f(0),
# not the population-averaged oracle above. (In the real trial the 40/60
# IU/kg groups are heavier than average, so their total dose -- and effect --
# is larger; the bootstrap carries that into every replicate.) Groups are
# recovered exactly as in comparator_sims.R. Named by iteration.
# ------------------------------------------------------------------------------

sim_current <- readRDS(file.path(input.dir, "simdata_current.rds"))
sim_current$dosegrp <- round(sim_current$dose * 1000 / sim_current$weight)

group_true_vals <- lapply(oracle_effect, function(eff) {
  per_iter <- split(sim_current, sim_current$iteration)
  list(
    ATE_4 = vapply(per_iter, function(d) mean(eff(d$dose[d$dosegrp == 40])), numeric(1)),
    ATE_6 = vapply(per_iter, function(d) mean(eff(d$dose[d$dosegrp == 60])), numeric(1))
  )
})
group_truth_methods <- c("ancova", "anova")
cat("Group-specific truth (mean over replicates):\n")
print(lapply(group_true_vals, function(g) sapply(g, mean, na.rm = TRUE)))

# ------------------------------------------------------------------------------
# Load grids and the 36 shared scenarios
# ------------------------------------------------------------------------------

grid_gppp <- readRDS(file.path(input.dir, "grid_gppp.rds"))
grid_comp <- readRDS(file.path(input.dir, "grid_comparator.rds"))

scen_cols <- c("q", "prob.exch", "nonexch.type", "dose.resp", "n0")
scenarios <- unique(grid_gppp[, scen_cols])
rownames(scenarios) <- NULL

match_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)

# Build the id_..._doseresp_....rds filename gppp_sims.R / comparator_sims.R /
# comparator_npp_sims.R use, given a grid row (its own row index is its SLURM_ARRAY_TASK_ID)
build_filename <- function(save.dir, id, grid.row) {
  file.path(save.dir, paste0(
    'id_', id, '_n0_', grid.row$n0, '_probexch_', grid.row$prob.exch,
    '_nonexch_', grid.row$nonexch.type, '_q_', grid.row$q,
    '_doseresp_', grid.row$dose.resp, '.rds'
  ))
}

# Read and row-bind every sub-job's $simres for one scenario from one grid
read_scenario_files <- function(grid, save.dir, scen) {
  match_idx <- which(
    match_na(grid$q, scen$q) &
      grid$prob.exch == scen$prob.exch &
      match_na(as.character(grid$nonexch.type), as.character(scen$nonexch.type)) &
      grid$dose.resp == scen$dose.resp &
      grid$n0 == scen$n0
  )
  dfs <- list()
  for (id in match_idx) {
    f <- build_filename(save.dir, id, grid[id, ])
    if (file.exists(f)) {
      dfs[[as.character(id)]] <- readRDS(f)$simres
    }
  }
  if (length(dfs) == 0) return(NULL)
  dplyr::bind_rows(dfs)
}

# Count replicates comparator_sims.R skipped due to a fitting error (its
# per-sub-job $failed log; see comparator_sims.R's run_one_iteration/
# tryCatch), summed across every sub-job matching this scenario. A failed
# replicate contributes zero rows to EVERY comparator method for that
# iteration, so this is what a lower-than-expected n_reps in the compiled
# output for those methods traces back to. GPPP files have no $failed
# element (gppp_sims.R has no equivalent try/catch), so this is always 0
# there -- the column reflects only the comparator pipeline.
# comparator_npp_sims.R keeps its own $failed log (only the NPP row is lost
# for those replicates); it's read the same way, counted in n_failed_npp.
read_scenario_failures <- function(grid, save.dir, scen) {
  match_idx <- which(
    match_na(grid$q, scen$q) &
      grid$prob.exch == scen$prob.exch &
      match_na(as.character(grid$nonexch.type), as.character(scen$nonexch.type)) &
      grid$dose.resp == scen$dose.resp &
      grid$n0 == scen$n0
  )
  failed <- list()
  for (id in match_idx) {
    f <- build_filename(save.dir, id, grid[id, ])
    if (file.exists(f)) {
      lst <- readRDS(f)
      if (!is.null(lst$failed) && nrow(lst$failed) > 0) {
        failed[[as.character(id)]] <- lst$failed
      }
    }
  }
  if (length(failed) == 0) return(data.frame(iteration = integer(0), error = character(0)))
  dplyr::bind_rows(failed)
}

# ------------------------------------------------------------------------------
# Bias / MSE / coverage / width for one method x variable, against a fixed
# true value, PLUS the Monte Carlo error of each of those summary metrics
# (i.e. how precisely bias/MSE/coverage/width are themselves estimated,
# given only n_reps replicates -- not the posterior CI width metric itself).
# Standard formulas (Morris, White & Crowther 2019, Stat Med):
#   MCSE(bias)     = sd(estimate_i) / sqrt(n_reps)
#   MCSE(MSE)      = sd((estimate_i - true)^2) / sqrt(n_reps)
#   MCSE(coverage) = sqrt(cov_prob * (1 - cov_prob) / n_reps)   [proportion SE]
#   MCSE(width)    = sd(width_i) / sqrt(n_reps)
#   MCSE(int. score) = sd(IS_i) / sqrt(n_reps)
# 95% bounds are estimate +/- 1.96*MCSE (coverage's clipped to [0, 1] since
# it's a proportion). Returns one row PER METRIC (long format: bias, MSE,
# cov_prob, width, interval_score), which is what a faceted point-range plot needs directly
# (facet by `metric`, y = estimate, ymin = LB95, ymax = UB95).
#
# df is the subset of rows for that method/variable across all
# replicates/sub-jobs. true_val is either a single value (the oracle) or a
# vector named by iteration (the per-replicate group-specific truth), in
# which case each replicate is scored against its own truth and the
# reported true_val is their mean.
# ------------------------------------------------------------------------------

summarize_metric <- function(df, true_val, method_label, variable_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  n <- nrow(df)

  if (length(true_val) > 1) {
    true_val <- unname(true_val[as.character(df$iteration)])
    stopifnot(!anyNA(true_val))
  }

  errs    <- df$mean - true_val
  sqerrs  <- errs^2
  covind  <- df$X2.5. <= true_val & df$X97.5. >= true_val
  widths  <- df$X97.5. - df$X2.5.

  bias     <- mean(errs, na.rm = TRUE)
  bias_se  <- sd(errs, na.rm = TRUE) / sqrt(n)

  mse      <- mean(sqerrs, na.rm = TRUE)
  mse_se   <- sd(sqerrs, na.rm = TRUE) / sqrt(n)

  covprob  <- mean(covind, na.rm = TRUE)
  cov_se   <- sqrt(covprob * (1 - covprob) / n)

  width    <- mean(widths, na.rm = TRUE)
  width_se <- sd(widths, na.rm = TRUE) / sqrt(n)

  # Interval score for the 95% CI (Gneiting & Raftery 2007, JASA, Sec 6.2):
  # width + (2/alpha) * distance by which the truth falls outside the CI.
  # Rewards narrow intervals but penalizes misses, so it combines
  # coverage and width into one number (lower = better).
  alpha    <- 0.05
  iscores  <- widths +
    (2 / alpha) * (df$X2.5. - true_val) * (true_val < df$X2.5.) +
    (2 / alpha) * (true_val - df$X97.5.) * (true_val > df$X97.5.)
  iscore    <- mean(iscores, na.rm = TRUE)
  iscore_se <- sd(iscores, na.rm = TRUE) / sqrt(n)

  metrics <- data.frame(
    metric   = c("bias", "MSE", "cov_prob", "width", "interval_score"),
    estimate = c(bias, mse, covprob, width, iscore),
    mcse     = c(bias_se, mse_se, cov_se, width_se, iscore_se)
  )
  metrics$LB95 <- metrics$estimate - 1.96 * metrics$mcse
  metrics$UB95 <- metrics$estimate + 1.96 * metrics$mcse
  is_cov <- metrics$metric == "cov_prob"
  metrics$LB95[is_cov] <- pmax(metrics$LB95[is_cov], 0)
  metrics$UB95[is_cov] <- pmin(metrics$UB95[is_cov], 1)

  metrics$method     <- method_label
  metrics$variable   <- variable_label
  metrics$n_reps     <- n
  metrics$true_val   <- mean(true_val)
  metrics$post_mean  <- mean(df$mean, na.rm = TRUE)
  metrics$post_sd    <- mean(df$sd,   na.rm = TRUE)
  metrics$post_LB    <- mean(df$X2.5.,  na.rm = TRUE)
  metrics$post_UB    <- mean(df$X97.5., na.rm = TRUE)

  metrics
}

# Relative MSE vs a reference method (ANOVA): MSE_method / MSE_ref, each
# scored against its own truth. Every method is fit to the same current-data
# replicate (matched by iteration), so the two MSEs are positively correlated;
# the MCSE uses the paired delta-method formula for a ratio of means
# (cf. Morris, White & Crowther 2019, Table 6, relative precision):
#   Var(R) ~= R^2/n * [Var(a)/mean(a)^2 + Var(b)/mean(b)^2
#                      - 2 Cov(a, b)/(mean(a) mean(b))]
# with a_i, b_i the per-replicate squared errors. Only iterations present for
# both methods are used. Same long-format columns as summarize_metric().
summarize_rel_mse <- function(df, true_val, ref_df, ref_true_val,
                              method_label, variable_label) {
  if (is.null(df) || nrow(df) == 0 || is.null(ref_df) || nrow(ref_df) == 0) return(NULL)

  sq_errs <- function(d, tv) {
    if (length(tv) > 1) tv <- unname(tv[as.character(d$iteration)])
    setNames((d$mean - tv)^2, d$iteration)
  }
  a <- sq_errs(df, true_val)
  b <- sq_errs(ref_df, ref_true_val)
  common <- intersect(names(a), names(b))
  a <- a[common]
  b <- b[common]
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]
  b <- b[ok]
  n <- length(a)

  ma  <- mean(a)
  mb  <- mean(b)
  rel <- ma / mb
  rel_var <- rel^2 / n * (var(a) / ma^2 + var(b) / mb^2 - 2 * cov(a, b) / (ma * mb))
  rel_se  <- sqrt(max(rel_var, 0))   # exactly 0 for the reference vs itself

  data.frame(
    metric = "rel_MSE", estimate = rel, mcse = rel_se,
    LB95 = rel - 1.96 * rel_se, UB95 = rel + 1.96 * rel_se,
    method = method_label, variable = variable_label, n_reps = n
  )
}

# Descriptive (no true value/bias) MC summary of a posterior quantity's mean
# across replicates -- e.g. GPPP's learned a0. Same long format/column names
# as summarize_metric() (with NA for the metrics that don't apply) so it
# rbinds cleanly.
summarize_descriptive <- function(df, method_label, variable_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  n <- nrow(df)
  est <- mean(df$mean, na.rm = TRUE)
  se  <- sd(df$mean, na.rm = TRUE) / sqrt(n)

  data.frame(
    metric = paste0(variable_label, "_mean"), estimate = est, mcse = se,
    LB95 = est - 1.96 * se, UB95 = est + 1.96 * se,
    method = method_label, variable = variable_label, n_reps = n, true_val = NA,
    post_mean = est, post_sd = mean(df$sd, na.rm = TRUE),
    post_LB = mean(df$X2.5., na.rm = TRUE), post_UB = mean(df$X97.5., na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Main loop: one scenario at a time, GPPP + comparator results merged
# ------------------------------------------------------------------------------

compiled.results  <- data.frame()
compiled.failures <- data.frame()   # every skipped comparator replicate, with scenario + error message

for (i in seq_len(nrow(scenarios))) {
  scen <- scenarios[i, ]
  label <- sprintf("dose.resp=%s n0=%s prob.exch=%s nonexch=%s q=%s",
                    scen$dose.resp, scen$n0, scen$prob.exch,
                    ifelse(is.na(scen$nonexch.type), "NA", scen$nonexch.type),
                    ifelse(is.na(scen$q), "NA", scen$q))

  gppp_res <- read_scenario_files(grid_gppp, gppp.dir, scen)
  comp_res <- read_scenario_files(grid_comp, comp.dir, scen)
  npp_res  <- read_scenario_files(grid_comp, npp.dir,  scen)

  if (is.null(gppp_res) && is.null(comp_res) && is.null(npp_res)) {
    cat("Scenario", label, ": no files found, skipping.\n")
    next
  }

  comp_failed <- read_scenario_failures(grid_comp, comp.dir, scen)
  n_failed    <- nrow(comp_failed)
  if (n_failed > 0) {
    cat("Scenario", label, ":", n_failed, "comparator replicate(s) skipped due to a fitting error.\n")
    comp_failed$pipeline     <- "comparator"
    comp_failed$q            <- scen$q
    comp_failed$prob.exch    <- scen$prob.exch
    comp_failed$nonexch.type <- scen$nonexch.type
    comp_failed$dose.resp    <- scen$dose.resp
    comp_failed$n0           <- scen$n0
    compiled.failures <- rbind(compiled.failures, comp_failed)
  }

  npp_failed   <- read_scenario_failures(grid_comp, npp.dir, scen)
  n_failed_npp <- nrow(npp_failed)
  if (n_failed_npp > 0) {
    cat("Scenario", label, ":", n_failed_npp, "NPP replicate(s) skipped due to a fitting error.\n")
    npp_failed$pipeline     <- "npp"
    npp_failed$q            <- scen$q
    npp_failed$prob.exch    <- scen$prob.exch
    npp_failed$nonexch.type <- scen$nonexch.type
    npp_failed$dose.resp    <- scen$dose.resp
    npp_failed$n0           <- scen$n0
    compiled.failures <- rbind(compiled.failures, npp_failed)
  }

  results.all <- dplyr::bind_rows(gppp_res, comp_res, npp_res)
  tv  <- true_vals[[as.character(scen$dose.resp)]]
  gtv <- group_true_vals[[as.character(scen$dose.resp)]]

  # ATE_4 / ATE_6 for every method present; ancova/anova are scored against
  # the per-replicate group-specific truth, everything else the oracle
  methods <- unique(results.all$method[results.all$variable %in% c("ATE_4", "ATE_6")])
  scen_summary <- list()
  for (m in methods) {
    for (v in c("ATE_4", "ATE_6")) {
      df_mv <- results.all[results.all$method == m & results.all$variable == v, ]
      truth <- if (m %in% group_truth_methods) gtv[[v]] else tv[[v]]
      scen_summary[[paste(m, v)]] <- summarize_metric(df_mv, truth, m, v)

      # Relative MSE vs ANOVA fit to the same replicates in this scenario
      ref_df <- results.all[results.all$method == "anova" & results.all$variable == v, ]
      scen_summary[[paste(m, v, "rel")]] <-
        summarize_rel_mse(df_mv, truth, ref_df, gtv[[v]], m, v)
    }
  }

  # a0 (GPPP random-a0's and NPP's learned borrowing weight): no true value /
  # bias, just a descriptive posterior summary across replicates
  for (m in c("gppp", "npp")) {
    a0_df <- results.all[results.all$method == m & results.all$variable == "a0", ]
    scen_summary[[paste(m, "a0")]] <- summarize_descriptive(a0_df, m, "a0")
  }

  sim.summary <- dplyr::bind_rows(scen_summary)
  if (nrow(sim.summary) == 0) {
    cat("Scenario", label, ": files found but no matching rows, skipping.\n")
    next
  }

  sim.summary$q             <- scen$q
  sim.summary$prob.exch     <- scen$prob.exch
  sim.summary$nonexch.type  <- scen$nonexch.type
  sim.summary$dose.resp     <- scen$dose.resp
  sim.summary$n0            <- scen$n0
  sim.summary$n_failed      <- n_failed
  sim.summary$n_failed_npp  <- n_failed_npp

  compiled.results <- rbind(compiled.results, sim.summary)
  cat("Scenario", label, ": processed (", nrow(results.all), "rows across",
      length(unique(results.all$iteration)), "replicates ).\n")
}

head(compiled.results)

if (nrow(compiled.failures) > 0) {
  cat("\nTotal replicates skipped across all scenarios, by pipeline:\n")
  print(table(compiled.failures$pipeline))
}

saveRDS(compiled.results,  file = file.path(results.dir, "compiled_sim_results.rds"))
saveRDS(compiled.failures, file = file.path(results.dir, "compiled_sim_failures.rds"))
