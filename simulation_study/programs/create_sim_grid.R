# ------------------------------------------------------------------------------
# Program name: create_sim_grid.R
# Description.: Build the grid of simulation scenarios for the GPPP method.
#               Produces two expanded grids sharing the same 36 scenarios but
#               split into different numbers of cluster sub-jobs: GPPP is
#               ~99% of per-dataset runtime and needs a finer split to fit
#               within cluster job time limits; the hdbayes comparators are
#               cheap enough to need only a handful of sub-jobs.
# ------------------------------------------------------------------------------

## Directory to store files
save.dir <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3/Simulation_Study/Inputs"

# --- Simulation Parameters ---
n         <- 80                       # Current trial sample size (fixed)
n0_list   <- c(10, 20, 40)            # Historical sample sizes
q_bias    <- c(-1, -0.5, 0.5, 1)      # Mean shift multipliers for Scenario A (nonexchangeable)
                                       # Applied downstream as shift = q * sd(outcome), i.e. an
                                       # effect-size (SD) scale rather than a raw-unit shift
prob_exch <- 0                        # Fully nonexchangeable level (prob.exch = 0)
nonexch_types <- c("A", "B")          # A = mean shift, B = structural difference
dose_resp <- c("linear", "emax")      # Functional truths for the dose-response curve

n_reps_total <- 10000                 # Total Monte Carlo repetitions per scenario

n_subjobs_gppp       <- 40            # GPPP: 250 reps/sub-job
n_subjobs_comparator <- 4             # hdbayes comparators: 2,500 reps/sub-job

# --- Create Unified Grid ---

# Grid 1: Fully Exchangeable (prob.exch = 1)
# Nonexch.type/shift are irrelevant here since all historical data follow the
# same mean function as the current data, so these rows are not duplicated
# across nonexch.type or q.
grid1 <- expand.grid(
  q             = NA,
  prob.exch     = 1,
  nonexch.type  = NA,
  dose.resp     = dose_resp,
  n0            = n0_list
)

# Grid 2a: Fully Nonexchangeable (prob.exch = 0), Type A (mean shift)
# Crossed with all shift magnitudes q
grid2a <- expand.grid(
  q             = q_bias,
  prob.exch     = prob_exch,
  nonexch.type  = "A",
  dose.resp     = dose_resp,
  n0            = n0_list
)

# Grid 2b: Fully Nonexchangeable (prob.exch = 0), Type B (structural difference)
# No shift magnitude applies since nonexchangeability comes from a different
# functional form rather than a mean shift
grid2b <- expand.grid(
  q             = NA,
  prob.exch     = prob_exch,
  nonexch.type  = "B",
  dose.resp     = dose_resp,
  n0            = n0_list
)

# Combine and Sort
grid <- rbind(grid1, grid2a, grid2b)
# Sorting by dose-response type, then prob.exch, then nonexch.type, then bias, then sample size
grid <- grid[order(grid$dose.resp, -grid$prob.exch, grid$nonexch.type, grid$q, grid$n0), ]
rownames(grid) <- NULL

# Add fixed current sample size
grid$n <- n

# --- Expansion for Parallel Processing ---
# Split each scenario's n_reps_total repetitions across n_subjobs cluster
# jobs, so each sub-job processes n_reps_total / n_subjobs repetitions, and
# assign each sub-job an independent seed.
expand_grid <- function(grid, n_subjobs, seed) {
  n_iter_per_job <- n_reps_total / n_subjobs
  grid_expanded <- grid[rep(seq_len(nrow(grid)), each = n_subjobs), ]
  grid_expanded$subjob     <- rep(seq_len(n_subjobs), times = nrow(grid))
  grid_expanded$iter.start <- (grid_expanded$subjob - 1) * n_iter_per_job + 1
  grid_expanded$iter.end   <- grid_expanded$subjob * n_iter_per_job
  rownames(grid_expanded) <- NULL

  set.seed(seed)
  grid_expanded$seed <- sample(seq_len(1000 * nrow(grid_expanded)), nrow(grid_expanded), replace = FALSE)
  grid_expanded
}

grid_gppp       <- expand_grid(grid, n_subjobs_gppp, seed = 4030)
grid_comparator <- expand_grid(grid, n_subjobs_comparator, seed = 4031)

saveRDS(grid_gppp,       file = file.path(save.dir, "grid_gppp.rds"))
saveRDS(grid_comparator, file = file.path(save.dir, "grid_comparator.rds"))
