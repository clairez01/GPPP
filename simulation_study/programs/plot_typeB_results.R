# ------------------------------------------------------------------------------
# Program name: plot_typeB_results.R
# Description.: Point-range plots comparing each method's performance under
#               Type B (structural difference) historical nonexchangeability
#               vs the fully-exchangeable (prob.exch = 1) reference. For one
#               (n0, dose.resp, ATE contrast) combination at a time: x = method,
#               two point-ranges per method (exchangeable vs Type B), faceted
#               by metric (|bias|, MSE relative to ANOVA, coverage, width,
#               interval score). Points = summary
#               estimate, whiskers = 95% Monte Carlo error (from
#               compile_sim_results.R's mcse/LB95/UB95 columns).
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

results.dir <- "C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/Paper 3/Simulation_Study/Results"
fig.dir     <- file.path(results.dir, "Figures/TypeB")
dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)

compiled <- readRDS(file.path(results.dir, "compiled_sim_results.rds"))

# Method display order/labels, matching what gppp_sims.R/comparator_sims.R/
# comparator_npp_sims.R produce
method_levels <- c("gppp", "gppp_a0_0", "gppp_a0_0.5", "gppp_a0_1",
                    "pp_a0_0", "pp_a0_0.5", "pp_a0_1", "npp",
                    "rmap", "ancova", "anova")
method_labels <- c("GPPP", "GP", "GPPP (a0=0.5)", "GPPP (a0=1)",
                    "BLR", "PP (a0=0.5)", "PP (a0=1)", "NPP",
                    "RMAP", "ANCOVA", "ANOVA")

# ATE contrasts as point shapes: circle = 40 IU/kg, triangle = 60 IU/kg
dose_levels <- c("ATE_4", "ATE_6")
dose_labels <- c("40 IU/kg vs. placebo", "60 IU/kg vs. placebo")
dose_shapes <- c(16, 17)

# Okabe-Ito (colorblind-safe): black for the exchangeable reference
exch_cols <- c("Fully Exchangeable" = "#000000",
               "Non-exchangeable (Functional Difference)" = "#CC79A7")

metric_levels <- c("bias", "rel_MSE", "cov_prob", "width", "interval_score")
metric_labels <- c("Absolute Bias", "Relative MSE", "95% CI Coverage", "95% CI Width",
                   "95% Interval Score")

# Reference lines: 0 for |bias|, 1 (= ANOVA) for relative MSE, 0.95 nominal coverage
ref_lines <- data.frame(
  metric     = factor(c("Absolute Bias", "Relative MSE", "95% CI Coverage"), levels = metric_labels),
  yintercept = c(0, 1, 0.95)
)

plot_caption <- paste0(
  "Relative MSE = MSE of each method / MSE of ANOVA (current data only, no baseline adjustment) ",
  "on the same replicates; values < 1 indicate a gain in efficiency over ANOVA.\n",
  "Interval score (Gneiting & Raftery, 2007) = 95% CI width + (2/0.05) x distance by which the truth ",
  "falls outside the CI; lower is better.\n",
  "Points = Monte Carlo estimates; whiskers = 95% Monte Carlo error intervals."
)

# Methods that never use the historical data (a0=0 for GPPP/PP, and the
# current-data-only comparators): their result is identical whether the
# (unused) historical data happens to be Type B or exchangeable -- showing
# both would just be the same current-data analysis re-run under two
# scenario labels, not a genuine comparison. Keep only the exchangeable
# run's version for these; every other method keeps both.
no_borrow_methods <- c("gppp_a0_0", "pp_a0_0", "ancova", "anova")

# Type B rows + the fully-exchangeable reference (nonexch.type is NA there)
typeB_data <- compiled %>%
  filter(variable %in% c("ATE_4", "ATE_6"),
         metric %in% metric_levels,
         method %in% method_levels,
         (nonexch.type == "B" & !is.na(nonexch.type)) | prob.exch == 1) %>%
  filter(!(method %in% no_borrow_methods & prob.exch == 0)) %>%
  mutate(
    exch_group = factor(ifelse(prob.exch == 1, "Fully Exchangeable", "Non-exchangeable (Functional Difference)"),
                         levels = c("Fully Exchangeable", "Non-exchangeable (Functional Difference)")),
    method = factor(method, levels = method_levels, labels = method_labels),
    metric = factor(metric, levels = metric_levels, labels = metric_labels),
    dose   = factor(variable, levels = dose_levels, labels = dose_labels)
  )

# Present Bias as |bias|. Folding a normal-approximation CI through abs():
# if the original 95% interval excludes 0, |bias|'s interval is
# [min(|LB|,|UB|), max(|LB|,|UB|)]; if it straddles 0 (bias could plausibly
# be exactly 0), the folded lower bound is 0 and the upper bound is
# max(|LB|,|UB|). The signed bias/CI stay untouched in compiled_sim_results
# itself -- this transform is presentation-only, applied here.
is_bias <- typeB_data$metric == "Absolute Bias"
straddles_zero <- is_bias & typeB_data$LB95 < 0 & typeB_data$UB95 > 0
folded_LB <- ifelse(straddles_zero, 0, pmin(abs(typeB_data$LB95), abs(typeB_data$UB95)))
folded_UB <- pmax(abs(typeB_data$LB95), abs(typeB_data$UB95))
typeB_data$LB95[is_bias]     <- folded_LB[is_bias]
typeB_data$UB95[is_bias]     <- folded_UB[is_bias]
typeB_data$estimate[is_bias] <- abs(typeB_data$estimate[is_bias])

# Posterior mean of a0 (GPPP random-a0 and NPP), averaged over replicates, for
# the footnote. a0 doesn't depend on the ATE contrast, only on the scenario.
a0_summary <- compiled %>%
  filter(variable == "a0", method %in% c("gppp", "npp"),
         (nonexch.type == "B" & !is.na(nonexch.type)) | prob.exch == 1)

a0_caption <- function(n0_val, dose_resp_val) {
  a <- a0_summary %>% filter(n0 == n0_val, dose.resp == dose_resp_val)
  if (nrow(a) == 0) return(NULL)
  get_a0 <- function(m, exch) {
    sprintf("%.2f", a$estimate[a$method == m & (a$prob.exch == 1) == exch])
  }
  paste0("Posterior mean of a0 (averaged over replicates): GPPP ", get_a0("gppp", TRUE),
         " (fully exchangeable), ", get_a0("gppp", FALSE), " (functional difference); NPP ",
         get_a0("npp", TRUE), ", ", get_a0("npp", FALSE), ". ",
         sprintf("MCSE ≤ %.3f.", ceiling(max(a$mcse) * 1000) / 1000))
}

# One figure per (n0, dose.resp): x = method, point-range per exch_group
# (colour) x ATE contrast (shape), dodged, faceted by metric.
plot_typeB <- function(df, n0_val, dose_resp_val, out_file = NULL) {
  d <- df %>% filter(n0 == n0_val, dose.resp == dose_resp_val)
  if (nrow(d) == 0) {
    warning("No data for n0=", n0_val, ", dose.resp=", dose_resp_val)
    return(NULL)
  }

  subtitle_expr <- bquote(
    "For scenario " * n[0] * " = " * .(n0_val) * ", dose response = " * .(dose_resp_val) * "."
  )

  p <- ggplot(d, aes(x = method, y = estimate, ymin = LB95, ymax = UB95,
                     color = exch_group, shape = dose,
                     group = interaction(dose, exch_group))) +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.5, fatten = 4,
                    linewidth = 0.4) +
    geom_hline(data = ref_lines, aes(yintercept = yintercept),
               linetype = "dashed", color = "grey50") +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_color_manual(values = exch_cols) +
    scale_shape_manual(values = dose_shapes) +
    guides(color = guide_legend(order = 1, nrow = 2), shape = guide_legend(order = 2, nrow = 2)) +
    labs(x = NULL, y = NULL, color = "Historical Data", shape = "ATE",
         title = "Operating Characteristics",
         subtitle = subtitle_expr,
         caption = paste(c(a0_caption(n0_val, dose_resp_val), plot_caption), collapse = "\n")) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          plot.caption = element_text(hjust = 0, size = 10))

  if (!is.null(out_file)) ggsave(file.path(fig.dir, out_file), p, width = 12, height = 10, dpi = 150)
  p
}

# Example: n0 = 10, dose.resp = linear (display only; not saved -- the loop
# below writes every combination)
p_example <- plot_typeB(typeB_data, n0_val = 10, dose_resp_val = "linear")
print(p_example)

# To generate every (n0 x dose.resp) combination:
 for (n0_val in c(10, 20, 40)) {
   for (dr in c("linear", "emax")) {
     plot_typeB(typeB_data, n0_val, dr,
                out_file = sprintf("typeB_n0-%d_%s.png", n0_val, dr))
   }
 }
