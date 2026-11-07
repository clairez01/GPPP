#' @export
print.gppp_fit <- function(x, ...) {
  ns           <- length(x$samples$eta2)
  d            <- x$d
  sigma2_fixed <- isTRUE(x$sigma2_fixed)
  outcome_type <- if (!is.null(x$outcome_type)) x$outcome_type else "continuous"

  cat("GP Power Prior MCMC Fit\n")
  cat("=======================\n")
  cat(sprintf("Outcome    : %s\n", outcome_type))
  cat(sprintf("Kernel     : %s%s\n", x$kernel,
              if (d > 1L) sprintf(" (ARD, d = %d)", d) else ""))
  cat(sprintf("Iterations : %d total | %d warmup | %d saved (thin = %d)\n",
              x$n_iter, x$n_warmup, ns, x$thin))
  cat(sprintf("a0         : %s\n",
              if (x$random_a0) "random (Uniform prior)" else
              sprintf("fixed = %.4f", x$a0_fixed)))
  cat("\nPosterior means:\n")
  cat(sprintf("  %-12s %.4f\n", "eta2", mean(x$samples$eta2)))
  if (d == 1L) {
    cat(sprintf("  %-12s %.4f\n", "ell", mean(x$samples$ell)))
  } else {
    for (j in seq_len(d))
      cat(sprintf("  %-12s %.4f\n", paste0("ell_", j), mean(x$samples$ell[, j])))
  }
  if (sigma2_fixed) {
    cat(sprintf("  %-12s fixed = 1.0000\n", "tau"))
  } else {
    cat(sprintf("  %-12s %.4f\n", "tau", mean(x$samples$tau)))
  }
  if (x$random_a0)
    cat(sprintf("  %-12s %.4f\n", "a0", mean(x$samples$a0)))
  cat(sprintf("  %-12s %.4f\n", "mu", mean(x$samples$mu)))

  cat("\nMH acceptance rates:\n")
  for (nm in names(x$acceptance_rates)) {
    rate <- x$acceptance_rates[nm]
    label <- if (is.na(rate)) {
      if (nm == "tau") {
        if (sigma2_fixed) "N/A (fixed)" else "N/A (conjugate gamma)"
      } else {
        "N/A (slice)"
      }
    } else {
      sprintf("%.1f%%", 100 * rate)
    }
    cat(sprintf("  %-8s %s\n", nm, label))
  }
  invisible(x)
}

#' @export
summary.gppp_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  d            <- object$d
  sigma2_fixed <- isTRUE(object$sigma2_fixed)

  summarise_vec <- function(s) {
    c(mean = mean(s), sd = sd(s), quantile(s, probs))
  }

  rows <- list(eta2 = summarise_vec(object$samples$eta2))

  if (d == 1L) {
    rows[["ell"]] <- summarise_vec(as.vector(object$samples$ell))
  } else {
    for (j in seq_len(d))
      rows[[paste0("ell_", j)]] <- summarise_vec(object$samples$ell[, j])
  }

  if (!sigma2_fixed)
    rows[["tau"]] <- summarise_vec(object$samples$tau)

  if (object$random_a0)
    rows[["a0"]] <- summarise_vec(object$samples$a0)

  rows[["mu"]] <- summarise_vec(object$samples$mu)

  res <- do.call(rbind, rows)
  cat("GP Power Prior Posterior Summary\n")
  cat("=================================\n")
  if (sigma2_fixed)
    cat("(binary outcomes: logit link, Polya-Gamma-augmented power prior; tau fixed at 1)\n")
  print(round(res, 4))
  invisible(res)
}
