# =============================================================================
# Final model-quality summary report
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (FINAL SUMMARY block).

print_final_summary <- function(all_results, corr_cols, bias) {

  selection_bias_detected <- bias$selection_bias_detected
  selection_bias_info     <- bias$selection_bias_info

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  FINAL MODEL QUALITY SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  for (outcome in corr_cols) {
    result <- all_results[[outcome]]

    cat(sprintf("Outcome: %s\n", outcome))
    cat(strrep("-", 60), "\n")
    cat(sprintf("  Random Effects Structure: %s\n", result$random_effects_formula))
    if (result$used_random_slopes) {
      cat("    → Random slopes for logT were justified and included\n")
    } else {
      cat("    → Random intercept only (slopes not needed)\n")
    }
    cat(sprintf("  ICC: %.3f (%.1f%% between-network variance)\n", result$icc, result$icc * 100))

    if (result$convergence_ok) {
      cat("  Convergence: ✓ All models converged successfully\n")
    } else {
      cat("  Convergence: ⚠ Warnings detected (see above)\n")
    }

    if (result$singular_fit) {
      cat("  Singularity: ⚠ Primary model has singular fit\n")
    } else {
      cat("  Singularity: ✓ No singularity issues\n")
    }

    delta_aic_display <- AIC(result$two_way) - AIC(result$three_way)
    cat(sprintf("  Primary Model: %s (ΔAIC 2-way vs 3-way = %.1f)\n",
                result$primary_label, delta_aic_display))
    cat("\n")
  }

  cat("═══════════════════════════════════════════════════════════════\n\n")

  # Report selection bias status in summary
  if (selection_bias_detected) {
    cat("⚠⚠⚠ CRITICAL LIMITATION ⚠⚠⚠\n")
    cat("Selection bias was detected in the data.\n")
    cat(sprintf("Biased predictors: %s\n",
                paste(selection_bias_info$biased_predictors, collapse = ", ")))
    cat("See detailed warning above for implications and recommendations.\n")
    cat("Results should be interpreted with caution.\n\n")
  }

  cat("Analysis completed!\n")
  cat("Review console output above for:\n")
  cat("  1. Selection bias warnings (CRITICAL if present)\n")
  cat("  2. Model convergence diagnostics\n")
  cat("  3. Random effects structure decisions\n")
  cat("  4. Singularity warnings\n\n")

  if (selection_bias_detected) {
    cat("IMPORTANT: Due to selection bias, consider:\n")
    cat("  • Reporting this limitation in your paper\n")
    cat("  • Conducting sensitivity analysis excluding high-failure conditions\n")
    cat("  • Collecting more data to improve estimation success rates\n\n")
  }
}
