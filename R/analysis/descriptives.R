# =============================================================================
# Descriptive statistics (z-aggregated, then back-transformed)
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PART 3.1).

# Returns list(table = descriptive_stats_table) and writes HTML/LaTeX tables.
descriptive_stats <- function(dat_sim, corr_cols, available_other_metrics) {

  descriptive_stats_table <- data.frame(
    Metric = character(),
    Mean = numeric(),
    SD = numeric(),
    Median = numeric(),
    Min = numeric(),
    Max = numeric(),
    Range = numeric(),
    N = integer(),
    stringsAsFactors = FALSE
  )

  # 1. CORRELATIONS: Back-transform from z-scale
  for (outcome in corr_cols) {
    outcome_z <- paste0(outcome, "_z")

    if (outcome_z %in% names(dat_sim)) {
      # Back-transform z-values to r-scale
      r_values <- tanh(dat_sim[[outcome_z]])

      descriptive_stats_table <- rbind(
        descriptive_stats_table,
        data.frame(
          Metric = outcome,
          Mean = round(mean(r_values, na.rm = TRUE), 2),
          SD = round(sd(r_values, na.rm = TRUE), 2),
          Median = round(median(r_values, na.rm = TRUE), 2),
          Min = round(min(r_values, na.rm = TRUE), 2),
          Max = round(max(r_values, na.rm = TRUE), 2),
          Range = round(max(r_values, na.rm = TRUE) - min(r_values, na.rm = TRUE), 2),
          N = sum(!is.na(r_values)),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # 2. OTHER METRICS: Already aggregated on original scale
  for (metric in available_other_metrics) {
    if (metric %in% names(dat_sim)) {
      values <- dat_sim[[metric]]

      descriptive_stats_table <- rbind(
        descriptive_stats_table,
        data.frame(
          Metric = metric,
          Mean = round(mean(values, na.rm = TRUE), 2),
          SD = round(sd(values, na.rm = TRUE), 2),
          Median = round(median(values, na.rm = TRUE), 2),
          Min = round(min(values, na.rm = TRUE), 2),
          Max = round(max(values, na.rm = TRUE), 2),
          Range = round(max(values, na.rm = TRUE) - min(values, na.rm = TRUE), 2),
          N = sum(!is.na(values)),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  cat("\n=== DESCRIPTIVE STATISTICS (Simulation-Level, z-aggregated) ===\n\n")
  print(descriptive_stats_table, digits = 2, row.names = FALSE)
  cat("\n")

  cat("Note: Correlation metrics aggregated on Fisher-z scale (consistent with\n")
  cat("      regression analysis), then back-transformed to r-scale. Other metrics\n")
  cat("      aggregated on original scale. All aggregation at simulation level.\n\n")

  # Create formatted tables
  if (require("kableExtra", quietly = TRUE)) {

    descriptive_stats_formatted <- descriptive_stats_table %>%
      select(Metric, Mean, SD, Median, Min, Max, Range) %>%
      mutate(across(where(is.numeric), ~round(., 3)))

    # HTML table
    descriptive_stats_formatted %>%
      kable(caption = "Descriptive Statistics (Simulation-Level, z-aggregated)",
            format = "html",
            align = c("l", rep("r", 6))) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      save_kable(file = "descriptive_statistics_table.html")

    # LaTeX table
    print(
      xtable(descriptive_stats_formatted,
             caption = "Descriptive Statistics (Simulation-Level, z-aggregated)",
             digits = 3),
      file = "descriptive_statistics_table.tex",
      include.rownames = FALSE
    )
  }

  # Store for later use
  list(table = descriptive_stats_table)
}
