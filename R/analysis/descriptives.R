# =============================================================================
# Descriptive statistics (z-aggregated, then back-transformed)
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PART 3.1).

# Returns list(table = descriptive_stats_table) and writes HTML/LaTeX tables.
descriptive_stats <- function(dat_sim, corr_cols, available_other_metrics) {
  
  descriptive_stats_table <- data.frame(
    Metric = character(),
    Mean = numeric(),
    TrimMean = numeric(),   # 5% trimmed mean: robust to skew/outliers
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
          TrimMean = round(mean(r_values, trim = 0.05, na.rm = TRUE), 2),
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
          TrimMean = round(mean(values, trim = 0.05, na.rm = TRUE), 2),
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
  
  # Order rows by outcome block (not by metric type): first the temporal
  # network Beta (correlation, MAE, sensitivity, specificity), then the average
  # controllability of Beta, then the contemporaneous network Kappa. Metrics not
  # listed keep their original order at the end.
  metric_order <- c(
    "Beta_corr", "NRMSE_Beta", "Beta_sen", "Beta_spec",   # Beta block
    "Beta_ac_corr_pearson", "Beta_ac_corr_spearman",      # AC(Beta) block (both variants)
    "Kappa_corr", "NRMSE_Kappa", "Kappa_sen", "Kappa_spec"  # Kappa block
  )
  descriptive_stats_table <- descriptive_stats_table[
    order(match(descriptive_stats_table$Metric, metric_order)), , drop = FALSE
  ]
  rownames(descriptive_stats_table) <- NULL
  
  cat("\n=== DESCRIPTIVE STATISTICS (Simulation-Level, z-aggregated) ===\n\n")
  print(descriptive_stats_table, digits = 2, row.names = FALSE)
  cat("\n")
  
  cat("Note: Correlation metrics aggregated on Fisher-z scale (consistent with\n")
  cat("      regression analysis), then back-transformed to r-scale. Other metrics\n")
  cat("      aggregated on original scale. All aggregation at simulation level.\n\n")
  
  # Create formatted tables
  if (require("kableExtra", quietly = TRUE)) {
    
    # Keep N (the number of non-NA observations that actually entered each
    # metric) in the exported table -- with the Kappa condition-number guard,
    # this now differs across metrics (Kappa metrics have fewer valid obs).
    descriptive_stats_formatted <- descriptive_stats_table %>%
      select(Metric, Mean, TrimMean, SD, Median, Min, Max, Range, N) %>%
      mutate(across(c(Mean, TrimMean, SD, Median, Min, Max, Range), ~round(., 3)))
    
    # HTML table
    descriptive_stats_formatted %>%
      kable(caption = "Descriptive Statistics (Simulation-Level, z-aggregated)",
            format = "html",
            align = c("l", rep("r", 8))) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      save_kable(file = results_file("descriptive_statistics_table.html"))
    
    # LaTeX table: map the Metric column to math labels (subscript-per-row).
    # Recovery correlations use r_<network>; NRMSE/Sensitivity/Specificity
    # carry the network as a subscript. Applied only to the .tex export so
    # console and HTML keep the readable raw names. Requires \usepackage{amsmath}.
    tex_metric_labels <- c(
      Beta_corr    = "$r_{\\text{Beta}}$",
      NRMSE_Beta   = "$\\text{NRMSE}_{\\text{Beta}}$",
      Beta_sen     = "$\\text{Sens}_{\\text{Beta}}$",
      Beta_spec    = "$\\text{Spec}_{\\text{Beta}}$",
      Beta_ac_corr_pearson  = "$r_{\\text{AC, Pearson}}$",
      Beta_ac_corr_spearman = "$r_{\\text{AC, Spearman}}$",
      Kappa_corr   = "$r_{\\text{Kappa}}$",
      NRMSE_Kappa  = "$\\text{NRMSE}_{\\text{Kappa}}$",
      Kappa_sen    = "$\\text{Sens}_{\\text{Kappa}}$",
      Kappa_spec   = "$\\text{Spec}_{\\text{Kappa}}$"
    )
    descriptive_stats_tex <- descriptive_stats_formatted
    mapped <- tex_metric_labels[descriptive_stats_tex$Metric]
    descriptive_stats_tex$Metric <- ifelse(is.na(mapped),
                                           descriptive_stats_tex$Metric, mapped)
    
    # N shown as integer via per-column digits; sanitize.* = identity so the
    # LaTeX math in the Metric column is written verbatim (not escaped).
    print(
      xtable(descriptive_stats_tex,
             caption = "Descriptive Statistics (Simulation-Level, z-aggregated)",
             digits = c(0, 0, 3, 3, 3, 3, 3, 3, 3, 0)),
      file = results_file("descriptive_statistics_table.tex"),
      include.rownames = FALSE,
      sanitize.text.function = identity
    )
  }
  
  # Store for later use
  list(table = descriptive_stats_table)
}