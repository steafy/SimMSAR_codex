# =============================================================================
# Statistical Analysis Script for MSAR Models (FINAL VERSION)
# =============================================================================
# 
# RANDOM EFFECTS STRUCTURE:
# This script adaptively selects random effects structure for each outcome:
# 1. Tests whether random slopes for logT improve fit (LR test)
# 2. If yes, checks for boundary singularity in correlated RE model
# 3. If singular (ρ ≈ 1.0), uses uncorrelated RE to ensure stability
# 4. Otherwise, uses correlated RE (best fit)
#
# This outcome-specific approach follows standard practice (Bates et al., 2015)
# and ensures both optimal fit and numerical stability.
# =============================================================================
# -----------------------------------------------------------------------------
# SETUP: Load required packages
# -----------------------------------------------------------------------------

# Core packages
required_packages <- c("dplyr", "ggplot2", "kableExtra", "xtable", "lmtest", "sandwich",
                       "lme4", "lmerTest", "performance", "dunn.test", "cowplot")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}


# -----------------------------------------------------------------------------
# PART 1: DATA PREPARATION
# -----------------------------------------------------------------------------

# Verify msar_results object
if (!inherits(MSAR_dynamics_list, "msar_results")) {
  stop("MSAR_dynamics_list is not an msar_results object. Please re-run estimate_MSAR().")
}


# -----------------------------------------------------------------------------
# 1.1: Create corr_results with grouping variables
# -----------------------------------------------------------------------------

corr_results <- MSAR_dynamics_list %>%
  rename(
    Timesteps = timesteps,
    Density = density,
    Nodes = nodes,
    Regimes = regimes
  ) %>%
  # Add SimID and RegimeIndex
  group_by(Timesteps, Density, Nodes, Regimes) %>%
  mutate(
    SimID = ts_id,
    RegimeIndex = regime_id,
    N = n_distinct(SimID)
  ) %>%
  ungroup() %>%
  # Create Condition identifier
  mutate(Condition = interaction(Timesteps, Density, Nodes, Regimes, drop = TRUE)) %>%
  # Create SimUID (globally unique simulation identifier)
  mutate(SimUID = interaction(Condition, SimID, drop = TRUE)) %>%
  # Reorder columns for clarity
  select(SimUID, Condition, SimID, Timesteps, Density, Nodes, Regimes, RegimeIndex, N, everything())

# Transform to factors
corr_results <- corr_results %>%
  mutate(
    Timesteps = factor(Timesteps),
    Density = factor(Density),
    Nodes = factor(Nodes),
    Regimes = factor(Regimes),
    RegimeIndex = factor(RegimeIndex)
  )


# -----------------------------------------------------------------------------
# 1.2: Calculate estimation process statistics
# -----------------------------------------------------------------------------

# Create a summary table of unique factor level combinations and their N
aggr_factors <- corr_results %>%
  select(Timesteps, Density, Nodes, Regimes, N) %>%
  distinct()

# Calculate expected total simulations and omissions
n_ts_per_condition <- n_distinct(corr_results$SimID)
n_conditions <- nrow(aggr_factors)
expected_total <- n_ts_per_condition * n_conditions

omissions <- expected_total - sum(aggr_factors$N)
omissions_per <- round(((omissions / expected_total) * 100), 2)

cat("\n=== SELECTION BIAS ANALYSIS ===\n")
cat(sprintf("Total expected observations: %d\n", expected_total))
cat(sprintf("Total successful fits: %d\n", sum(aggr_factors$N)))
cat(sprintf("Total failures: %d (%.2f%%)\n\n", omissions, omissions_per))

# Calculate failure rates by condition
failure_analysis <- aggr_factors %>%
  mutate(
    expected_n = n_ts_per_condition,
    failures = expected_n - N,
    failure_rate = failures / expected_n
  ) %>%
  arrange(desc(failure_rate))

# Report conditions with highest failure rates
has_any_failures <- any(failure_analysis$failure_rate > 0)

if (has_any_failures) {
  cat("Conditions with failures (sorted by failure rate):\n")
  high_failure <- failure_analysis %>% filter(failure_rate > 0)
  print(high_failure[, c("Timesteps", "Density", "Nodes", "Regimes", "N", "failures", "failure_rate")],
        n = min(20, nrow(high_failure)))
  cat("\n")
  cat("Systematic bias test reported after the N-distribution analysis below.\n\n")
} else {
  cat("✓ No failures detected - all estimations successful.\n\n")
  selection_bias_detected <- FALSE
  selection_bias_info <- NULL
}

# Build condition-level success rates for sensitivity analyses.
selection_weights <- aggr_factors %>%
  mutate(
    expected_n = n_ts_per_condition,
    success_rate = pmax(N / expected_n, 1e-6)  # avoid division by zero
  ) %>%
  mutate(
    Timesteps_chr = as.character(Timesteps),
    Density_chr = as.character(Density),
    Nodes_chr = as.character(Nodes),
    Regimes_chr = as.character(Regimes)
  ) %>%
  select(Timesteps_chr, Density_chr, Nodes_chr, Regimes_chr, success_rate)

cat("Condition-level success rate summary (for sensitivity analyses):\n")
print(summary(selection_weights$success_rate))
cat("\n")


# Calculate mean no. of estimated models (N) for factorlevels
n_means <- list()

# Calculate dunn-test to compare N across factorlevels
dunn_results <- list()

# First check: Is there any variance in N at all?
if (length(unique(aggr_factors$N)) == 1) {
  message(sprintf("Note: All conditions have the same number of successful fits (N = %d).",
                  unique(aggr_factors$N)))
  message("Kruskal-Wallis and Dunn tests are not applicable - no variance to test.")
  
  # Still calculate means
  for (i in 1:4) {
    var <- colnames(aggr_factors)[i]
    val <- unique(aggr_factors[[var]])
    for (j in seq_along(val)) {
      x <- val[j]
      subset <- aggr_factors %>% filter(.data[[var]] == x)
      n_means[[paste(x, var)]] <- mean(subset$N)
    }
    dunn_results[[var]] <- list(
      Kruskal = "Not applicable - no variance in N",
      Dunn = "Not applicable - no variance in N"
    )
  }
} else {
  # There is variance - proceed with tests
  for (i in 1:4) {
    var <- colnames(aggr_factors)[i]
    val <- unique(aggr_factors[[var]])
    
    # Calculate means for each level
    for (j in seq_along(val)) {
      x <- val[j]
      subset <- aggr_factors %>% filter(.data[[var]] == x)
      n_means[[paste(x, var)]] <- mean(subset$N)
    }
    
    # Only perform statistical tests if there are multiple groups to compare
    if (length(val) > 1) {
      # Check if there's variance within this specific factor
      group_means <- tapply(aggr_factors$N, aggr_factors[[var]], mean)
      if (length(unique(group_means)) == 1) {
        message(sprintf("Skipping tests for '%s' - all groups have same mean N (%.1f)",
                        var, group_means[1]))
        dunn_results[[var]] <- list(
          Kruskal = "Not applicable - no between-group variance",
          Dunn = "Not applicable - no between-group variance"
        )
        next
      }
      
      # Run tests with error handling
      kw_result <- tryCatch({
        kruskal.test(aggr_factors$N, aggr_factors[[var]])
      }, error = function(e) {
        message(sprintf("Error in Kruskal-Wallis for '%s': %s", var, e$message))
        return(NULL)
      }, warning = function(w) {
        message(sprintf("Warning in Kruskal-Wallis for '%s': %s", var, w$message))
        return(NULL)
      })
      
      dunn_result <- tryCatch({
        dunn.test(aggr_factors$N, aggr_factors[[var]], method = "bonferroni")
      }, error = function(e) {
        message(sprintf("Error in Dunn test for '%s': %s", var, e$message))
        return(NULL)
      }, warning = function(w) {
        message(sprintf("Warning in Dunn test for '%s': %s", var, w$message))
        return(NULL)
      })
      
      # Store results (even if NULL, for debugging)
      if (!is.null(kw_result) && !is.null(dunn_result)) {
        dunn_matrix <- data.frame(
          comp = dunn_result$comparisons,
          Z_val = dunn_result$Z,
          p_val = dunn_result$P,
          p.adj = round(dunn_result$P.adjusted, 4)
        )
        dunn_matrix <- dunn_matrix[order(dunn_matrix$p.adj), ]
        dunn_results[[var]] <- list(
          Kruskal = kw_result,
          Dunn = dunn_matrix
        )
      } else {
        dunn_results[[var]] <- list(
          Kruskal = if (is.null(kw_result)) "Test failed" else kw_result,
          Dunn = if (is.null(dunn_result)) "Test failed" else dunn_result
        )
      }
    } else {
      # Only one group - tests not applicable
      message(sprintf("Skipping Kruskal-Wallis and Dunn tests for '%s' - only one group (%s)",
                      var, val))
      dunn_results[[var]] <- list(
        Kruskal = "Not applicable - only one group",
        Dunn = "Not applicable - only one group"
      )
    }
  }
}

# -----------------------------------------------------------------------------
# PART 1.3: SELECTION BIAS WARNING (reuses KW results from N-distribution above)
# -----------------------------------------------------------------------------

if (has_any_failures) {
  cat("Testing if failure rates differ systematically by predictor (KW on N):\n")
  bias_p_values <- list()
  for (predictor in c("Timesteps", "Density", "Nodes", "Regimes")) {
    kw <- dunn_results[[predictor]]$Kruskal
    if (is.list(kw) && !is.null(kw$p.value)) {
      p_val <- kw$p.value
      cat(sprintf("  %s: χ²(%.0f) = %.3f, p = %.4f %s\n",
                  predictor, kw$parameter, kw$statistic, p_val,
                  ifelse(p_val < 0.05, "**SIGNIFICANT**", "")))
    } else {
      p_val <- 1.0
      cat(sprintf("  %s: %s\n", predictor,
                  if (is.character(kw)) kw else "test failed"))
    }
    bias_p_values[[predictor]] <- p_val
  }
  cat("\n")
  
  systematic_bias <- any(unlist(bias_p_values) < 0.05)
  biased_predictors <- names(bias_p_values)[unlist(bias_p_values) < 0.05]
  
  if (systematic_bias) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("⚠⚠⚠ CRITICAL WARNING: SELECTION BIAS DETECTED ⚠⚠⚠\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")
    
    cat(sprintf("Failure rates differ significantly for: %s\n\n",
                paste(biased_predictors, collapse = ", ")))
    
    cat("IMPLICATIONS:\n")
    cat("  • Data are NOT missing completely at random (MCAR)\n")
    cat("  • Statistical models below use ONLY successful estimations\n")
    cat("  • Mixed models below are fit without IPW (lmer weights are variance weights,\n")
    cat("    not sampling/IPW weights)\n")
    cat("  • Residual bias may remain if missingness depends on unmodeled factors\n")
    cat(sprintf("  • Effects of %s may be particularly biased\n\n",
                paste(biased_predictors, collapse = ", ")))
    
    cat("THIS LIMITATION IS NOT FULLY ADDRESSED IN THE PRIMARY MODELS\n\n")
    
    cat("REQUIRED ACTIONS:\n")
    cat("  1. Report this limitation prominently in your Methods/Limitations\n")
    cat("  2. Interpret results for affected predictors with caution\n")
    cat("  3. Consider sensitivity analysis (see below)\n\n")
    
    cat("SENSITIVITY ANALYSIS OPTION:\n")
    cat("  Re-run analysis excluding conditions with failure rate > 10%%:\n")
    high_fail_conditions <- failure_analysis %>%
      filter(failure_rate > 0.1) %>%
      select(Timesteps, Density, Nodes, Regimes)
    
    if (nrow(high_fail_conditions) > 0) {
      cat("  Exclude these conditions:\n")
      print(high_fail_conditions, n = min(10, nrow(high_fail_conditions)))
      cat("\n  Then check if main conclusions change.\n")
    } else {
      cat("  (No conditions with >10%% failure rate)\n")
    }
    
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")
    
    selection_bias_detected <- TRUE
    selection_bias_info <- list(
      biased_predictors = biased_predictors,
      p_values = bias_p_values,
      high_failure_conditions = high_fail_conditions
    )
  } else {
    cat("✓ No systematic selection bias detected (failures appear random).\n\n")
    selection_bias_detected <- FALSE
    selection_bias_info <- NULL
  }
}

# -----------------------------------------------------------------------------
# PART 2: FISHER-Z TRANSFORMATION
# -----------------------------------------------------------------------------

corr_cols <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")
eps <- 1e-6

for (cc in corr_cols) {
  r <- corr_results[[cc]]
  # Clip extreme values to avoid Inf/-Inf
  r <- pmin(pmax(r, -1 + eps), 1 - eps)
  # Fisher-z transformation: z = atanh(r)
  corr_results[[paste0(cc, "_z")]] <- atanh(r)
}

# Check for extreme values
z_cols <- paste0(corr_cols, "_z")
for (zc in z_cols) {
  n_extreme <- sum(abs(corr_results[[zc]]) > 5)
  if (n_extreme > 0) {
    cat("  ⚠", zc, "has", n_extreme, "extreme values (|z| > 5)\n")
  }
}

# -----------------------------------------------------------------------------
# PART 3: AGGREGATION TO SIMULATION-LEVEL
# -----------------------------------------------------------------------------

# Convert factors to numeric for aggregation
corr_results_for_agg <- corr_results %>%
  mutate(
    Timesteps_num = as.numeric(as.character(Timesteps)),
    Density_num = as.numeric(as.character(Density)),
    Nodes_num = as.numeric(as.character(Nodes)),
    Regimes_num = as.numeric(as.character(Regimes))
  )

# Aggregate: Mean across RegimeIndex for each SimUID
# CORRELATIONS: Aggregate on z-scale
dat_sim <- corr_results_for_agg %>%
  group_by(SimUID, Condition, SimID, 
           Timesteps_num, Density_num, Nodes_num, Regimes_num) %>%
  summarise(
    across(all_of(z_cols), \(x) mean(x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  ) %>%
  rename(
    Timesteps = Timesteps_num,
    Density = Density_num,
    Nodes = Nodes_num,
    Regimes_fac = Regimes_num  
  )

# OTHER METRICS: Aggregate on original scale
other_metrics <- c("MAE_Wtemp", "MAE_Wcont", 
                   "Wtemp_sen", "Wcont_sen",
                   "Wtemp_spec", "Wcont_spec")

available_other_metrics <- other_metrics[other_metrics %in% names(corr_results_for_agg)]

if (length(available_other_metrics) > 0) {
  cat("  Aggregating other metrics:", paste(available_other_metrics, collapse = ", "), "\n")
  
  dat_sim_other <- corr_results_for_agg %>%
    group_by(SimUID, Condition, SimID,
             Timesteps_num, Density_num, Nodes_num, Regimes_num) %>%
    summarise(
      across(all_of(available_other_metrics), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    select(SimUID, all_of(available_other_metrics))
  
  # Merge with dat_sim
  dat_sim <- dat_sim %>%
    left_join(dat_sim_other, by = "SimUID")
}

# Check: weighting problems solved?
original_regime_dist <- table(corr_results_for_agg$Regimes_num)
cat("Original distribution (row-level):\n")
print(original_regime_dist)
cat("\nPercentages:\n")
print(round(100 * original_regime_dist / sum(original_regime_dist), 1))

sim_regime_dist <- table(dat_sim$Regimes_fac)
cat("\nAggregated distribution (simulation-level):\n")
print(sim_regime_dist)
cat("\nPercentages:\n")
print(round(100 * sim_regime_dist / sum(sim_regime_dist), 1))
cat("\n")



# -----------------------------------------------------------------------------
# PART 3.1: CALCULATE DESCRIPTIVE STATISTICS (z-aggregated, then back-transformed)
# -----------------------------------------------------------------------------

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
descript_stats <- list(
  table = descriptive_stats_table
)


# -----------------------------------------------------------------------------
# PART 4: TRANSFORM & SCALE PREDICTORS
# -----------------------------------------------------------------------------

dat_sim <- dat_sim %>%
  mutate(
    # log(Timesteps) for diminishing returns
    logT = as.numeric(scale(log(Timesteps))),
    
    # Density and Nodes: centered and scaled
    Density_s = as.numeric(scale(Density)),
    Nodes_s = as.numeric(scale(Nodes)),
    
    # Regimes as FACTOR (not scaled!)
    Regimes = factor(Regimes_fac),
    
    # Create network_id: unique identifier for each network (same network measured at different T)
    # NOTE: Does NOT include Timesteps because same network is used across T values
    network_id = interaction(Density, Nodes, Regimes, SimID, drop = TRUE)
  ) %>%
  select(-Regimes_fac)  # Remove temporary variable

# Attach condition-level success rates to each retained simulation row.
dat_sim <- dat_sim %>%
  mutate(
    Timesteps_chr = as.character(Timesteps),
    Density_chr = as.character(Density),
    Nodes_chr = as.character(Nodes),
    Regimes_chr = as.character(Regimes)
  ) %>%
  left_join(
    selection_weights,
    by = c("Timesteps_chr", "Density_chr", "Nodes_chr", "Regimes_chr")
  ) %>%
  mutate(
    success_rate = ifelse(is.na(success_rate), 1, success_rate)
  ) %>%
  select(-Timesteps_chr, -Density_chr, -Nodes_chr, -Regimes_chr)



# -----------------------------------------------------------------------------
# PART 5: FIT MODELS (MAIN, 2-WAY, 3-WAY) - MIXED EFFECTS
# -----------------------------------------------------------------------------

cat("\n=== Fitting Mixed-Effects Models ===\n")
cat("Note: Testing random slopes for logT to determine appropriate random-effects structure\n\n")

all_results <- list()

for (outcome in corr_cols) {
  
  outcome_z <- paste0(outcome, "_z")
  
  cat("Fitting models for:", outcome, "\n")
  cat(strrep("-", 60), "\n")
  
  # 5.0: Test random slopes for logT
  # Different networks may have different trajectories across T
  cat("Testing random-effects structure...\n")
  
  # Model with random intercept only
  formula_ri <- as.formula(
    paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes + (1|network_id)")
  )
  model_ri <- lmer(formula_ri, data = dat_sim, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa"))
  
  # Model with random intercept + random slope for logT (correlated)
  formula_rs_corr <- as.formula(
    paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes + (1 + logT|network_id)")
  )
  model_rs_corr <- tryCatch({
    lmer(formula_rs_corr, data = dat_sim, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa"))
  }, error = function(e) {
    cat("  Random slopes (correlated) model failed to converge\n")
    return(NULL)
  })
  
  # Model with random intercept + random slope for logT (uncorrelated)
  formula_rs_uncorr <- as.formula(
    paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes + (1|network_id) + (0 + logT|network_id)")
  )
  model_rs_uncorr <- tryCatch({
    lmer(formula_rs_uncorr, data = dat_sim, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa"))
  }, error = function(e) {
    cat("  Random slopes (uncorrelated) model failed to converge\n")
    return(NULL)
  })
  
  # Compare models and choose structure
  use_random_slopes <- FALSE
  use_uncorrelated <- FALSE
  
  if (!is.null(model_rs_corr)) {
    lr_test <- anova(model_ri, model_rs_corr)
    p_value <- lr_test$`Pr(>Chisq)`[2]
    is_singular <- isSingular(model_rs_corr)
    
    cat(sprintf("  Random Intercept AIC: %.1f\n", AIC(model_ri)))
    cat(sprintf("  Random Slopes (Correlated) AIC: %.1f\n", AIC(model_rs_corr)))
    cat(sprintf("  LR test: χ²(%.0f) = %.2f, p = %.4f\n",
                lr_test$Df[2], lr_test$Chisq[2], p_value))
    cat(sprintf("  Singular: %s\n", is_singular))
    
    # Extract correlation if available
    if (!is_singular) {
      vc <- VarCorr(model_rs_corr)
      if ("network_id" %in% names(vc)) {
        corr_mat <- attr(vc$network_id, "correlation")
        if (!is.null(corr_mat) && length(corr_mat) >= 4) {
          rho <- corr_mat[1, 2]
          cat(sprintf("  Intercept-Slope correlation: %.3f\n", rho))
        }
      }
    }
    
    if (p_value < 0.05) {
      use_random_slopes <- TRUE
      
      # Check for singularity - if singular, use uncorrelated
      if (is_singular) {
        cat("  → Boundary singularity detected! Testing uncorrelated random effects...\n")
        
        if (!is.null(model_rs_uncorr)) {
          cat(sprintf("  Random Slopes (Uncorrelated) AIC: %.1f\n", AIC(model_rs_uncorr)))
          cat(sprintf("  ΔAIC (Corr vs Uncorr): %.1f\n", AIC(model_rs_corr) - AIC(model_rs_uncorr)))
          use_uncorrelated <- TRUE
          cat("  → Using uncorrelated random slopes (avoids singularity)\n\n")
        } else {
          cat("  → Uncorrelated model failed, falling back to random intercept only\n\n")
          use_random_slopes <- FALSE
        }
      } else {
        cat("  → Using correlated random slopes (significantly better fit, no singularity)\n\n")
      }
    } else {
      cat("  → Using random intercept only (slopes not justified)\n\n")
    }
  } else {
    cat("  → Using random intercept only (random slopes failed to fit)\n\n")
  }
  
  # Select random effects formula based on tests
  if (!use_random_slopes) {
    re_formula <- "(1|network_id)"
  } else if (use_uncorrelated) {
    re_formula <- "(1|network_id) + (0 + logT|network_id)"
  } else {
    re_formula <- "(1 + logT|network_id)"
  }
  
  # 5.1: Main effects model (baseline)
  formula_main <- as.formula(
    paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes + ", re_formula)
  )
  
  model_main <- lmer(formula_main, data = dat_sim, REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
  
  # 5.2: 2-way interactions (PRIMARY MODEL)
  formula_2way <- as.formula(
    paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2 + ", re_formula)
  )
  
  model_2way <- lmer(formula_2way, data = dat_sim, REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
  
  # 5.2.1: CHECK FOR SINGULARITY IN 2-WAY MODEL AND REFIT IF NEEDED
  # The 2-way model may be singular even if the main effects test was not
  if (use_random_slopes && !use_uncorrelated && isSingular(model_2way)) {
    cat("\n⚠ SINGULARITY DETECTED in 2-way model!\n")
    cat("  The 2-way model is singular even though main effects test was not.\n")
    cat("  Refitting all models with uncorrelated random effects...\n\n")
    
    # Update flags and formula
    use_uncorrelated <- TRUE
    re_formula <- "(1|network_id) + (0 + logT|network_id)"
    
    # Refit main effects model
    formula_main <- as.formula(
      paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes + ", re_formula)
    )
    model_main <- lmer(formula_main, data = dat_sim, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa"))
    
    # Refit 2-way model
    formula_2way <- as.formula(
      paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2 + ", re_formula)
    )
    model_2way <- lmer(formula_2way, data = dat_sim, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa"))
    
    cat("  ✓ Models refitted with uncorrelated random effects\n")
    cat(sprintf("  New formula: %s\n\n", re_formula))
  }
  
  # 5.3: 3-way interactions
  formula_3way <- as.formula(
    paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^3 + ", re_formula)
  )
  
  model_3way <- lmer(formula_3way, data = dat_sim, REML = FALSE,
                     control = lmerControl(optimizer = "bobyqa"))
  
  # 5.3.1: CHECK FOR SINGULARITY IN 3-WAY MODEL AND REFIT IF NEEDED
  if (use_random_slopes && !use_uncorrelated && isSingular(model_3way)) {
    cat("\n⚠ SINGULARITY DETECTED in 3-way model!\n")
    cat("  Refitting all models with uncorrelated random effects...\n\n")
    
    use_uncorrelated <- TRUE
    re_formula <- "(1|network_id) + (0 + logT|network_id)"
    
    formula_main <- as.formula(
      paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes + ", re_formula)
    )
    model_main <- lmer(formula_main, data = dat_sim, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa"))
    
    formula_2way <- as.formula(
      paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2 + ", re_formula)
    )
    model_2way <- lmer(formula_2way, data = dat_sim, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa"))
    
    formula_3way <- as.formula(
      paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^3 + ", re_formula)
    )
    model_3way <- lmer(formula_3way, data = dat_sim, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa"))
    
    cat("  ✓ All models refitted with uncorrelated random effects\n")
    cat(sprintf("  New formula: %s\n\n", re_formula))
  }
  
  # 5.3.1: CONVERGENCE AND SINGULARITY CHECKS
  cat("\nModel diagnostics:\n")
  
  # Check main model
  if (isSingular(model_main)) {
    cat("  ⚠ WARNING: Main effects model has singular fit (overparameterized)\n")
    cat("    → Random effects variance may be degenerate\n")
  }
  if (length(model_main@optinfo$conv$lme4$messages) > 0) {
    cat("  ⚠ WARNING: Main effects model convergence issues:\n")
    cat(paste0("    ", model_main@optinfo$conv$lme4$messages, collapse = "\n"), "\n")
  }
  
  # Check 2-way model
  if (isSingular(model_2way)) {
    cat("  ⚠ WARNING: 2-way model has singular fit (overparameterized)\n")
    cat("    → Random effects variance may be degenerate\n")
    cat("    → Consider simplifying random effects structure\n")
  }
  conv_msgs_2way <- model_2way@optinfo$conv$lme4$messages
  if (length(conv_msgs_2way) > 0) {
    cat("  ⚠ WARNING: 2-way model convergence issues:\n")
    cat(paste0("    ", conv_msgs_2way, collapse = "\n"), "\n")
  }
  
  # Check 3-way model
  if (isSingular(model_3way)) {
    cat("  ⚠ WARNING: 3-way model has singular fit (overparameterized)\n")
    cat("    → Random effects variance may be degenerate\n")
  }
  conv_msgs_3way <- model_3way@optinfo$conv$lme4$messages
  if (length(conv_msgs_3way) > 0) {
    cat("  ⚠ WARNING: 3-way model convergence issues:\n")
    cat(paste0("    ", conv_msgs_3way, collapse = "\n"), "\n")
  }
  
  # Report if all models OK
  all_ok <- !isSingular(model_main) && !isSingular(model_2way) && !isSingular(model_3way) &&
    length(model_main@optinfo$conv$lme4$messages) == 0 &&
    length(conv_msgs_2way) == 0 &&
    length(conv_msgs_3way) == 0
  
  if (all_ok) {
    cat("  ✓ All models converged successfully without singularity\n")
  } else {
    cat("\n  NOTE: Singularity or convergence warnings detected.\n")
    cat("  Results should be interpreted with caution.\n")
    cat("  Consider simplifying model or increasing sample size.\n")
  }
  cat("\n")
  
  # 5.4: Model comparison
  # Get R² values using performance package
  r2_main <- performance::r2_nakagawa(model_main)
  r2_2way <- performance::r2_nakagawa(model_2way)
  r2_3way <- performance::r2_nakagawa(model_3way)
  
  # Likelihood ratio tests (appropriate for nested mixed models)
  anova_2v1 <- anova(model_main, model_2way)
  anova_3v2 <- anova(model_2way, model_3way)
  
  # Primary model selection based on ΔAIC (threshold: 10)
  delta_aic_3way <- AIC(model_2way) - AIC(model_3way)  # positive => 3-way better
  
  if (delta_aic_3way > 10) {
    primary_model   <- model_3way
    primary_label   <- "3-Way (primary)"
    secondary_label <- "2-Way (comparison)"
    cat(sprintf("\n→ 3-Way model selected as primary (ΔAIC = %.1f)\n", delta_aic_3way))
  } else {
    primary_model   <- model_2way
    primary_label   <- "2-Way (primary)"
    secondary_label <- "3-Way (comparison)"
    cat(sprintf("\n→ 2-Way model selected as primary (ΔAIC = %.1f, below threshold of 10)\n",
                delta_aic_3way))
  }
  
  comparison <- data.frame(
    Model = c("Main Effects", secondary_label, primary_label),
    df_fixed = c(
      length(fixef(model_main)),
      length(fixef(model_2way)),
      length(fixef(model_3way))
    ),
    AIC = c(AIC(model_main), AIC(model_2way), AIC(model_3way)),
    BIC = c(BIC(model_main), BIC(model_2way), BIC(model_3way)),
    R2_marginal = c(
      r2_main$R2_marginal,
      r2_2way$R2_marginal,
      r2_3way$R2_marginal
    ),
    R2_conditional = c(
      r2_main$R2_conditional,
      r2_2way$R2_conditional,
      r2_3way$R2_conditional
    )
  )
  
  comparison$Delta_AIC <- comparison$AIC - min(comparison$AIC)
  comparison$Delta_BIC <- comparison$BIC - min(comparison$BIC)
  
  # Print random effects variance for primary model
  cat("\nRandom Effects (primary model):\n")
  vc <- VarCorr(primary_model)
  network_var <- as.numeric(vc$network_id[1])
  residual_var <- attr(vc, "sc")^2
  icc <- network_var / (network_var + residual_var)
  cat(sprintf("  Network variance: %.4f\n", network_var))
  cat(sprintf("  Residual variance: %.4f\n", residual_var))
  cat(sprintf("  ICC: %.3f (%.1f%% of variance is between-network)\n", icc, icc*100))
  cat("\n")
  
  # Store results
  all_results[[outcome]] <- list(
    main        = model_main,
    two_way     = model_2way,
    three_way   = model_3way,
    comparison  = comparison,
    anova = list(
      two_vs_main  = anova_2v1,
      three_vs_two = anova_3v2
    ),
    primary_model          = primary_model,
    primary_label          = primary_label,
    icc                    = icc,
    random_effects_formula = re_formula,
    used_random_slopes     = use_random_slopes,
    used_uncorrelated      = use_uncorrelated,
    convergence_ok         = all_ok,
    singular_fit           = isSingular(primary_model)
  )
}

# -----------------------------------------------------------------------------
# PART 5.5: SENSITIVITY ANALYSIS — exclude high-failure conditions (>10%)
# -----------------------------------------------------------------------------

if (exists("selection_bias_detected") && selection_bias_detected &&
    exists("high_fail_conditions") && nrow(high_fail_conditions) > 0) {
  
  cat("\n=== SENSITIVITY ANALYSIS: Excluding high-failure conditions (>10%) ===\n\n")
  cat("Excluded conditions:\n")
  print(high_fail_conditions)
  cat("\n")
  
  # Build exclusion filter on dat_sim (numeric Timesteps/Density/Nodes, factor Regimes)
  excl <- high_fail_conditions %>%
    mutate(
      Timesteps_num = as.numeric(as.character(Timesteps)),
      Density_num   = as.numeric(as.character(Density)),
      Nodes_num     = as.numeric(as.character(Nodes)),
      Regimes_chr   = as.character(Regimes)
    ) %>%
    select(Timesteps_num, Density_num, Nodes_num, Regimes_chr)
  
  dat_sim_sens <- dat_sim %>%
    mutate(Regimes_chr = as.character(Regimes)) %>%
    anti_join(excl,
              by = c("Timesteps" = "Timesteps_num",
                     "Density"   = "Density_num",
                     "Nodes"     = "Nodes_num",
                     "Regimes_chr")) %>%
    select(-Regimes_chr)
  
  cat(sprintf("Rows in full dat_sim:        %d\n", nrow(dat_sim)))
  cat(sprintf("Rows after exclusion:        %d\n", nrow(dat_sim_sens)))
  cat(sprintf("Rows removed:                %d\n\n", nrow(dat_sim) - nrow(dat_sim_sens)))
  
  sens_results <- list()
  
  for (outcome in corr_cols) {
    outcome_z   <- paste0(outcome, "_z")
    re_formula  <- all_results[[outcome]]$random_effects_formula
    primary_lbl <- all_results[[outcome]]$primary_label
    
    # Mirror the primary model's fixed-effects structure
    if (grepl("3-Way", primary_lbl)) {
      sens_formula <- as.formula(
        paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^3 + ", re_formula)
      )
    } else {
      sens_formula <- as.formula(
        paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2 + ", re_formula)
      )
    }
    
    model_sens <- tryCatch(
      lmer(sens_formula, data = dat_sim_sens, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) { cat("  Sensitivity model failed for", outcome, ":", e$message, "\n"); NULL }
    )
    
    if (is.null(model_sens)) next
    
    # Compare primary vs sensitivity coefficients
    coef_full <- fixef(all_results[[outcome]]$primary_model)
    coef_sens <- fixef(model_sens)
    common    <- intersect(names(coef_full), names(coef_sens))
    
    comparison <- data.frame(
      Predictor  = common,
      Full       = round(coef_full[common], 4),
      Sensitivity = round(coef_sens[common], 4),
      Delta      = round(coef_sens[common] - coef_full[common], 4),
      Delta_pct  = round(100 * (coef_sens[common] - coef_full[common]) /
                           ifelse(abs(coef_full[common]) < 1e-10, NA, coef_full[common]), 1),
      row.names  = NULL
    )
    comparison$Predictor <- gsub(":", " × ", comparison$Predictor)
    
    cat(sprintf("Outcome: %s\n", outcome))
    cat(strrep("-", 70), "\n")
    print(comparison, row.names = FALSE)
    cat("\n")
    
    sens_results[[outcome]] <- list(model = model_sens, comparison = comparison)
  }
  
  cat("Delta = Sensitivity - Full. Large |Delta| or |Delta_pct| indicate bias.\n\n")
  
} else {
  cat("\nNo sensitivity analysis needed (no high-failure conditions or no bias detected).\n\n")
  sens_results <- NULL
}

# -----------------------------------------------------------------------------
# PART 6: EXPORT RESULTS
# -----------------------------------------------------------------------------

cat("\n=== Exporting Results ===\n\n")

# 7.1: Coefficients for primary models - Mixed Effects Version
for (outcome in corr_cols) {
  
  model       <- all_results[[outcome]]$primary_model
  model_label <- all_results[[outcome]]$primary_label  # e.g. "3-Way (primary)"
  
  # Derive a clean file tag: "3way" or "2way"
  model_tag <- ifelse(grepl("3-Way", model_label), "3way", "2way")
  
  # Get coefficient summary from lmerTest (includes p-values via Satterthwaite approximation)
  coef_summary <- summary(model)$coefficients
  
  coef_df <- data.frame(
    Predictor = rownames(coef_summary),
    Estimate = coef_summary[, "Estimate"],
    `Std. Error` = coef_summary[, "Std. Error"],
    df = coef_summary[, "df"],
    `t value` = coef_summary[, "t value"],
    p_raw = coef_summary[, "Pr(>|t|)"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # FDR correction within each outcome model
  coef_df$`p (BH-adj.)` <- p.adjust(coef_df$p_raw, method = "BH")
  
  # Significance stars based on FDR-adjusted p-values
  coef_df$sig <- ifelse(coef_df$`p (BH-adj.)` < 0.001, "***",
                        ifelse(coef_df$`p (BH-adj.)` < 0.01, "**",
                               ifelse(coef_df$`p (BH-adj.)` < 0.05, "*", "")))
  
  # Final table: only adjusted p-values
  coef_df_export <- coef_df[, c("Predictor", "Estimate", "Std. Error", "df", "t value",
                                "p (BH-adj.)", "sig")]
  
  # Clean predictor names: interaction symbol + remove standardisation suffixes
  coef_df_export$Predictor <- gsub(":",         " × ", coef_df_export$Predictor)
  coef_df_export$Predictor <- gsub("_s\\b",     "",    coef_df_export$Predictor)
  coef_df_export$Predictor <- gsub("_num\\b",   "",    coef_df_export$Predictor)
  
  # File names reflect actual primary model
  file_base <- paste0("Results_", model_tag, "_", outcome, "_mixed")
  
  # Export to HTML
  kable(coef_df_export,
        caption = paste0(model_label, ": ", outcome, " (Mixed Effects)"),
        digits = 3,
        format = "html",
        row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    save_kable(file = paste0(file_base, ".html"))
  
  # Export to LaTeX (digits=3)
  print(
    xtable(coef_df_export,
           caption = paste0(model_label, ": ", outcome, " (Mixed Effects)"),
           digits = 3,
           row.names = FALSE),
    file = paste0(file_base, ".tex"),
    include.rownames = FALSE
  )
  
  cat("Exported results for", outcome, "(", model_label, ")\n")
}

# 7.2: Model comparison table (one row per outcome)
# ─────────────────────────────────────────────────────────────────────────────
# Each row summarises the full model-selection process and the selected primary
# model's fit indices for one outcome variable.

# LaTeX-formatted outcome labels (fully written out for readability)
# AC = average controllability; corr. = Pearson correlation with true matrix
outcome_tex_labels <- c(
  Wtemp_corr    = "$W_{\\text{temp}}$ correlation",
  Wtemp_ac_corr = "$\\mathrm{AC}(W_{\\text{temp}})$ correlation",
  Wcont_corr    = "$W_{\\text{cont}}$ correlation",
  Wcont_ac_corr = "$\\mathrm{AC}(W_{\\text{cont}})$ correlation"
)

# Human-readable random-effects structure labels
# RS always refers to a random slope for log-transformed timesteps (log T)
format_re <- function(re) {
  if (grepl("0 \\+ logT", re)) return("RI + RS($\\log T$, uncorr.)")
  if (grepl("1 \\+ logT", re)) return("RI + RS($\\log T$, corr.)")
  return("RI only")
}

comparison_all <- do.call(rbind, lapply(corr_cols, function(outcome) {
  res <- all_results[[outcome]]
  
  # R² for the *primary* model only
  r2p <- performance::r2_nakagawa(res$primary_model)
  
  # ΔAIC: positive values favour the 3-way model
  delta_aic <- round(AIC(res$two_way) - AIC(res$three_way), 1)
  
  # Clean primary-model label
  primary_clean <- ifelse(grepl("3-Way", res$primary_label), "3-Way", "2-Way")
  
  # Flag any convergence / singularity issues with the primary model
  # (used for post-hoc console check; column dropped from LaTeX output)
  conv_note <- dplyr::case_when(
    res$convergence_ok  & !res$singular_fit ~ "OK",
    !res$convergence_ok & res$singular_fit  ~ "Conv. warn.; Singular",
    !res$convergence_ok                     ~ "Conv. warning",
    TRUE                                    ~ "Singular fit"
  )
  
  data.frame(
    Outcome       = outcome_tex_labels[[outcome]],
    RE_structure  = format_re(res$random_effects_formula),
    AIC_main      = round(AIC(res$main),             1),
    AIC_2way      = round(AIC(res$two_way),          1),
    AIC_3way      = round(AIC(res$three_way),        1),
    Primary_model = primary_clean,
    Delta_AIC     = delta_aic,
    R2_m          = round(r2p$R2_marginal,           3),
    R2_c          = round(r2p$R2_conditional,        3),
    ICC           = round(res$icc,                   3),
    Conv_note     = conv_note,
    stringsAsFactors = FALSE
  )
}))

# Console-friendly column names (Conv_note kept for console; dropped in LaTeX)
colnames(comparison_all) <- c(
  "Outcome", "Random Effects", "AIC (Main)", "AIC (2-Way)", "AIC (3-Way)",
  "Primary Model", "dAIC (2w-3w)", "R2_m", "R2_c", "ICC", "Conv./Sing."
)

cat("\n=== MODEL COMPARISON SUMMARY (one row per outcome) ===\n\n")
print(comparison_all, row.names = FALSE)
cat("\nNote: dAIC (2w-3w) = AIC(2-Way) - AIC(3-Way). Positive values favour the 3-Way model.\n")
cat("R2_m = marginal R2; R2_c = conditional R2 (primary model only).\n")

# Summarise convergence status in one sentence (replaces table column)
any_issues <- any(comparison_all[["Conv./Sing."]] != "OK")
if (any_issues) {
  cat("\n⚠ Convergence/singularity issues detected in at least one primary model — see Conv./Sing. column.\n\n")
} else {
  cat("✓ No convergence or singularity problems were observed for any primary model.\n\n")
}

# ── LaTeX export ──────────────────────────────────────────────────────────────
# Requires \usepackage{booktabs} and \usepackage{adjustbox} in the preamble.
# Conv./Sing. column is intentionally omitted from the LaTeX table:
#   - If all models are OK, a single sentence in the caption suffices.
#   - If issues exist, replace the caption sentence with a specific note.

# Caption sentence on convergence (auto-generated)
conv_caption_sentence <- if (!any_issues) {
  "No convergence warnings or singular fits were observed for any primary model."
} else {
  issues <- comparison_all[comparison_all[["Conv./Sing."]] != "OK",
                           c("Outcome", "Conv./Sing.")]
  paste0("Note: ",
         paste(apply(issues, 1, function(r) paste0(r[1], " (", r[2], ")")),
               collapse = "; "), ".")
}

make_tex_row <- function(row) {
  paste0(
    "    ",
    row["Outcome"],          " & ",
    row["Random Effects"],   " & ",
    row["AIC (Main)"],       " & ",
    row["AIC (2-Way)"],      " & ",
    row["AIC (3-Way)"],      " & ",
    row["Primary Model"],    " & ",
    row["dAIC (2w-3w)"],     " & ",
    row["R2_m"],             " & ",
    row["R2_c"],             " & ",
    row["ICC"],              " \\\\"
    # Conv./Sing. intentionally omitted
  )
}

tex_rows <- apply(comparison_all, 1, make_tex_row)

tex_out <- c(
  "% Model Comparison Table --- generated by stat_analysis.R",
  "% Requires in preamble: \\usepackage{booktabs}  \\usepackage{adjustbox}",
  "\\begin{table}[htbp]",
  "  \\centering",
  "  \\caption{%",
  "    Model comparison summary across all outcome variables.",
  "    Primary model selection was based on the comparison between the 2-way and 3-way",
  "    interaction models; AIC values for the main-effects model are reported for completeness.",
  "    $\\Delta\\text{AIC}_{\\text{2w}-\\text{3w}} = \\text{AIC}_{\\text{2-way}} - \\text{AIC}_{\\text{3-way}}$;",
  "    positive values favour the 3-way model.",
  "    The 3-way model was selected as primary when $\\Delta\\text{AIC}_{\\text{2w}-\\text{3w}} > 10$.",
  "    $R^2_m$ and $R^2_c$ denote marginal and conditional $R^2$ \\citep{Nakagawa2013}",
  "    of the selected primary model; ICC is the intraclass correlation coefficient.",
  "    RI~=~random intercept; RS~=~random slope for log-transformed timesteps ($\\log T$);",
  "    uncorr.\\ indicates that intercept and slope were modelled without a correlation parameter.",
  paste0("    ", conv_caption_sentence),
  "  }",
  "  \\label{tab:model_comparison}",
  "  \\begin{adjustbox}{max width=\\textwidth}",
  "  \\begin{tabular}{@{}llrrrcrrrr@{}}",
  "    \\toprule",
  "    Outcome & RE structure &",
  "    $\\text{AIC}_{\\text{main}}$ & $\\text{AIC}_{\\text{2w}}$ & $\\text{AIC}_{\\text{3w}}$ &",
  "    Primary & $\\Delta\\text{AIC}_{\\text{2w}-\\text{3w}}$ &",
  "    $R^2_m$ & $R^2_c$ & ICC \\\\",
  "    \\midrule",
  tex_rows,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\end{adjustbox}",
  "\\end{table}"
)

writeLines(tex_out, "comparison_all.tex")
cat("Exported: comparison_all.tex\n\n")


# 7.3: Diagnostics (PDF) - Mixed Model Diagnostics
cat("\nGenerating diagnostic plots...\n")

pdf("diagnostics_all_models_mixed.pdf", width = 12, height = 10)

for (outcome in corr_cols) {
  model <- all_results[[outcome]]$primary_model
  
  cat("  Plotting diagnostics for", outcome, "\n")
  
  # Set up 2x3 layout for mixed model diagnostics
  par(mfrow = c(2, 3))
  
  # Standard lmer diagnostic plots
  plot(model, main = paste(outcome, "- Residuals vs Fitted"))
  plot(model, sqrt(abs(resid(.))) ~ fitted(.),
       main = paste(outcome, "- Scale-Location"),
       ylab = expression(sqrt("|Standardized residuals|")))
  
  # Q-Q plot
  qqnorm(resid(model), main = paste(outcome, "- Q-Q Plot"))
  qqline(resid(model))
  
  # Residuals vs leverage (manual for lmer)
  plot(hatvalues(model), resid(model),
       xlab = "Leverage", ylab = "Residuals",
       main = paste(outcome, "- Residuals vs Leverage"))
  abline(h = 0, lty = 2)
  
  # Random effects Q-Q plot
  re <- ranef(model)$network_id[[1]]
  qqnorm(re, main = paste(outcome, "- Random Effects Q-Q"))
  qqline(re)
  
  # Histogram of residuals
  hist(resid(model), breaks = 50,
       main = paste(outcome, "- Residual Distribution"),
       xlab = "Residuals")
}

dev.off()

cat("Diagnostic plots saved to: diagnostics_all_models_mixed.pdf\n\n")

# -----------------------------------------------------------------------------
# PART 7: Graphical Visualization (Line Plots with Facets)
# -----------------------------------------------------------------------------

# Set variables to plot
cols <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")

for (col_name in cols) {
  
  title <- switch(
    col_name,
    "Wtemp_corr"      = "Mean correlations for Wtemp",
    "Wtemp_ac_corr"   = "Mean correlations for Wtemp average controllability",
    "Wcont_corr"      = "Mean correlations for Wcont",
    "Wcont_ac_corr"   = "Mean correlations for Wcont average controllability",
  )
  
  # Summarize data
  summary_data <- corr_results %>%
    group_by(Timesteps, Density, Nodes, Regimes) %>%
    summarise(
      mean_val = mean(.data[[col_name]], na.rm = TRUE),
      sd_val   = sd(.data[[col_name]], na.rm = TRUE),
      .groups  = "drop"
    )
  
  # Make hover info
  summary_data <- summary_data %>%
    mutate(HoverInfo = paste("Timesteps:", Timesteps,
                             "<br>Nodes:", Nodes,
                             "<br>Mean:", round(mean_val, 2),
                             "<br>Sd:",   round(sd_val, 2)))
  
  # Position for shifted plots
  dodge <- position_dodge(width = 0.5)
  
  # Make plot
  p <- ggplot() +
    # Singular values as scatterplot
    geom_jitter(
      data = corr_results,
      aes(x = Timesteps, y = .data[[col_name]], color = Nodes),
      position = dodge, alpha = 0.2, size = 0.1
    ) +
    # Aggregated data as line plots
    geom_line(
      data = summary_data,
      aes(x = as.numeric(as.character(Timesteps)), y = mean_val, color = Nodes, group = Nodes),
      position = dodge
    ) +
    # Aggregated data as scatterplots
    geom_point(
      data = summary_data,
      aes(x = as.numeric(as.character(Timesteps)), y = mean_val, color = Nodes, text = HoverInfo),
      position = dodge, size = 1
    ) +
    # Facets
    facet_grid(Density ~ Regimes, labeller = label_value) +
    labs(
      title = title,
      x = "Timesteps",
      y = "Mean correlations",
      color = "Nodes"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
      strip.text.y = element_text(vjust = -0.25),
      panel.spacing = unit(0.2, "in"),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Additional labels
  label_x_right <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Density", angle = -90, size = 4, hjust = 0)
  
  label_y_top <- ggplot() +
    theme_void() +
    xlim(0, 1) +
    ylim(0, 1) +
    annotate("text", x = 0.4725, y = 0.5, label = "Regimes", size = 4, hjust = 0.5)
  
  # Combine plots with `cowplot`
  final_plot <- ggdraw() +
    draw_plot(p, 0, 0, 1, 1) +
    draw_plot(label_x_right, 0.96, 0.08, 0.03, 0.8) +
    draw_plot(label_y_top,   0.1,  0.875, 0.84, 0.05)
  
  # Store plot
  output_file <- paste0("Plots/", col_name, "_plot_with_labels.pdf")
  
  ggsave(
    filename = output_file,
    plot     = final_plot,
    width    = 10,
    height   = 8,
    units    = "in",
    dpi      = 600
  )
  
  message("saved: ", output_file)
}

# -----------------------------------------------------------------------------
# FINAL SUMMARY: Model Quality Report
# -----------------------------------------------------------------------------

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
if (exists("selection_bias_detected") && selection_bias_detected) {
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

if (exists("selection_bias_detected") && selection_bias_detected) {
  cat("IMPORTANT: Due to selection bias, consider:\n")
  cat("  • Reporting this limitation in your paper\n")
  cat("  • Conducting sensitivity analysis excluding high-failure conditions\n")
  cat("  • Collecting more data to improve estimation success rates\n\n")
}



# -----------------------------------------------------------------------------
# PART 8: COEFFICIENT PLOT
# -----------------------------------------------------------------------------

# Readable labels for outcomes
outcome_labels <- c(
  Wtemp_corr    = "W_temp (Correlation)",
  Wtemp_ac_corr = "W_temp (Avg. Controllability)",
  Wcont_corr    = "W_cont (Correlation)",
  Wcont_ac_corr = "W_cont (Avg. Controllability)"
)

# Collect coefficients + CIs from all primary models
coef_plot_data <- lapply(corr_cols, function(outcome) {
  model        <- all_results[[outcome]]$primary_model
  coef_summary <- summary(model)$coefficients
  
  data.frame(
    outcome   = outcome_labels[outcome],
    predictor = rownames(coef_summary),
    estimate  = coef_summary[, "Estimate"],
    se        = coef_summary[, "Std. Error"],
    p_adj     = p.adjust(coef_summary[, "Pr(>|t|)"], method = "BH"),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Clean predictor names (same as export)
coef_plot_data$predictor <- gsub(":",    " × ", coef_plot_data$predictor)
coef_plot_data$predictor <- gsub("_s\\b", "",   coef_plot_data$predictor)
coef_plot_data$predictor <- gsub("_num\\b", "",  coef_plot_data$predictor)

# Compute 95% CI and significance flag
coef_plot_data <- coef_plot_data %>%
  mutate(
    ci_lo  = estimate - 1.96 * se,
    ci_hi  = estimate + 1.96 * se,
    sig    = ifelse(p_adj < 0.05, "p < .05", "n.s."),
    # Order predictors: intercept last, main effects first, then interactions by n of ×
    n_cross = stringr::str_count(predictor, "×")
  )

# Exclude intercept from plot (usually not of interest)
coef_plot_data <- coef_plot_data %>%
  filter(predictor != "(Intercept)")

# Order predictors by complexity then alphabetically
predictor_order <- coef_plot_data %>%
  arrange(n_cross, predictor) %>%
  pull(predictor) %>%
  unique()

coef_plot_data$predictor <- factor(coef_plot_data$predictor, levels = rev(predictor_order))

# --- 8.1: Full coefficient plot (all terms) ---
p_coef_full <- ggplot(coef_plot_data,
                      aes(x = estimate, y = predictor,
                          color = outcome, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1)) +
  scale_color_viridis_d(option = "D", end = 0.85) +
  facet_wrap(~ n_cross, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(
               "0" = "Main Effects",
               "1" = "2-Way Interactions",
               "2" = "3-Way Interactions"
             ))) +
  labs(
    x      = "Regression Coefficient (95% CI)",
    y      = NULL,
    color  = "Outcome",
    shape  = "Significance"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "bottom",
    legend.box       = "vertical",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Plots/coefficient_plot_full.pdf",
       plot = p_coef_full, width = 10, height = 14, dpi = 600)

# --- 8.2: Main effects only (compact version for main text) ---
p_coef_main <- coef_plot_data %>%
  filter(n_cross == 0) %>%
  ggplot(aes(x = estimate, y = predictor,
             color = outcome, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1)) +
  scale_color_viridis_d(option = "D", end = 0.85) +
  labs(
    x      = "Regression Coefficient (95% CI)",
    y      = NULL,
    color  = "Outcome",
    shape  = "Significance",
    title  = "Main Effects"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.box       = "vertical",
    panel.grid.minor = element_blank()
  )

ggsave("Plots/coefficient_plot_main_effects.pdf",
       plot = p_coef_main, width = 8, height = 5, dpi = 600)

message("Coefficient plots saved to Plots/")