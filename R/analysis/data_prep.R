# =============================================================================
# Data preparation, estimation-process statistics, and selection-bias reporting
# =============================================================================
# Functions extracted from the former monolithic stat_analysis.R (PARTs 1.1-1.3).
# Each function keeps the original console output so the orchestrator produces
# identical logs when calling them in sequence.

# -----------------------------------------------------------------------------
# 1.1: Create corr_results with grouping variables
# -----------------------------------------------------------------------------
prepare_corr_results <- function(MSAR_dynamics_list) {

  # Verify msar_results object
  if (!inherits(MSAR_dynamics_list, "msar_results")) {
    stop("MSAR_dynamics_list is not an msar_results object. Please re-run estimate_MSAR().")
  }

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

  corr_results
}

# -----------------------------------------------------------------------------
# 1.2: Calculate estimation process statistics (failure rates, KW/Dunn on N)
# -----------------------------------------------------------------------------
# Returns a list with everything PART 1.3 and later sensitivity analyses need.
analyze_estimation_process <- function(corr_results) {

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

  list(
    aggr_factors       = aggr_factors,
    n_ts_per_condition = n_ts_per_condition,
    n_conditions       = n_conditions,
    expected_total     = expected_total,
    failure_analysis   = failure_analysis,
    has_any_failures   = has_any_failures,
    selection_weights  = selection_weights,
    n_means            = n_means,
    dunn_results       = dunn_results
  )
}

# -----------------------------------------------------------------------------
# 1.3: Selection bias warning (reuses KW results from analyze_estimation_process)
# -----------------------------------------------------------------------------
# Returns list(selection_bias_detected, selection_bias_info, high_fail_conditions).
report_selection_bias <- function(est) {

  failure_analysis <- est$failure_analysis
  dunn_results     <- est$dunn_results

  selection_bias_detected <- FALSE
  selection_bias_info     <- NULL
  high_fail_conditions    <- NULL

  if (est$has_any_failures) {
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

  list(
    selection_bias_detected = selection_bias_detected,
    selection_bias_info     = selection_bias_info,
    high_fail_conditions    = high_fail_conditions
  )
}
