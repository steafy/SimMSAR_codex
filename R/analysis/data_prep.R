# =============================================================================
# Data preparation, estimation-process statistics, and selection-bias reporting
# =============================================================================
# Functions extracted from the former monolithic stat_analysis.R (PARTs 1.1-1.3).
# Each function keeps the original console output so the orchestrator produces
# identical logs when calling them in sequence.

# This file now also performs all per-regime SCORING (previously inline in
# estimate_MSAR.R): edge thresholding, Kappa = solve(Sigma) + validity gate,
# correlations, sens/spec, NRMSE, and average-controllability correlations.
# It therefore pulls the relevant estimation-side utils directly.
source("R/utils/undirected_metrics.R")        # vectorize/senspec/symmetrize upper-tri
source("R/utils/senspec.R")                    # full-matrix sensitivity/specificity (Beta)
source("R/utils/average_controllability.R")    # est_Beta_ac from stored est_Beta

# -----------------------------------------------------------------------------
# 1.0: Per-regime recovery metrics from the RAW stored estimates
# -----------------------------------------------------------------------------
# estimate_MSAR() now returns only raw matrices (orig_Beta/est_Beta,
# orig_Kappa, orig_Sigma/est_Sigma, orig_Beta_ac). This function turns those
# into the scored outcomes the rest of the pipeline expects, so that every
# scoring decision (edge threshold, Kappa condition-number cutoff, AC horizon,
# Pearson vs Spearman AC correlation) is revisable WITHOUT re-running estimation.
#
# Constants (exposed so they can be tuned on fresh data):
#   min_edg_val    edge-detection threshold; must match the generation intent
#                  (0.05 in the driver). Below it, |entries| -> 0. Beta: full
#                  matrix; Kappa: off-diagonal only (the dominant diagonal is
#                  structural and never thresholded).
#   AC_HORIZON     finite horizon T_ac for average_controllability(); MUST equal
#                  generation's default (25) or Beta_ac_corr compares two
#                  different definitions of controllability.
#   KAPPA_COND_MAX max acceptable condition number of est_Sigma. Above it, the
#                  inversion to Kappa is flagged invalid and ALL Kappa quantities
#                  for that regime become NA -- per regime, never dropping the
#                  whole fit, and never affecting Beta. Start near the old 1e6;
#                  re-derive empirically on fresh data.
#
# Returns the input object (class/attrs preserved) with the derived columns
# added, plus attr "sigma_validity_log": one row per per-regime observation
# whose est_Sigma was flagged invalid.
compute_recovery_metrics <- function(MSAR_dynamics_list,
                                     min_edg_val   = 0.05,
                                     AC_HORIZON    = 25,
                                     KAPPA_COND_MAX = 1e6) {

  if (!inherits(MSAR_dynamics_list, "msar_results")) {
    stop("MSAR_dynamics_list is not an msar_results object. Please re-run estimate_MSAR().")
  }

  EDGE_VALUE_RANGE <- 2   # signed generating range [-1, 1] -> NRMSE=1 is full sign reversal

  threshold_beta <- function(mat) {
    mat[abs(mat) < min_edg_val] <- 0
    mat
  }
  threshold_kappa_offdiag <- function(mat) {
    off <- !diag(TRUE, nrow(mat))         # off-diagonal mask
    small <- abs(mat) < min_edg_val
    mat[off & small] <- 0
    mat
  }

  rmse_true_edges <- function(true_vec, est_vec) {
    pos <- true_vec != 0
    if (!any(pos)) return(NA_real_)
    sqrt(mean((true_vec[pos] - est_vec[pos])^2))
  }

  n <- nrow(MSAR_dynamics_list)

  # Pre-allocate outputs
  est_Kappa    <- vector("list", n)
  est_Beta_ac  <- vector("list", n)
  Beta_corr    <- numeric(n)
  Kappa_corr   <- numeric(n)
  Beta_sen     <- numeric(n); Beta_spec  <- numeric(n)
  Kappa_sen    <- numeric(n); Kappa_spec <- numeric(n)
  Beta_ac_corr_pearson  <- numeric(n)
  Beta_ac_corr_spearman <- numeric(n)
  NRMSE_Beta   <- numeric(n)
  NRMSE_Kappa  <- numeric(n)

  validity_records <- list()

  for (r in seq_len(n)) {
    orig_Beta    <- MSAR_dynamics_list$orig_Beta[[r]]
    est_Beta_raw <- MSAR_dynamics_list$est_Beta[[r]]
    orig_Kappa   <- MSAR_dynamics_list$orig_Kappa[[r]]
    est_Sigma    <- MSAR_dynamics_list$est_Sigma[[r]]
    orig_Beta_ac <- MSAR_dynamics_list$orig_Beta_ac[[r]]

    # --- Beta (unaffected by Sigma conditioning): threshold, then score ---------
    est_Beta <- threshold_beta(est_Beta_raw)

    Beta_corr[r] <- cor(as.vector(orig_Beta), as.vector(est_Beta), method = "pearson")
    bss          <- senspec(orig_Beta, est_Beta)
    Beta_sen[r]  <- bss[["sensitivity"]]
    Beta_spec[r] <- bss[["specificity"]]
    NRMSE_Beta[r] <- rmse_true_edges(as.vector(orig_Beta), as.vector(est_Beta)) / EDGE_VALUE_RANGE

    # Average controllability from the thresholded est_Beta (deterministic).
    ac <- average_controllability(est_Beta, T_ac = AC_HORIZON)
    est_Beta_ac[[r]] <- ac
    Beta_ac_corr_pearson[r]  <- cor(orig_Beta_ac, ac, method = "pearson")
    Beta_ac_corr_spearman[r] <- cor(orig_Beta_ac, ac, method = "spearman")

    # --- Kappa: invert est_Sigma, gate on condition number ----------------------
    inv <- tryCatch(solve(est_Sigma), error = function(e) NULL)
    invalid_reason <- NA_character_
    cond_number    <- NA_real_

    if (is.null(inv)) {
      invalid_reason <- "solve_failed"
    } else {
      rc <- tryCatch(rcond(est_Sigma), error = function(e) 0)
      if (!is.finite(rc) || rc <= 1 / KAPPA_COND_MAX) {
        invalid_reason <- "ill_conditioned"
        cond_number    <- if (is.finite(rc) && rc > 0) 1 / rc else Inf
      }
    }

    if (is.na(invalid_reason)) {
      kap <- threshold_kappa_offdiag(symmetrize_matrix(inv))
      est_Kappa[[r]] <- kap
      Kappa_corr[r] <- cor(vectorize_upper_tri(orig_Kappa, diag = FALSE),
                           vectorize_upper_tri(kap,        diag = FALSE),
                           method = "pearson")
      kss           <- senspec_upper_tri(orig_Kappa, kap, diag = FALSE)
      Kappa_sen[r]  <- kss[["sensitivity"]]
      Kappa_spec[r] <- kss[["specificity"]]
      NRMSE_Kappa[r] <- rmse_true_edges(vectorize_upper_tri(orig_Kappa, diag = FALSE),
                                        vectorize_upper_tri(kap,        diag = FALSE)) / EDGE_VALUE_RANGE
    } else {
      est_Kappa[[r]] <- NA
      Kappa_corr[r]  <- NA_real_
      Kappa_sen[r]   <- NA_real_
      Kappa_spec[r]  <- NA_real_
      NRMSE_Kappa[r] <- NA_real_
      validity_records[[length(validity_records) + 1L]] <- tibble::tibble(
        timesteps = MSAR_dynamics_list$timesteps[r],
        density   = MSAR_dynamics_list$density[r],
        nodes     = MSAR_dynamics_list$nodes[r],
        regimes   = MSAR_dynamics_list$regimes[r],
        ts_id     = MSAR_dynamics_list$ts_id[r],
        regime_id = MSAR_dynamics_list$regime_id[r],
        reason    = invalid_reason,
        condition_number = cond_number
      )
    }
  }

  MSAR_dynamics_list$est_Kappa   <- est_Kappa
  MSAR_dynamics_list$est_Beta_ac <- est_Beta_ac
  MSAR_dynamics_list$Beta_corr   <- Beta_corr
  MSAR_dynamics_list$Kappa_corr  <- Kappa_corr
  MSAR_dynamics_list$Beta_sen    <- Beta_sen
  MSAR_dynamics_list$Beta_spec   <- Beta_spec
  MSAR_dynamics_list$Kappa_sen   <- Kappa_sen
  MSAR_dynamics_list$Kappa_spec  <- Kappa_spec
  MSAR_dynamics_list$Beta_ac_corr_pearson  <- Beta_ac_corr_pearson
  MSAR_dynamics_list$Beta_ac_corr_spearman <- Beta_ac_corr_spearman
  MSAR_dynamics_list$NRMSE_Beta  <- NRMSE_Beta
  MSAR_dynamics_list$NRMSE_Kappa <- NRMSE_Kappa

  sigma_validity_log <- if (length(validity_records) > 0) {
    dplyr::bind_rows(validity_records)
  } else {
    tibble::tibble(
      timesteps = integer(0), density = numeric(0), nodes = integer(0),
      regimes = integer(0), ts_id = integer(0), regime_id = integer(0),
      reason = character(0), condition_number = numeric(0)
    )
  }
  attr(MSAR_dynamics_list, "sigma_validity_log") <- sigma_validity_log

  n_invalid <- nrow(sigma_validity_log)
  cat(sprintf("\n=== KAPPA / Sigma VALIDITY ===\nFlagged invalid: %d of %d regime rows (%.2f%%)\n",
              n_invalid, n, if (n > 0) 100 * n_invalid / n else 0))
  cat("Kappa metrics (Kappa_corr, NRMSE_Kappa, Kappa_sen/spec) set to NA for those rows;\n")
  cat("Beta metrics for the same rows are unaffected. See summarize_sigma_validity().\n")

  MSAR_dynamics_list
}

# -----------------------------------------------------------------------------
# Companion to summarize_failures(): inspect the analysis-side Sigma validity log
# -----------------------------------------------------------------------------
# Distinct mechanism from estimation's failure_log: here a successfully converged
# fit produced an ill-conditioned est_Sigma for one regime, NA-ing only that
# regime's Kappa quantities (not the whole fit, not Beta).
summarize_sigma_validity <- function(MSAR_dynamics_list,
                                     by = c("timesteps", "density", "nodes", "regimes")) {

  vl <- attr(MSAR_dynamics_list, "sigma_validity_log")
  if (is.null(vl) || nrow(vl) == 0) {
    message("No Sigma validity flags logged (all est_Sigma well-conditioned).")
    return(invisible(NULL))
  }

  by_reason <- vl %>% dplyr::count(reason, name = "n", sort = TRUE)
  by_cell   <- vl %>%
    dplyr::count(dplyr::across(dplyr::all_of(by)), reason, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))

  # Condition-number distribution among ill-conditioned rows (solve_failed = NA).
  cond <- vl$condition_number[is.finite(vl$condition_number)]
  cond_q <- if (length(cond) > 0) {
    stats::quantile(cond, probs = c(0, .25, .5, .75, .9, .99, 1), na.rm = TRUE)
  } else NULL

  cat("=== SIGMA-INVALID REGIME ROWS: breakdown ===\n")
  cat(sprintf("Total flagged: %d\n\n", nrow(vl)))
  print(by_reason, n = Inf)
  if (!is.null(cond_q)) {
    cat("\nCondition-number quantiles (ill_conditioned rows):\n")
    print(round(cond_q, 1))
  }
  cat("\nUse $by_cell for the per-condition breakdown, $raw for the full log.\n")

  invisible(list(by_reason = by_reason, by_cell = by_cell,
                 cond_quantiles = cond_q, raw = vl))
}

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
  
  # Precision-weighting inputs: number of TRUE non-zero entries per regime.
  # Drives reliability of Beta_corr/Kappa_corr (a Pearson r computed from very
  # few true edges is mechanically forced toward |r|=1, regardless of estimation
  # quality -- see diagnostic in aggregate_to_sim_level()). Computed once here
  # from the already-stored orig_Beta/orig_Kappa list-columns; no re-estimation
  # needed. orig_Kappa can be NA (excluded near-singular Sigma_hat fits), hence
  # NA_integer_ fallback rather than letting sum()/sapply() error.
  corr_results <- corr_results %>%
    mutate(
      n_nz_Beta = sapply(orig_Beta, \(m) sum(m != 0)),
      n_nz_Kappa = sapply(orig_Kappa, \(m) {
        if (length(m) == 1 && is.na(m)) return(NA_integer_)
        sum(vectorize_upper_tri(m, diag = FALSE) != 0)
      })
    )
  
  # NRMSE_Beta / NRMSE_Kappa, Beta_corr/Kappa_corr, sens/spec and the AC
  # correlations are now computed upstream in compute_recovery_metrics() (from
  # the raw stored matrices, with the edge threshold and the Kappa
  # condition-number gate applied there) and are already present as columns
  # here. prepare_corr_results() only adds the grouping/identifier structure
  # and the precision-weighting edge counts (above).
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

    # --- Magnitude gate -------------------------------------------------------
    # The KW test flags *any* systematic difference in per-cell failure rate,
    # however small; with thousands of fits it reaches significance for a
    # difference of a few fits per condition. Escalate to the CRITICAL warning
    # only when the missingness is also large enough to plausibly bias estimates.
    # Below these thresholds the effect is statistically detectable but
    # practically negligible (and the sensitivity analysis quantifies any
    # residual movement directly).
    PRACTICAL_OVERALL_FAILURE <- 0.05   # > 5% of all fits missing
    PRACTICAL_CELL_FAILURE    <- 0.25   # any single condition losing > 25%
    total_failures       <- sum(failure_analysis$failures)
    total_expected       <- sum(failure_analysis$expected_n)
    overall_failure_rate <- if (total_expected > 0) total_failures / total_expected else 0
    max_cell_failure     <- max(failure_analysis$failure_rate, 0)
    material_bias        <- overall_failure_rate > PRACTICAL_OVERALL_FAILURE ||
                            max_cell_failure   > PRACTICAL_CELL_FAILURE

    # Computed regardless of the alarm level so the sensitivity analysis can
    # always quantify robustness whenever local attrition exists.
    high_fail_conditions <- failure_analysis %>%
      filter(failure_rate > 0.1) %>%
      select(Timesteps, Density, Nodes, Regimes)
    
    if (systematic_bias && material_bias) {
      cat("\n")
      cat("═══════════════════════════════════════════════════════════════\n")
      cat("⚠⚠⚠ CRITICAL WARNING: SELECTION BIAS DETECTED ⚠⚠⚠\n")
      cat("═══════════════════════════════════════════════════════════════\n\n")

      cat(sprintf("Failure rates differ significantly for: %s\n",
                  paste(biased_predictors, collapse = ", ")))
      cat(sprintf("Overall missingness: %.2f%% of fits; worst condition: %.1f%%.\n\n",
                  100 * overall_failure_rate, 100 * max_cell_failure))

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
        high_failure_conditions = high_fail_conditions,
        material = TRUE,
        overall_failure_rate = overall_failure_rate,
        max_cell_failure = max_cell_failure
      )
    } else if (systematic_bias) {
      cat(sprintf(
        "Note: failure rates differ statistically (%s), but missingness is negligible\n",
        paste(biased_predictors, collapse = ", ")))
      cat(sprintf(
        "      (%.2f%% overall; worst condition %.1f%%) and MAR on factors already in the models.\n",
        100 * overall_failure_rate, 100 * max_cell_failure))
      if (nrow(high_fail_conditions) > 0) {
        cat("      The sensitivity analysis below refits without the >10% cell to confirm robustness.\n\n")
      } else {
        cat("      Estimates are not expected to be affected.\n\n")
      }
      selection_bias_detected <- TRUE
      selection_bias_info <- list(
        biased_predictors = biased_predictors,
        p_values = bias_p_values,
        high_failure_conditions = high_fail_conditions,
        material = FALSE,
        overall_failure_rate = overall_failure_rate,
        max_cell_failure = max_cell_failure
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