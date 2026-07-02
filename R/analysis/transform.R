# =============================================================================
# Fisher-z transformation, simulation-level aggregation, predictor scaling
# =============================================================================
# Functions extracted from the former monolithic stat_analysis.R (PARTs 2-4).

# -----------------------------------------------------------------------------
# Shared helper: eps-clipped Fisher-z (atanh), used by both add_fisher_z() and
# aggregate_to_sim_level() so the clipping threshold has one source of truth.
# -----------------------------------------------------------------------------
safe_atanh <- function(r, eps = 1e-6) {
  atanh(pmin(pmax(r, -1 + eps), 1 - eps))
}

# -----------------------------------------------------------------------------
# PART 2: Fisher-z transformation of correlation outcomes
# -----------------------------------------------------------------------------
# Adds <outcome>_z columns to corr_results and warns about extreme z-values.
add_fisher_z <- function(corr_results, corr_cols) {

  for (cc in corr_cols) {
    # Fisher-z transformation: z = atanh(r), clipped to avoid Inf/-Inf
    corr_results[[paste0(cc, "_z")]] <- safe_atanh(corr_results[[cc]])
  }
  
  # Check for extreme values
  z_cols <- paste0(corr_cols, "_z")
  for (zc in z_cols) {
    n_extreme <- sum(abs(corr_results[[zc]]) > 5, na.rm = TRUE)
    n_na      <- sum(is.na(corr_results[[zc]]))
    if (n_extreme > 0) {
      cat("  ⚠", zc, "has", n_extreme, "extreme values (|z| > 5)\n")
    }
    if (n_na > 0) {
      cat("  ⚠", zc, "has", n_na, "NA values\n")
    }
  }
  
  corr_results
}

# -----------------------------------------------------------------------------
# PART 3: Aggregate to simulation level
# -----------------------------------------------------------------------------
# Correlations aggregated on the z-scale; other metrics on the original scale.
# Returns list(dat_sim, available_other_metrics, z_cols).
aggregate_to_sim_level <- function(corr_results, corr_cols) {
  
  z_cols <- paste0(corr_cols, "_z")
  
  # Precision-weighting: outcomes computed from very few TRUE non-zero edges
  # (n_nz_Kappa, added in prepare_corr_results()) have a Pearson r that is
  # mechanically forced toward |r|=1, regardless of estimation quality.
  # Weight = pmax(n_nz - 3, 1), the standard Fisher-z sampling-variance weight
  # (Var(z) ~ 1/(n-3)), floored at 1 so n_nz in {2,3,4} doesn't go <=0.
  # Alternative if more differentiation at small n_nz is wanted: weight = n_nz.
  #
  # Beta_ac_corr_pearson has the analogous artifact, but with a different
  # "n": it correlates the AC vector across Nodes (4/6/8), not across edges, so
  # its precision weight uses Nodes instead of n_nz. Decided 2026-06 after the
  # pilot run showed 1.5% |z|>5 for Pearson AC (vs. 18.1% for the now-dropped
  # Spearman variant -- discreteness at N=4 makes rank correlation collapse
  # onto a handful of values regardless of estimation quality).
  #
  # Beta_corr is deliberately NOT weighted (was, until 2026-06-29; see git
  # history). This now matches modeling.R's OUTCOME_WEIGHT_MAP, which already
  # excluded Beta_corr from LMM-level weighting on the documented grounds that
  # the full N x N matrix gives it enough degrees of freedom -- but the
  # regime-level aggregation here used to weight it anyway, which was an
  # unresolved code/text inconsistency. Checked empirically (regime-level mean
  # |Beta_corr_z| by Density x Nodes cell): the lowest-n_nz cell (Density=0.25,
  # Nodes=4, n_nz~4) had the LOWEST mean |z| in the design, not the highest --
  # the opposite of what the "few true edges forces r toward 1" mechanism
  # predicts. What the cross-tabulation actually shows is a real (if small,
  # see the I(Density_s^2) check in modeling.R) Density recovery effect, which
  # a precision weight would suppress/misrepresent rather than correctly
  # reflect. Both levels (LMM and regime-aggregation) are now unweighted for
  # Beta_corr, not just one.
  weighted_outcomes <- c(Kappa_corr = "weight_Kappa", Beta_ac_corr_pearson = "weight_BetaAC")
  weighted_z_cols   <- paste0(names(weighted_outcomes), "_z")
  
  # Convert factors to numeric for aggregation
  corr_results_for_agg <- corr_results %>%
    mutate(
      Timesteps_num = as.numeric(as.character(Timesteps)),
      Density_num = as.numeric(as.character(Density)),
      Nodes_num = as.numeric(as.character(Nodes)),
      Regimes_num = as.numeric(as.character(Regimes)),
      weight_Kappa  = pmax(n_nz_Kappa - 3, 1),
      weight_BetaAC = pmax(Nodes_num - 3, 1)
    )
  
  # Aggregate: Mean across RegimeIndex for each SimUID
  # CORRELATIONS: aggregate ON THE RAW r-SCALE (precision-weighted for
  # Beta_corr/Kappa_corr/Beta_ac_corr_pearson), THEN Fisher-z-transform the
  # resulting per-simulation value ONCE. Any remaining z-outcome not in
  # weighted_outcomes keeps a simple mean of the already-computed _z column
  # (currently none are routed that way; kept for outcomes added later that
  # are deliberately excluded from the weighting, e.g. Beta_ac_corr_spearman).
  #
  # CHANGED 2026-06-28 (was: weighted.mean on the _z columns directly, i.e.
  # transform-then-average). atanh(r) -> Inf as r -> 1, so with only 1-4
  # regimes per simulation, a single near-perfect regime could drag the
  # back-transformed simulation-level mean far above where the regimes
  # actually sit on average -- verified empirically: the old approach showed
  # a SYSTEMATIC upward bias that grows with regime count (M=1: 0; M=2:
  # +0.026; M=3: +0.062; M=4: +0.083 on the r-scale, averaged over all
  # affected simulations), i.e. it specifically understated the recovery
  # decline for high M, the design's dominant effect. Averaging on the raw
  # scale first removes the mechanism entirely (no value can leverage the
  # mean toward +-Inf), and the eps-clipped atanh() applied once afterwards
  # still gives the LMMs the same bounded, variance-stabilised outcome scale.
  plain_z_cols <- setdiff(z_cols, weighted_z_cols)
  
  dat_sim <- corr_results_for_agg %>%
    group_by(SimUID, Condition, SimID,
             Timesteps_num, Density_num, Nodes_num, Regimes_num) %>%
    summarise(
      # Beta_corr: unweighted mean on the raw scale, then atanh once -- same
      # "average-raw-first" principle as the weighted outcomes below (see the
      # 2026-06-28 fix note above this block), just without a precision
      # weight, per the rationale in the comment above weighted_outcomes.
      Beta_corr_z            = safe_atanh(mean(Beta_corr, na.rm = TRUE)),
      Kappa_corr_z           = safe_atanh(stats::weighted.mean(Kappa_corr,          w = weight_Kappa,  na.rm = TRUE)),
      Beta_ac_corr_pearson_z = safe_atanh(stats::weighted.mean(Beta_ac_corr_pearson, w = weight_BetaAC, na.rm = TRUE)),
      across(all_of(plain_z_cols), \(x) mean(x, na.rm = TRUE), .names = "{.col}"),
      w_Kappa = sum(weight_Kappa[!is.na(Kappa_corr)], na.rm = TRUE),
      w_BetaAC = sum(weight_BetaAC[!is.na(Beta_ac_corr_pearson)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(
      Timesteps = Timesteps_num,
      Density = Density_num,
      Nodes = Nodes_num,
      Regimes_fac = Regimes_num
    )
  
  # OTHER METRICS: Aggregate on original scale
  other_metrics <- c("NRMSE_Beta", "NRMSE_Kappa",
                     "Beta_sen", "Kappa_sen",
                     "Beta_spec", "Kappa_spec")
  
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
  
  list(
    dat_sim                 = dat_sim,
    available_other_metrics = available_other_metrics,
    z_cols                  = z_cols
  )
}

# -----------------------------------------------------------------------------
# PART 4: Transform & scale predictors; attach condition-level success rates
# -----------------------------------------------------------------------------
scale_predictors <- function(dat_sim, selection_weights) {
  
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
  
  dat_sim
}