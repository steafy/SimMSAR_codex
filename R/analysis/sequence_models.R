# =============================================================================
# RQ4: Regime-sequence recovery -- accuracy / Cohen's kappa, single-RI LMM
# =============================================================================
# Consumes the fit-level table attached by estimate_MSAR() as
# attr(MSAR_dynamics_list, "fit_results"): one row per fit with the (already
# regime-matched) true and estimated regime sequences. Cohen's kappa and
# accuracy are computed HERE (not in estimation) so the scoring stays revisable.
#
# Nesting: the generating parameter set (regime_dynamics + transmat) is reused
# across all levels of T within a (density, nodes, regimes, ts_id) cell, so the
# T-level observations within such a cell are correlated. The grouping factor is
#   param_set_id <- interaction(density, nodes, regimes, ts_id, drop = TRUE)
# Within a param_set_id only timesteps varies; density/nodes/regimes are constant.

source("R/utils/regime_sequence_recovery.R")   # cohens_kappa_manual()

# Per-fit accuracy / Cohen's kappa + counts, restricted to regimes >= 2.
# Predictors are z-scaled (logT, Density_s, Nodes_s) with Regimes as a factor,
# matching the Beta/Kappa LMMs so coefficients are comparable across outcomes.
prepare_sequence_results <- function(fit_results) {

  seq_df <- fit_results %>%
    filter(regimes >= 2)

  if (nrow(seq_df) == 0) {
    return(seq_df[0, , drop = FALSE])
  }

  seq_df %>%
    rowwise() %>%
    mutate(
      n_total     = length(orig_regime_sequence),
      n_correct   = sum(orig_regime_sequence == est_regime_sequence),
      n_incorrect = n_total - n_correct,
      accuracy    = n_correct / n_total,
      cohens_kappa = cohens_kappa_manual(orig_regime_sequence, est_regime_sequence)
    ) %>%
    ungroup() %>%
    mutate(
      param_set_id = interaction(density, nodes, regimes, ts_id, drop = TRUE),
      logT      = as.numeric(scale(log(timesteps))),
      Density_s = as.numeric(scale(density)),
      Nodes_s   = as.numeric(scale(nodes)),
      Regimes   = factor(regimes)
    )
}

# Fit the RQ4 models and report. Primary = Gaussian LMM on Cohen's kappa with a
# single random intercept for param_set_id; secondary = binomial GLMM on
# accuracy (length-weighted) as a robustness check. Also reports the empty-model
# ICC as the exploratory justification for the random intercept.
fit_sequence_recovery_model <- function(fit_results) {

  cat("\n=== RQ4: REGIME-SEQUENCE RECOVERY (Cohen's kappa LMM) ===\n")

  seq_results <- prepare_sequence_results(fit_results)

  if (nrow(seq_results) == 0) {
    cat("No multi-regime fits (regimes >= 2) available; skipping RQ4.\n")
    return(invisible(NULL))
  }

  n_clusters <- nlevels(seq_results$param_set_id)
  cat(sprintf("Fits: %d   param_set_id clusters: %d   (avg %.1f fits/cluster)\n",
              nrow(seq_results), n_clusters, nrow(seq_results) / n_clusters))

  # --- Empty model: param_set_id ICC ------------------------------------------
  m0 <- lmer(cohens_kappa ~ 1 + (1 | param_set_id), data = seq_results, REML = TRUE)
  vc <- as.data.frame(VarCorr(m0))
  var_cluster  <- vc$vcov[vc$grp == "param_set_id"]
  var_resid    <- vc$vcov[vc$grp == "Residual"]
  icc0 <- var_cluster / (var_cluster + var_resid)
  cat(sprintf("\nEmpty-model param_set_id ICC: %.3f (%.1f%% between-cluster)\n",
              icc0, 100 * icc0))
  if (icc0 < 0.01) {
    cat("  Note: ICC ~ 0 -- a plain lm would essentially suffice; the random intercept\n")
    cat("        is retained for principled clustering but adds little here.\n")
  }

  # --- Primary: Gaussian LMM on Cohen's kappa ---------------------------------
  m_primary <- lmer(
    cohens_kappa ~ logT + Density_s + Nodes_s + Regimes + (1 | param_set_id),
    data = seq_results, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  cat("\n--- Primary LMM: cohens_kappa ~ logT + Density_s + Nodes_s + Regimes + (1|param_set_id) ---\n")
  print(summary(m_primary)$coefficients)
  cat(sprintf("\nMarginal R2: %.3f   Conditional R2: %.3f\n",
              performance::r2_nakagawa(m_primary)$R2_marginal,
              performance::r2_nakagawa(m_primary)$R2_conditional))

  # --- Secondary: binomial GLMM on accuracy (length-weighted) -----------------
  m_secondary <- tryCatch(
    glmer(
      cbind(n_correct, n_incorrect) ~ logT + Density_s + Nodes_s + Regimes + (1 | param_set_id),
      data = seq_results, family = binomial,
      control = glmerControl(optimizer = "bobyqa")
    ),
    error = function(e) { cat("  Binomial GLMM failed:", e$message, "\n"); NULL }
  )
  if (!is.null(m_secondary)) {
    cat("\n--- Secondary binomial GLMM on accuracy (robustness check) ---\n")
    print(summary(m_secondary)$coefficients)
    # Quick overdispersion check (Pearson chi-square / residual df).
    rdf <- df.residual(m_secondary)
    od  <- sum(residuals(m_secondary, type = "pearson")^2) / rdf
    cat(sprintf("Overdispersion ratio: %.2f", od))
    if (od > 1.5) {
      cat("  (>1.5: consider an observation-level random effect rather than switching outcome)\n")
    } else {
      cat("  (acceptable)\n")
    }
  }

  invisible(list(
    seq_results   = seq_results,
    empty_model   = m0,
    icc           = icc0,
    primary_model = m_primary,
    secondary_model = m_secondary
  ))
}
