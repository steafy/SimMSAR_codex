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
      # Observation-level random effect (OLRE): one unique level per fit. Added
      # as a grouping factor here so the secondary binomial GLMM can carry a
      # (1|obs_id) term, which absorbs the heavy extra-binomial overdispersion
      # (length-weighted accuracy is far more variable across fits than a pure
      # binomial allows). See fit_sequence_recovery_model() for the rationale.
      obs_id    = factor(dplyr::row_number()),
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

  # --- Interaction check: does logT x Regimes matter for RQ4? ----------------
  # RQ1--3 report a logT x Regimes interaction as a central finding ("T helps
  # recovery more strongly at higher M"). RQ4 deliberately uses a main-effects-
  # only LMM (separately justified -- simpler model, Cohen's kappa as a more
  # holistic outcome than the per-parameter correlations), but a reviewer can
  # reasonably ask whether that interaction was checked here too rather than
  # silently skipped because it didn't fit the simpler-model narrative. This
  # fits a single additional interaction -- not the full 3-way structure used
  # for RQ1--3 -- and applies the SAME decision rule already used throughout
  # modeling.R (Delta AIC > 10 favours the more complex model), so RQ4 model
  # selection stays methodologically consistent with RQ1--3.
  m_interact <- lmer(
    cohens_kappa ~ logT * Regimes + Density_s + Nodes_s + (1 | param_set_id),
    data = seq_results, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  delta_aic_interact <- AIC(m_primary) - AIC(m_interact)  # positive favours interaction model
  used_interaction    <- delta_aic_interact > 10           # single source of truth, reused below

  cat(sprintf("\nlogT x Regimes interaction check: AIC(main effects) = %.1f, AIC(+ logT x Regimes) = %.1f\n",
              AIC(m_primary), AIC(m_interact)))
  cat(sprintf("Delta AIC (main effects vs. + interaction) = %.1f\n", delta_aic_interact))

  if (used_interaction) {
    m_rq4_primary      <- m_interact
    rq4_primary_label  <- "Main effects + logT x Regimes (primary)"
    cat("-> logT x Regimes interaction selected as primary RQ4 model (Delta AIC > 10).\n")
    cat("\n--- Primary LMM (with interaction): cohens_kappa ~ logT * Regimes + Density_s + Nodes_s + (1|param_set_id) ---\n")
    print(summary(m_interact)$coefficients)
    cat(sprintf("\nMarginal R2: %.3f   Conditional R2: %.3f\n",
                performance::r2_nakagawa(m_interact)$R2_marginal,
                performance::r2_nakagawa(m_interact)$R2_conditional))
  } else {
    m_rq4_primary      <- m_primary
    rq4_primary_label  <- "Main effects (primary)"
    cat("-> No evidence for a logT x Regimes interaction (Delta AIC <= 10); main-effects model remains primary.\n")
    cat(sprintf("   (Note: with n = %d observations, the interaction's own p-value would almost certainly\n", nrow(seq_results)))
    cat("    be significant regardless of practical relevance -- Delta AIC, not significance, is the\n")
    cat("    deciding criterion here, for consistency with the RQ1-3 model-selection rule.)\n")
  }

  # --- Secondary: binomial GLMM on accuracy (length-weighted) -----------------
  # An ordinary binomial GLMM on per-fit accuracy is massively overdispersed
  # (Pearson/df ~ 60): sequence-recovery accuracy varies across fits far more
  # than a binomial with these large n_total values permits, so the binomial SEs
  # are anti-conservative. Rather than switch outcome family, we add an
  # observation-level random effect (1|obs_id) -- one latent normal deviate per
  # fit -- which is the standard lme4 remedy for binomial overdispersion: it
  # soaks up the extra-binomial variance into a quantified random term, leaving
  # the fixed-effect inference (the robustness check we actually care about)
  # valid. The conditional-Pearson overdispersion ratio is recomputed below and
  # should now sit near 1.
  #
  # Formula mirrors whichever model won the interaction check above
  # (used_interaction), so the robustness check always tests the same spec
  # that ends up reported as primary -- a GLMM without logT x Regimes would
  # not actually be a robustness check FOR the interaction finding when that
  # interaction is what's being reported.
  secondary_formula <- if (used_interaction) {
    cbind(n_correct, n_incorrect) ~ logT * Regimes + Density_s + Nodes_s +
      (1 | param_set_id) + (1 | obs_id)
  } else {
    cbind(n_correct, n_incorrect) ~ logT + Density_s + Nodes_s + Regimes +
      (1 | param_set_id) + (1 | obs_id)
  }
  m_secondary <- tryCatch(
    glmer(
      secondary_formula,
      data = seq_results, family = binomial,
      control = glmerControl(optimizer = "bobyqa")
    ),
    error = function(e) { cat("  Binomial GLMM failed:", e$message, "\n"); NULL }
  )
  if (!is.null(m_secondary)) {
    cat(sprintf("\n--- Secondary binomial GLMM on accuracy + OLRE (robustness check%s) ---\n",
                if (used_interaction) ", + logT x Regimes" else ""))
    print(summary(m_secondary)$coefficients)
    # Overdispersion check (Pearson chi-square / residual df). With the OLRE in
    # the model this is the CONDITIONAL ratio: the per-observation random effect
    # has absorbed the extra-binomial variance, so a value near 1 indicates the
    # remaining residual dispersion is well-behaved.
    rdf <- df.residual(m_secondary)
    od  <- sum(residuals(m_secondary, type = "pearson")^2) / rdf
    cat(sprintf("Overdispersion ratio (with OLRE): %.2f", od))
    if (od > 1.5) {
      cat("  (still >1.5: residual overdispersion remains beyond the OLRE)\n")
    } else {
      cat("  (acceptable -- OLRE absorbed the overdispersion)\n")
    }
  }

  invisible(list(
    seq_results        = seq_results,
    empty_model        = m0,
    icc                = icc0,
    # main_effects_model: always the simple logT + Density_s + Nodes_s + Regimes
    # spec. Kept separate from primary_model so the secondary GLMM and the
    # sensitivity refit (run_sequence_sensitivity_analysis) can keep comparing
    # against this exact spec regardless of which model wins the interaction
    # check below -- both of those were built around an interaction-free
    # coefficient set, and extending them to a logT x Regimes term is a
    # separate decision, not an automatic consequence of this check.
    main_effects_model = m_primary,
    interaction_model   = m_interact,
    delta_aic_interact  = delta_aic_interact,
    used_interaction    = used_interaction,  # single source of truth for downstream consumers
    # primary_model/primary_label: whichever of the two above won the Delta
    # AIC > 10 comparison -- this is what export_sequence_recovery_table()
    # writes out as "the" RQ4 coefficient table, and what the sensitivity
    # refit below mirrors via used_interaction.
    primary_model = m_rq4_primary,
    primary_label = rq4_primary_label,
    secondary_model = m_secondary
  ))
}

# =============================================================================
# RQ4 sensitivity analysis: exclude high-missingness conditions (>10%)
# =============================================================================
# Mirror of run_sensitivity_analysis() (sensitivity.R) for the regime-sequence
# recovery LMM. The question is identical -- do MNAR estimation failures bias the
# reported recovery estimates? -- but the relevant missingness differs from the
# correlation outcomes', so it is recomputed here rather than reusing
# bias$high_fail_conditions (which is built from the corr-outcome failures only):
#
#   * For Beta_corr/Kappa_corr/Beta_ac the unit is a per-REGIME row (regime_rows);
#     a fit is "missing" if estimation failed (fit_null / zero_var_beta / no_data).
#   * For RQ4 the unit is a per-FIT row (fit_results), absent under ANY of those
#     same estimation failures, under seq_length_mismatch (the EM fit + regime
#     matching succeeded, but the true and reconstructed regime sequences had
#     different lengths, so no fit-level sequence row could be built), AND under
#     partial_zero_var_beta (SOME but not all estimated Betas were degenerate: the
#     healthy regimes are salvaged into the main table for RQ1-3, but no fit-level
#     sequence row is built, so the fit is deliberately kept OUT of RQ4 -- a
#     partial regime map cannot support a well-defined full-sequence kappa). Any
#     of these removes a replicate from the RQ4 sample, so the per-cell RQ4
#     missingness is the COMBINED set. Because the five mechanisms are mutually
#     exclusive per replicate (fit_one_replicate returns at the first one it
#     hits), a missing fit_results row corresponds to exactly one of them -- so
#     the combined missingness is measured directly as (n_ts - number of
#     fit_results rows) per design cell, with no need to parse the failure log to
#     DEFINE it. The log is used only to ATTRIBUTE that missingness by stage
#     (printed below), which in particular makes the seq_length_mismatch and
#     partial_zero_var_beta contributions explicit.
#
# Exclusion threshold (>10% per cell, fail_threshold) and the Full-vs-Sensitivity
# coefficient comparison match sensitivity.R so the two slot together in the
# pipeline's reporting.
run_sequence_sensitivity_analysis <- function(seq_recovery, fit_results,
                                              failure_log, n_ts_per_condition,
                                              fail_threshold = 0.10) {

  cat("\n=== RQ4 SENSITIVITY ANALYSIS: regime-sequence recovery ===\n\n")

  if (is.null(seq_recovery) || is.null(seq_recovery$seq_results) ||
      nrow(seq_recovery$seq_results) == 0) {
    cat("No RQ4 model available (no regimes >= 2 fits); skipping.\n\n")
    return(NULL)
  }

  seq_results   <- seq_recovery$seq_results
  # primary_model here is whichever model fit_sequence_recovery_model() reported
  # as primary (main effects, or + logT x Regimes if that won Delta AIC > 10).
  # The refit below mirrors the SAME spec via used_interaction, so Full vs.
  # Sensitivity stays an apples-to-apples comparison of the model that's
  # actually reported -- not a comparison against a simpler model nobody is
  # claiming as the result.
  primary_model     <- seq_recovery$primary_model
  used_interaction  <- isTRUE(seq_recovery$used_interaction)

  # --- Per-cell combined RQ4 missingness -------------------------------------
  # Build the FULL M>=2 design grid first, so a cell where every replicate
  # failed (zero fit_results rows) still appears as 100% missing rather than
  # being silently absent from a plain count().
  reg_levels <- sort(unique(fit_results$regimes[fit_results$regimes >= 2]))
  grid <- expand.grid(
    timesteps = sort(unique(fit_results$timesteps)),
    density   = sort(unique(fit_results$density)),
    nodes     = sort(unique(fit_results$nodes)),
    regimes   = reg_levels,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  cell_counts <- fit_results %>%
    filter(regimes >= 2) %>%
    count(timesteps, density, nodes, regimes, name = "n_fits")

  rq4_failures <- grid %>%
    left_join(cell_counts, by = c("timesteps", "density", "nodes", "regimes")) %>%
    mutate(
      n_fits       = dplyr::coalesce(n_fits, 0L),
      expected_n   = n_ts_per_condition,
      failures     = expected_n - n_fits,
      failure_rate = failures / expected_n
    ) %>%
    arrange(desc(failure_rate))

  # --- Attribute the combined missingness by stage (transparency) ------------
  # Confirms WHICH mechanisms drive the RQ4 attrition and, in particular, makes
  # the seq_length_mismatch contribution explicit (it has been 0 in practice --
  # the fit-level sequence guard in estimate_MSAR() never fired across the runs).
  if (!is.null(failure_log) && nrow(failure_log) > 0) {
    fl_rq4    <- failure_log %>% filter(regimes >= 2)
    stage_tab <- sort(table(fl_rq4$stage), decreasing = TRUE)
    cat(sprintf("Combined RQ4 missingness: %d fits removed of %d expected (%.2f%%).\n",
                sum(rq4_failures$failures), sum(rq4_failures$expected_n),
                100 * sum(rq4_failures$failures) / sum(rq4_failures$expected_n)))
    if (length(stage_tab) > 0) {
      cat("  By stage: ",
          paste(names(stage_tab), as.integer(stage_tab), sep = "=", collapse = ", "),
          "\n", sep = "")
    }
    if (!"seq_length_mismatch" %in% names(stage_tab)) {
      cat("  (seq_length_mismatch = 0: the fit-level sequence guard never fired.)\n")
    }
    cat("\n")
  }

  high_fail_conditions <- rq4_failures %>%
    filter(failure_rate > fail_threshold) %>%
    select(timesteps, density, nodes, regimes, failures, failure_rate)

  if (nrow(high_fail_conditions) == 0) {
    cat(sprintf(paste0("No condition exceeded the %.0f%% combined-missingness ",
                       "threshold; no RQ4 sensitivity refit needed.\n\n"),
                100 * fail_threshold))
    return(invisible(list(rq4_failures = rq4_failures,
                          high_fail_conditions = high_fail_conditions)))
  }

  cat(sprintf("Excluded conditions (>%.0f%% combined missingness):\n",
              100 * fail_threshold))
  print(as.data.frame(high_fail_conditions), row.names = FALSE)
  cat("\n")

  # --- Exclude high-missingness cells and refit the primary RQ4 LMM ----------
  # Reuse the already-scaled predictors in seq_results (do NOT rescale on the
  # reduced data), exactly as sensitivity.R does, so Full vs. Sensitivity stays
  # an apples-to-apples coefficient comparison. The fixed/random-effects spec
  # mirrors whichever model fit_sequence_recovery_model() selected as primary
  # (see used_interaction above).
  excl <- high_fail_conditions %>% select(timesteps, density, nodes, regimes)
  seq_results_sens <- seq_results %>%
    anti_join(excl, by = c("timesteps", "density", "nodes", "regimes"))

  cat(sprintf("Rows in full seq_results: %d\n", nrow(seq_results)))
  cat(sprintf("Rows after exclusion:     %d\n", nrow(seq_results_sens)))
  cat(sprintf("Rows removed:             %d\n\n",
              nrow(seq_results) - nrow(seq_results_sens)))

  sens_formula <- if (used_interaction) {
    cohens_kappa ~ logT * Regimes + Density_s + Nodes_s + (1 | param_set_id)
  } else {
    cohens_kappa ~ logT + Density_s + Nodes_s + Regimes + (1 | param_set_id)
  }
  model_sens <- tryCatch(
    lmer(sens_formula,
         data = seq_results_sens, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { cat("  RQ4 sensitivity model failed:", e$message, "\n"); NULL }
  )
  if (is.null(model_sens)) return(NULL)

  # Compare primary (Full) vs. sensitivity coefficients
  coef_full <- fixef(primary_model)
  coef_sens <- fixef(model_sens)
  common    <- intersect(names(coef_full), names(coef_sens))

  comparison <- data.frame(
    Predictor   = common,
    Full        = round(coef_full[common], 4),
    Sensitivity = round(coef_sens[common], 4),
    Delta       = round(coef_sens[common] - coef_full[common], 4),
    Delta_pct   = round(100 * (coef_sens[common] - coef_full[common]) /
                          ifelse(abs(coef_full[common]) < 1e-10, NA, coef_full[common]), 1),
    row.names   = NULL
  )

  cat("Outcome: cohens_kappa (RQ4 regime-sequence recovery)\n")
  cat(strrep("-", 70), "\n")
  print(comparison, row.names = FALSE)
  cat("\nDelta = Sensitivity - Full. Large |Delta| or |Delta_pct| indicate bias.\n\n")

  invisible(list(
    model                = model_sens,
    comparison           = comparison,
    high_fail_conditions = high_fail_conditions,
    rq4_failures         = rq4_failures
  ))
}
