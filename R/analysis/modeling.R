# =============================================================================
# Mixed-effects model fitting (main, 2-way, 3-way) with adaptive RE structure
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PART 5).
#
# RANDOM EFFECTS STRUCTURE:
# For each outcome the random-effects structure is selected adaptively:
# 1. Tests whether random slopes for logT improve fit (LR test)
# 2. If yes, checks for boundary singularity in correlated RE model
# 3. If singular (rho ~ 1.0), uses uncorrelated RE to ensure stability
# 4. Otherwise, uses correlated RE (best fit)

# -----------------------------------------------------------------------------
# Precision-weighting lookup (single source of truth for modeling.R AND
# sensitivity.R): which outcome uses which sim-level weight column from
# aggregate_to_sim_level(). Kappa_corr is weighted by true non-zero edge count
# (n_nz); Beta_ac_corr_pearson by Nodes (the analogous "n" for a correlation
# computed across nodes rather than edges) -- both correlations are mechanically
# pulled toward |r|=1 at low n, regardless of estimation quality. Beta_corr is
# NOT weighted (see transform.R: the full N x N matrix gives it enough degrees
# of freedom that this artifact doesn't bite). Outcomes not listed here are fit
# unweighted.
# -----------------------------------------------------------------------------
OUTCOME_WEIGHT_MAP <- c(Kappa_corr = "w_Kappa", Beta_ac_corr_pearson = "w_BetaAC")

# -----------------------------------------------------------------------------
# Helpers (remove the triplicated formula construction / fitting boilerplate)
# -----------------------------------------------------------------------------

# term: "main" (main effects), "2way" (^2), "3way" (^3)
build_formula <- function(outcome_z, term, re_formula) {
  fixed <- switch(
    term,
    main   = paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes"),
    `2way` = paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2"),
    `3way` = paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^3")
  )
  as.formula(paste0(fixed, " + ", re_formula))
}

fit_lmer <- function(formula, data, weights = NULL) {
  lmer(formula, data = data, weights = weights, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa"))
}

# -----------------------------------------------------------------------------
# Fit all three models for a single outcome and return the results entry.
# -----------------------------------------------------------------------------
fit_outcome_models <- function(outcome, dat_sim, weight_col = NULL) {

  outcome_z <- paste0(outcome, "_z")

  # Precision weights (see aggregate_to_sim_level()): NULL for outcomes
  # without an analogous true-edge count (e.g. Beta_ac_corr), in which case
  # fit_lmer()'s weights=NULL default reproduces the original unweighted fit.
  # Normalized to mean 1: lmer's fixed-effect estimates are invariant to a
  # global rescaling of weights, but the residual-variance scale (and hence
  # ICC, R2_nakagawa) is NOT -- raw w_Beta/w_Kappa run into the hundreds for
  # some rows, which inflated residual variance ~30x and likely destabilized
  # the optimizer (observed: "negative eigenvalue" convergence warning).
  w <- if (!is.null(weight_col)) dat_sim[[weight_col]] / mean(dat_sim[[weight_col]]) else NULL

  cat("Fitting models for:", outcome, "\n")
  if (!is.null(weight_col)) {
    cat("  (precision-weighted by", weight_col, ")\n")
  }
  cat(strrep("-", 60), "\n")

  # 5.0: Test random slopes for logT
  # Different networks may have different trajectories across T
  cat("Testing random-effects structure...\n")

  # Model with random intercept only
  model_ri <- fit_lmer(build_formula(outcome_z, "main", "(1|network_id)"), dat_sim, weights = w)

  # Model with random intercept + random slope for logT (correlated)
  model_rs_corr <- tryCatch({
    fit_lmer(build_formula(outcome_z, "main", "(1 + logT|network_id)"), dat_sim, weights = w)
  }, error = function(e) {
    cat("  Random slopes (correlated) model failed to converge\n")
    return(NULL)
  })

  # Model with random intercept + random slope for logT (uncorrelated)
  model_rs_uncorr <- tryCatch({
    fit_lmer(build_formula(outcome_z, "main",
                           "(1|network_id) + (0 + logT|network_id)"), dat_sim, weights = w)
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
  model_main <- fit_lmer(build_formula(outcome_z, "main", re_formula), dat_sim, weights = w)

  # 5.2: 2-way interactions (PRIMARY MODEL)
  model_2way <- fit_lmer(build_formula(outcome_z, "2way", re_formula), dat_sim, weights = w)

  # 5.2.1: CHECK FOR SINGULARITY IN 2-WAY MODEL AND REFIT IF NEEDED
  # The 2-way model may be singular even if the main effects test was not
  if (use_random_slopes && !use_uncorrelated && isSingular(model_2way)) {
    cat("\n⚠ SINGULARITY DETECTED in 2-way model!\n")
    cat("  The 2-way model is singular even though main effects test was not.\n")
    cat("  Refitting all models with uncorrelated random effects...\n\n")

    # Update flags and formula
    use_uncorrelated <- TRUE
    re_formula <- "(1|network_id) + (0 + logT|network_id)"

    # Refit main effects + 2-way models
    model_main <- fit_lmer(build_formula(outcome_z, "main", re_formula), dat_sim, weights = w)
    model_2way <- fit_lmer(build_formula(outcome_z, "2way", re_formula), dat_sim, weights = w)

    cat("  ✓ Models refitted with uncorrelated random effects\n")
    cat(sprintf("  New formula: %s\n\n", re_formula))
  }

  # 5.3: 3-way interactions
  model_3way <- fit_lmer(build_formula(outcome_z, "3way", re_formula), dat_sim, weights = w)

  # 5.3.1: CHECK FOR SINGULARITY IN 3-WAY MODEL AND REFIT IF NEEDED
  if (use_random_slopes && !use_uncorrelated && isSingular(model_3way)) {
    cat("\n⚠ SINGULARITY DETECTED in 3-way model!\n")
    cat("  Refitting all models with uncorrelated random effects...\n\n")

    use_uncorrelated <- TRUE
    re_formula <- "(1|network_id) + (0 + logT|network_id)"

    model_main <- fit_lmer(build_formula(outcome_z, "main", re_formula), dat_sim, weights = w)
    model_2way <- fit_lmer(build_formula(outcome_z, "2way", re_formula), dat_sim, weights = w)
    model_3way <- fit_lmer(build_formula(outcome_z, "3way", re_formula), dat_sim, weights = w)

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
    cat(sprintf("\n→ 2-Way model selected as primary (ΔAIC[3way-2way] = %.1f; 3-way improves AIC by less than the 10-point threshold, so the simpler model is kept)\n",
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
  list(
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
# Fit models for every outcome; returns the named all_results list.
# -----------------------------------------------------------------------------
fit_all_models <- function(dat_sim, corr_cols) {

  cat("\n=== Fitting Mixed-Effects Models ===\n")
  cat("Note: Testing random slopes for logT to determine appropriate random-effects structure\n\n")

  all_results <- list()
  for (outcome in corr_cols) {
    weight_col <- if (outcome %in% names(OUTCOME_WEIGHT_MAP)) OUTCOME_WEIGHT_MAP[[outcome]] else NULL
    all_results[[outcome]] <- fit_outcome_models(outcome, dat_sim, weight_col = weight_col)
  }
  all_results
}
