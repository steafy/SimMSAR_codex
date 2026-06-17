# =============================================================================
# Sensitivity analysis: exclude high-failure conditions (>10%)
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PART 5.5).

# Returns the sens_results list (or NULL if no analysis needed).
run_sensitivity_analysis <- function(dat_sim, all_results, corr_cols, bias) {

  selection_bias_detected <- bias$selection_bias_detected
  high_fail_conditions    <- bias$high_fail_conditions

  if (selection_bias_detected &&
      !is.null(high_fail_conditions) && nrow(high_fail_conditions) > 0) {

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
    sens_results

  } else {
    cat("\nNo sensitivity analysis needed (no high-failure conditions or no bias detected).\n\n")
    NULL
  }
}
