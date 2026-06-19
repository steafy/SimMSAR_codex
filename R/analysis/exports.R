# =============================================================================
# Export results: coefficient tables, model comparison table, diagnostics PDF
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PART 6).

# -----------------------------------------------------------------------------
# 6.1: Coefficient tables for the primary models (HTML + LaTeX)
# -----------------------------------------------------------------------------
export_coefficient_tables <- function(all_results, corr_cols) {

  cat("\n=== Exporting Results ===\n\n")

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
}

# -----------------------------------------------------------------------------
# 6.2: Model comparison table (one row per outcome) + LaTeX export
# -----------------------------------------------------------------------------
export_comparison_table <- function(all_results, corr_cols) {

  # LaTeX-formatted outcome labels (fully written out for readability)
  # AC = average controllability; corr. = Pearson correlation with true matrix
  # Recovery-correlation symbols (subscript style; requires amsmath for \text).
  # r_AC is the correlation of the (nodal) average controllability AC(Beta).
  outcome_tex_labels <- c(
    Beta_corr     = "$r_{\\text{Beta}}$",
    Kappa_corr    = "$r_{\\text{Kappa}}$",
    Beta_ac_corr  = "$r_{\\text{AC}}$"
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
      # Effective N actually entering the model (listwise-complete obs). With the
      # Kappa condition-number guard this is smaller for Kappa-based outcomes.
      N             = stats::nobs(res$primary_model),
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
    "Outcome", "Random Effects", "N", "AIC (Main)", "AIC (2-Way)", "AIC (3-Way)",
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
      row["N"],                " & ",
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
    "  \\begin{tabular}{@{}llrrrrcrrrr@{}}",
    "    \\toprule",
    "    Outcome & RE structure & $N$ &",
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

  comparison_all
}

# -----------------------------------------------------------------------------
# 6.3: Diagnostic plots (PDF) for the primary mixed models
# -----------------------------------------------------------------------------
export_diagnostics_pdf <- function(all_results, corr_cols) {

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
}
