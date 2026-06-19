# =============================================================================
# Feasibility model: probability of successful estimation (convergence)
# =============================================================================
# The "two-part" companion to the recovery LMMs: where the LMMs model recovery
# QUALITY *conditional on* a successful fit, this models whether a fit succeeds
# at all, as a function of the same design factors. Together they separate
# "can you estimate it?" (here) from "if so, how well?" (modeling.R).
#
# Unit of analysis: one row per design condition. The response is the grouped
# binomial cbind(successes, failures) -- equivalent to a per-time-series
# Bernoulli GLM since predictors are constant within a condition. Predictors are
# standardised on the condition grid, mirroring the LMM coding (logT, *_s,
# Regimes factor), so coefficients are log-odds of SUCCESS per grid-SD.

# Returns list(model, coef_table, family) or NULL if not applicable.
fit_feasibility_model <- function(failure_analysis) {

  cat("\n=== FEASIBILITY MODEL: P(successful estimation) ~ design factors ===\n\n")

  if (sum(failure_analysis$failures) == 0) {
    cat("All conditions converged (0 failures). Feasibility model not applicable.\n\n")
    return(NULL)
  }
  if (all(failure_analysis$failures > 0) && all(failure_analysis$N == 0)) {
    cat("No successful fits at all. Feasibility model not applicable.\n\n")
    return(NULL)
  }

  # Build modelling frame: numeric predictors (+ standardised) and Regimes factor.
  df <- failure_analysis %>%
    mutate(
      Timesteps_num = as.numeric(as.character(Timesteps)),
      Density_num   = as.numeric(as.character(Density)),
      Nodes_num     = as.numeric(as.character(Nodes)),
      Regimes_num   = as.numeric(as.character(Regimes)),
      n_success     = N,
      n_fail        = failures
    )
  df$Regimes <- factor(df$Regimes_num)

  # Guard against (quasi-)complete separation. A Regimes level that NEVER fails
  # perfectly predicts success and sends the GLM coefficients to +/-Inf
  # (huge estimates, p ~ 1). Single-regime (M=1) fits have no regime-allocation
  # failure mode and essentially always converge, so they typically trigger
  # this. Such all-success levels carry no information about *variation* in
  # feasibility: we report them as "always converged" and exclude them.
  fail_by_regime <- tapply(df$n_fail, df$Regimes, sum)
  always_ok <- names(fail_by_regime)[fail_by_regime == 0]
  if (length(always_ok) > 0) {
    cat(sprintf("Regimes level(s) with 0 failures (always converged, excluded to avoid separation): %s\n\n",
                paste0("M=", always_ok, collapse = ", ")))
    df <- df[!(as.character(df$Regimes) %in% always_ok), , drop = FALSE]
    df$Regimes <- droplevels(df$Regimes)
  }

  # Standardise predictors on the (post-exclusion) condition grid -> per-SD effects
  df$logT      <- as.numeric(scale(log(df$Timesteps_num)))
  df$Density_s <- as.numeric(scale(df$Density_num))
  df$Nodes_s   <- as.numeric(scale(df$Nodes_num))

  # Drop Regimes from the RHS if it now has only one level
  rhs <- c("logT", "Density_s", "Nodes_s",
           if (nlevels(df$Regimes) > 1) "Regimes")
  form <- as.formula(paste("cbind(n_success, n_fail) ~", paste(rhs, collapse = " + ")))

  model <- tryCatch(
    glm(form, data = df, family = binomial),
    error = function(e) { cat("Feasibility GLM failed:", e$message, "\n"); NULL }
  )
  if (is.null(model)) return(NULL)

  # Residual-separation check: even after dropping all-success levels, a steep
  # continuous predictor can quasi-separate. Flag it so the estimate is not
  # taken at face value (Firth's penalised logistic regression would be the fix).
  if (max(abs(coef(model)), na.rm = TRUE) > 15) {
    cat("⚠ Possible (quasi-)complete separation remains (extreme coefficients).\n")
    cat("  Treat the affected estimate with caution; consider Firth's penalised\n")
    cat("  logistic regression (e.g. logistf::logistf or brglm2).\n\n")
  }

  # Overdispersion check: residual deviance / residual df. If clearly > 1, refit
  # quasibinomial so the SEs (and p-values) are not anticonservative.
  rdf  <- df.residual(model)
  disp <- if (rdf > 0) sum(residuals(model, type = "pearson")^2) / rdf else NA_real_
  fam  <- "binomial"
  if (!is.na(disp) && disp > 1.5) {
    cat(sprintf("Overdispersion detected (Pearson dispersion = %.2f); refitting quasibinomial.\n\n",
                disp))
    model <- glm(form, data = df, family = quasibinomial)
    fam   <- "quasibinomial"
  } else if (!is.na(disp)) {
    cat(sprintf("Pearson dispersion = %.2f (no overdispersion correction needed).\n\n", disp))
  }

  # Coefficient table: log-odds, odds ratio, Wald 95% CI, p-value, stars.
  sm  <- summary(model)$coefficients
  est <- sm[, "Estimate"]
  se  <- sm[, "Std. Error"]
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  p   <- sm[, pcol]

  coef_table <- data.frame(
    Predictor = rownames(sm),
    logOR     = round(est, 3),
    OR        = round(exp(est), 3),
    CI_low    = round(exp(est - 1.96 * se), 3),
    CI_high   = round(exp(est + 1.96 * se), 3),
    p         = signif(p, 3),
    sig       = ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", ""))),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  # Clean predictor names for display (match the LMM coefficient-table style)
  coef_table$Predictor <- gsub("_s\\b", "", coef_table$Predictor)
  coef_table$Predictor <- gsub("logT", "log T", coef_table$Predictor)

  cat("Response: P(successful estimation). Positive logOR -> higher convergence odds.\n")
  cat("Predictors standardised on the condition grid (per-SD effects).\n\n")
  print(coef_table, row.names = FALSE)
  cat(sprintf("\nFamily: %s | conditions (rows): %d | total fits: %d (%d failed)\n\n",
              fam, nrow(df), sum(df$n_success + df$n_fail), sum(df$n_fail)))

  # LaTeX export (thesis-ready), consistent with the other .tex tables.
  if (requireNamespace("xtable", quietly = TRUE)) {
    tex_tbl <- coef_table
    # Thesis-friendly p display: "<0.001" instead of a rounded 0.0000
    tex_tbl$p <- ifelse(coef_table$p < 0.001, "$<0.001$",
                        formatC(coef_table$p, format = "f", digits = 3))
    colnames(tex_tbl) <- c("Predictor", "$\\log\\text{OR}$", "OR",
                           "CI$_{2.5}$", "CI$_{97.5}$", "$p$", "")
    print(
      xtable::xtable(
        tex_tbl,
        caption = paste0("Feasibility model: log-odds of successful estimation ",
                         "as a function of the design factors (", fam, ")."),
        digits = c(0, 0, 3, 3, 3, 3, 0, 0)
      ),
      file = "feasibility_model.tex",
      include.rownames = FALSE,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity
    )
    cat("Exported: feasibility_model.tex\n\n")
  }

  list(model = model, coef_table = coef_table, family = fam)
}
