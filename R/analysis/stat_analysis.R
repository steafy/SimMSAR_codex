# =============================================================================
# Statistical Analysis for MSAR Models -- Orchestrator
# =============================================================================
# This script ties together the analysis modules in R/analysis/. It expects a
# `MSAR_dynamics_list` object (class "msar_results", as returned by
# estimate_MSAR()) to be present in the calling environment, e.g.:
#
#   MSAR_dynamics_list <- readRDS("output/MSAR_models_<timestamp>.rds")
#   source("R/analysis/stat_analysis.R")
#
# RANDOM EFFECTS STRUCTURE (see R/analysis/modeling.R):
# The random-effects structure is selected adaptively per outcome:
# 1. Tests whether random slopes for logT improve fit (LR test)
# 2. If yes, checks for boundary singularity in the correlated RE model
# 3. If singular (rho ~ 1.0), uses uncorrelated RE to ensure stability
# 4. Otherwise, uses correlated RE (best fit)
# =============================================================================

# -----------------------------------------------------------------------------
# SETUP: Load required packages
# -----------------------------------------------------------------------------
required_packages <- c("dplyr", "ggplot2", "kableExtra", "xtable", "lmtest", "sandwich",
                       "lme4", "lmerTest", "performance", "dunn.test", "cowplot")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# -----------------------------------------------------------------------------
# Load analysis modules
# -----------------------------------------------------------------------------
source("R/analysis/data_prep.R")      # PART 1: scoring, data prep + selection-bias reporting
source("R/analysis/feasibility.R")    # PART 1.4: feasibility (convergence) model
source("R/analysis/transform.R")      # PART 2-4: Fisher-z, aggregation, scaling
source("R/analysis/descriptives.R")   # PART 3.1: descriptive statistics
source("R/analysis/modeling.R")       # PART 5: mixed-effects model fitting
source("R/analysis/sensitivity.R")    # PART 5.5: sensitivity analysis
source("R/analysis/sequence_models.R")# PART 5.6: RQ4 regime-sequence recovery LMM
source("R/analysis/exports.R")        # PART 6: coefficient/comparison tables, diagnostics
source("R/analysis/plots.R")          # PART 7-8: line, sens/spec, coefficient plots
source("R/analysis/reporting.R")      # FINAL: model-quality summary

options(warn = 1)

# Scoring constants (now that metrics are computed in analysis, not estimation).
# These are the revisable scoring decisions the refactor was about:
MIN_EDG_VAL    <- 0.05   # edge-detection threshold; must match the generation intent
AC_HORIZON     <- 25     # T_ac for average_controllability(); MUST equal generation's default
KAPPA_COND_MAX <- 1e6    # max condition number of est_Sigma before Kappa is flagged invalid

# Outcome variables analysed throughout. Decided (test run 2026-06): Pearson is
# the primary AC-recovery outcome. Spearman AC correlation is computed in
# compute_recovery_metrics() (column Beta_ac_corr_spearman stays available on
# MSAR_dynamics_list/corr_results for ad-hoc robustness reporting) but is no
# longer run through Fisher-z/LMM/exports/plots: at N=4 nodes Spearman's rho is
# discretised onto a handful of values, so |rho|=1 occurs from discreteness
# alone -- 18.1% of regime rows had |z|>5 vs. 1.5% for Pearson in the pilot run,
# and R2_m was correspondingly worse (0.404 vs 0.560). Pearson AC correlation
# now gets the same precision-weighting treatment as Kappa_corr (see
# weight_BetaAC in transform.R / OUTCOME_WEIGHT_MAP in modeling.R), since it
# inherits the analogous low-n artifact (n = Nodes instead of n_nz).
corr_cols <- c("Beta_corr", "Kappa_corr", "Beta_ac_corr_pearson")

# -----------------------------------------------------------------------------
# PIPELINE
# -----------------------------------------------------------------------------

# PART 1.0: Score the raw estimates -- edge threshold, Kappa = solve(Sigma) +
# condition-number validity gate, correlations, sens/spec, NRMSE, AC corrs.
MSAR_dynamics_list <- compute_recovery_metrics(
  MSAR_dynamics_list,
  min_edg_val    = MIN_EDG_VAL,
  AC_HORIZON     = AC_HORIZON,
  KAPPA_COND_MAX = KAPPA_COND_MAX
)
sigma_validity <- summarize_sigma_validity(MSAR_dynamics_list)

# PART 1: Data preparation and estimation-process / selection-bias analysis
corr_results <- prepare_corr_results(MSAR_dynamics_list)
est          <- analyze_estimation_process(corr_results)
bias         <- report_selection_bias(est)

# PART 1.4: Feasibility model -- P(successful estimation) ~ design factors.
# The "can you estimate it?" half of the two-part report (complements the
# recovery LMMs, which are conditional on a successful fit).
feasibility  <- fit_feasibility_model(est$failure_analysis)

# PART 2: Fisher-z transformation of correlation outcomes
corr_results <- add_fisher_z(corr_results, corr_cols)

# PART 3: Aggregation to simulation level
agg     <- aggregate_to_sim_level(corr_results, corr_cols)
dat_sim <- agg$dat_sim

# PART 3.1: Descriptive statistics
descript_stats <- descriptive_stats(dat_sim, corr_cols, agg$available_other_metrics)

# PART 4: Transform & scale predictors; attach success rates
dat_sim <- scale_predictors(dat_sim, est$selection_weights)

# PART 5: Fit mixed-effects models
all_results <- fit_all_models(dat_sim, corr_cols)

# PART 5.5: Sensitivity analysis
sens_results <- run_sensitivity_analysis(dat_sim, all_results, corr_cols, bias)

# PART 5.6: RQ4 -- regime-sequence recovery LMM (Cohen's kappa), from the
# fit-level table attached by estimate_MSAR().
seq_recovery <- fit_sequence_recovery_model(attr(MSAR_dynamics_list, "fit_results"))

# PART 6: Export results
export_coefficient_tables(all_results, corr_cols)
comparison_all <- export_comparison_table(all_results, corr_cols)
export_diagnostics_pdf(all_results, corr_cols)

# PART 7: Visualisations
make_line_plots(corr_results, corr_cols)
make_senspec_plots(corr_results)

# FINAL SUMMARY
print_final_summary(all_results, corr_cols, bias)

# PART 8: Coefficient plots
make_coefficient_plots(all_results, corr_cols)
