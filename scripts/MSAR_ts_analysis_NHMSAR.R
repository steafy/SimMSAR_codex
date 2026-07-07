# =============================================================================
# SimMSAR Time Series Analysis - Main Execution Script
# =============================================================================

# Load all package dependencies
source("R/dependencies.R")

# Load functions for data generation
source("R/generation/generate_timeseries.R")

# Load functions for model estimation
source("R/estimation/estimate_MSAR.R")


# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================
# Adjust these parameters to customize your simulation and analysis

# -----------------------------------------------------------------------------
# Network Structure Parameters
# -----------------------------------------------------------------------------
N <- c(4, 6, 8)                 # Number of nodes in the network

Density <- c(0.25, 0.5, 0.75)   # Network edge density (proportion of possible edges)
                                # Range: 0 to 1
                                # Can be a vector: c(0.2, 0.5, 0.8)

# -----------------------------------------------------------------------------
# Regime Parameters
# -----------------------------------------------------------------------------
M <- c(1, 2, 3, 4)            # Number of regimes (network states)
                              # Can be a vector: c(2, 3, 4)

remain_lower <- 0.85          # Lower bound for probability to stay in same regime
                              # Range: 0 to 1 (higher = more stable regimes)

remain_upper <- 0.85          # Upper bound for probability to stay in same regime
                              # Range: remain_lower to 1

# -----------------------------------------------------------------------------
# Edge Weight Parameters
# -----------------------------------------------------------------------------
min_edg_val <- 0.05           # Minimum absolute edge weight in networks
                              # Edges below this threshold are set to zero

max_edg_val <- 1              # Maximum absolute edge weight in networks

# -----------------------------------------------------------------------------
# Time Series Parameters
# -----------------------------------------------------------------------------
T <- c(200, 400, 800, 1600)   # Number of time steps per time series
                              # Can be a vector: c(1000, 3500, 5000)
 
warmup <- 50                  # Number of warmup time steps (discarded from analysis)
                              # Allows dynamics to stabilize

totTime <- T + warmup         # Total time steps including warmup

n_ts <- 150                   # Number of time series to generate per condition
                              # Higher = more statistical power but slower

# -----------------------------------------------------------------------------
# Model Estimation Parameters
# -----------------------------------------------------------------------------
order <- 1                    # Lag order for autoregression (AR order)
                              # Usually 1 for first-order VAR models

MaxIter <- 200                # Maximum iterations for EM algorithm
                              # Higher = more likely to converge but slower

eps <- 1e-5                   # Convergence criterion (epsilon)
                              # Smaller = stricter convergence

retry_attempts <- 5           # Number of retry attempts if fitting fails

# -----------------------------------------------------------------------------
# LASSO / Penalization Parameters (first M-step of the EM algorithm)
# -----------------------------------------------------------------------------
# Control how the sparse network (Beta) is estimated in the first EM M-step.
# Background & validation: docs/MSTEP_LASSO_CV_PENALIZATION.md.
#
# Two engines are available:
#   "bic"      Legacy lars-subset + BIC + weighted-OLS refit. The L1 penalty has
#              NO real effect: it returns a (near-)saturated network, and all
#              sparsity comes from the downstream min_edg_val threshold applied
#              in the analysis pipeline. Fast (~1x).
#   "cvglmnet" Genuinely penalized: glmnet::cv.glmnet with correctly weighted
#              observations. Produces a sparse network from the ESTIMATOR itself
#              (before any threshold). To match/beat the legacy engine's recovery
#              it MUST re-select the support every EM iteration (see below);
#              with that on it recovers Beta/Kappa/Sigma as well or better AND
#              gives threshold-free sparsity, at ~5x the legacy runtime.
#
# DEFAULT below = the validated "cvglmnet + re-select + fast CV" configuration.
# Set lasso_engine <- "bic" to reproduce the previous (legacy) behaviour exactly.

lasso_engine     <- "cvglmnet" # "cvglmnet" (genuine LASSO, default) | "bic" (legacy)

lasso_reselect   <- TRUE       # Re-run the penalized support selection EVERY EM
                               # iteration instead of freezing iteration 1's
                               # support. REQUIRED for good recovery with
                               # "cvglmnet": a frozen support (FALSE) is chosen
                               # from the poor initial E-step and locks in missed
                               # edges -> worse than legacy. Ignored for "bic".

lasso_lambda     <- "1se"      # CV lambda rule: "1se" (sparser, recommended) or
                               # "min" (denser, barely sparse). cvglmnet only.

lasso_refit      <- TRUE       # TRUE = relaxed LASSO: re-estimate coefficients on
                               # the selected support by unpenalized weighted OLS
                               # (removes L1 shrinkage bias in the edge weights).
                               # FALSE = use the shrunk glmnet coefficients
                               # directly. cvglmnet only.

lasso_adaptive   <- FALSE      # Adaptive LASSO penalty (penalty.factor=1/|OLS|).
                               # Tested and NOT helpful here (over-penalizes weak
                               # true edges near min_edg_val); keep FALSE.

# --- cv.glmnet cost knobs (cvglmnet only) ------------------------------------
# Because coefficients are refit by OLS (lasso_refit=TRUE), cv.glmnet only has to
# pick the SUPPORT, so a cheap CV suffices. FIXED folds are important under
# re-selection: random folds make the support flicker between iterations and
# blow up the iteration count (5-fold random ~28 EM iters vs ~7 fixed). The
# defaults below are ~2.5x faster than a naive 10-fold/100-lambda CV with
# identical recovery. Going below 50 lambda or 3 folds starts to lose recovery.
lasso_nfolds     <- 5          # Number of CV folds (5 recommended; >=3).
lasso_nlambda    <- 50         # Length of the lambda grid (50 recommended; >=50).
lasso_fixedfolds <- TRUE       # Deterministic (fixed) fold partition every EM
                               # iteration -> stable support -> fast convergence.

# --- re-selection schedule (advanced; cvglmnet + lasso_reselect only) --------
# Full re-selection (defaults below) is the most robust. Partial schedules were
# tested and are generally NOT worth it (they trade convergence stability for a
# small speedup); see docs §5.1. Leave at the defaults unless experimenting.
lasso_reselect_iters <- Inf    # Re-select only for the first k EM iterations,
                               # then freeze. Inf = re-select every iteration.
lasso_reselect_every <- 1      # Among re-selecting iterations, re-select every
                               # m-th one. 1 = every iteration.

# -----------------------------------------------------------------------------
# Residual-covariance (Sigma) stabilization
# -----------------------------------------------------------------------------
# Optional in-EM regularization of the regime-weighted residual covariance Sigma,
# to stop it becoming pathologically ill-conditioned when a regime is assigned very
# few effective observations (small postmix) relative to its d(d+1)/2 covariance
# parameters -- the under-determination that produces exploding Kappa = solve(Sigma)
# estimates and outright fit failures. Diagnosis & validation:
# docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md. Applied on EVERY M-step; keeps Sigma
# dense (never introduces zeros). Default "none" reproduces the previous behaviour.
#
#   "none"   no stabilization (default; production behaviour unchanged).
#   "floor"  eigenvalue floor: clip eigenvalues below floor*max(eig), which bounds
#            the condition number at exactly 1/sigma_stab_floor and leaves any Sigma
#            already better-conditioned than that completely untouched. RECOMMENDED.
#   "ridge"  Sigma + lambda*(tr Sigma/d)*I: lifts every eigenvalue; converges a bit
#            more reliably but perturbs even well-conditioned regimes slightly.
sigma_stab       <- "floor"     # "none" (default) | "floor" (recommended) | "ridge"

sigma_stab_floor <- 1e-2       # floor only: condition-number cap = 1/floor.
                               # 1e-3 -> cap 1000 (safe default: removes fit
                               # failures, converges reliably). 1e-2 -> cap 100
                               # (tightest control / best NRMSE_Kappa, but can hit
                               # MaxIter benignly on borderline fits).

sigma_stab_lambda <- 1e-2      # ridge only: relative ridge strength (fraction of
                               # the matrix's own average variance tr(Sigma)/d).

# Assemble the control list passed to estimate_MSAR(). These are applied via
# options() INSIDE each parallel worker (options do not propagate to future
# workers automatically), so this is the correct way to configure a parallel run.
lasso_control <- list(
  simmsar_lasso_engine         = lasso_engine,
  simmsar_lasso_reselect       = lasso_reselect,
  simmsar_lasso_lambda         = lasso_lambda,
  simmsar_lasso_refit          = lasso_refit,
  simmsar_lasso_adaptive       = lasso_adaptive,
  simmsar_lasso_nfolds         = lasso_nfolds,
  simmsar_lasso_nlambda        = lasso_nlambda,
  simmsar_lasso_fixedfolds     = lasso_fixedfolds,
  simmsar_lasso_reselect_iters = lasso_reselect_iters,
  simmsar_lasso_reselect_every = lasso_reselect_every,
  # residual-covariance stabilization (see block above)
  simmsar_sigma_stab           = sigma_stab,
  simmsar_sigma_stab_floor     = sigma_stab_floor,
  simmsar_sigma_stab_lambda    = sigma_stab_lambda
)
# Also apply in the main process (covers sequential runs / interactive fits;
# the parallel workers get their own copy via estimate_MSAR(lasso_control=...)).
do.call(options, lasso_control)

# -----------------------------------------------------------------------------
# Parallelisierung
# -----------------------------------------------------------------------------
workers <- 5   # NULL = automatisch (physische Kerne - 1); explizit z.B. workers <- 5

# -----------------------------------------------------------------------------
# Output Parameters
# -----------------------------------------------------------------------------
verbose <- FALSE              # Print detailed progress messages
                              # Set to TRUE for debugging

save_output <- TRUE           # Save results to .rds files
output_dir <- "output"        # Directory for saving results

# Create output directory if it doesn't exist
if (save_output && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------
# Set seed for reproducible results
# IMPORTANT: This must be set BEFORE generate_timeseries() is called
# All random operations (network generation, time series simulation,
# model initialization) will be reproducible with this single seed
set.seed(20260701)

# =============================================================================
# PARAMETER SUMMARY
# =============================================================================
# Print configuration for verification

CONFIG <- list(
  # Network structure
  nodes = N,
  density = Density,
  regimes = M,

  # Edge weights
  min_edge = min_edg_val,
  max_edge = max_edg_val,

  # Time series
  timesteps = T,
  warmup = warmup,
  n_timeseries = n_ts,

  # Model parameters
  ar_order = order,
  max_iterations = MaxIter,
  convergence_eps = eps,

  # Regime dynamics
  stay_prob_lower = remain_lower,
  stay_prob_upper = remain_upper,

  # LASSO / penalization (first M-step)
  lasso_control = lasso_control
)

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  SimMSAR Configuration\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("Network:      %d nodes, %.0f%% density, %d regime(s)\n",
            N[1], Density[1]*100, M[1]))
cat(sprintf("Time Series:  %d timesteps, %d series per condition\n",
            T[1], n_ts))
cat(sprintf("Estimation:   AR(%d), MaxIter=%d, eps=%.0e\n",
            order, MaxIter, eps))
cat(sprintf("Edge Range:   [%.2f, %.2f]\n",
            min_edg_val, max_edg_val))
cat(sprintf("Stay Prob:    [%.2f, %.2f]\n",
            remain_lower, remain_upper))
if (identical(lasso_engine, "cvglmnet")) {
  cat(sprintf("LASSO:        cvglmnet, re-select=%s, lambda.%s, refit=%s, %d folds x %d lambda, fixed folds=%s\n",
              lasso_reselect, lasso_lambda, lasso_refit, lasso_nfolds, lasso_nlambda, lasso_fixedfolds))
} else {
  cat("LASSO:        bic (legacy; sparsity via min_edg_val threshold only)\n")
}
if (identical(sigma_stab, "floor")) {
  cat(sprintf("Sigma stab:   floor (eigenvalue floor; condition cap = %.0f)\n", 1 / sigma_stab_floor))
} else if (identical(sigma_stab, "ridge")) {
  cat(sprintf("Sigma stab:   ridge (lambda = %.1e of tr(Sigma)/d)\n", sigma_stab_lambda))
} else {
  cat("Sigma stab:   none (raw residual covariance; Kappa = solve(Sigma) unregularized)\n")
}
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")

# =============================================================================
# END CONFIGURATION
# =============================================================================



# =============================================================================
# DATA GENERATION
# =============================================================================

cat("Starting data generation...\n")
Timeseries_data <- generate_timeseries(
  Density = Density,
  min_edg_val = min_edg_val,
  max_edg_val = max_edg_val,
  M = M,
  N = N,
  warmup = warmup,
  T = T,
  totTime = totTime,
  n_ts = n_ts,
  remain_lower = remain_lower,
  remain_upper = remain_upper,
  workers = workers
  )

cat("Data generation complete!\n\n")

# =============================================================================
# MODEL ESTIMATION
# =============================================================================

cat("Starting model estimation...\n")
MSAR_dynamics_list <- estimate_MSAR(
  Density = Density,
  M = M,
  N = N,
  T = T,
  n_ts = n_ts,
  order = order,
  MaxIter = MaxIter,
  verbose = verbose,
  min_edg_val = min_edg_val,
  Timeseries_data = Timeseries_data,
  workers = workers,
  lasso_control = lasso_control   # LASSO engine / re-selection settings (see above)
)

cat("Model estimation complete!\n\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

if (save_output) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Save time series data
  ts_filename <- file.path(output_dir, sprintf("Timeseries_data_%s.rds", timestamp))
  saveRDS(Timeseries_data, ts_filename)
  cat(sprintf("✓ Saved time series data: %s\n", ts_filename))

  # Save MSAR models
  model_filename <- file.path(output_dir, sprintf("MSAR_models_%s.rds", timestamp))
  saveRDS(MSAR_dynamics_list, model_filename)
  cat(sprintf("✓ Saved MSAR models: %s\n", model_filename))

  # Save configuration
  config_filename <- file.path(output_dir, sprintf("CONFIG_%s.rds", timestamp))
  saveRDS(CONFIG, config_filename)
  cat(sprintf("✓ Saved configuration: %s\n", config_filename))

  # Save session info for reproducibility
  session_filename <- file.path(output_dir, sprintf("sessionInfo_%s.txt", timestamp))
  writeLines(capture.output(sessionInfo()), session_filename)
  cat(sprintf("✓ Saved session info: %s\n", session_filename))

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Analysis Complete!\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat(sprintf("All results saved to: %s/\n", output_dir))
  cat(sprintf("Session info saved for reproducibility\n"))
  cat("\n")
} else {
  cat("Results not saved (save_output = FALSE)\n")
}
