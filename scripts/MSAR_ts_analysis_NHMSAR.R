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
N <- c(4)                     # Number of nodes in the network
                              # Can be a vector for multiple conditions: c(3, 4, 5)

Density <- c(0.25)            # Network edge density (proportion of possible edges)
                              # Range: 0 to 1
                              # Can be a vector: c(0.2, 0.5, 0.8)

# -----------------------------------------------------------------------------
# Regime Parameters
# -----------------------------------------------------------------------------
M <- c(2)                     # Number of regimes (network states)
                              # Can be a vector: c(2, 3, 4)

remain_lower <- 0.33          # Lower bound for probability to stay in same regime
                              # Range: 0 to 1 (higher = more stable regimes)

remain_upper <- 0.66          # Upper bound for probability to stay in same regime
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
T <- c(3500)                  # Number of time steps per time series
                              # Can be a vector: c(1000, 3500, 5000)

warmup <- 50                  # Number of warmup time steps (discarded from analysis)
                              # Allows dynamics to stabilize

totTime <- T + warmup         # Total time steps including warmup

n_ts <- 30                    # Number of time series to generate per condition
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
# Uncomment and set seed for reproducible results
# set.seed(83742)

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
  stay_prob_upper = remain_upper
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
  remain_upper = remain_upper
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
  Timeseries_data = Timeseries_data
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

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Analysis Complete!\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat(sprintf("All results saved to: %s/\n", output_dir))
  cat("\n")
} else {
  cat("Results not saved (save_output = FALSE)\n")
}
