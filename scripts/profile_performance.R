# =============================================================================
# Performance Profiling Script for SimMSAR
# =============================================================================
# This script profiles a small but representative simulation to identify
# performance bottlenecks in the generation and estimation pipeline.

# Install profiling tools if needed
if (!require("profvis", quietly = TRUE)) {
  install.packages("profvis")
}
library(profvis)

# Load dependencies
source("R/dependencies.R")

# Load functions
source("R/generation/generate_timeseries.R")
source("R/estimation/estimate_MSAR.R")

cat("=================================================================\n")
cat("  SimMSAR Performance Profiling\n")
cat("=================================================================\n\n")

# =============================================================================
# PROFILING CONFIGURATION
# =============================================================================
# Small but representative workload
# Adjust these to profile different scenarios

# Network parameters (small for quick profiling)
N <- 4              # Nodes
Density <- 0.3      # Edge density
M <- 2              # Regimes

# Time series parameters
T <- 500            # Timesteps (reduced from 3500 for quick profiling)
n_ts <- 3           # Number of time series (reduced from 30)
warmup <- 50

# Edge weights
min_edg_val <- 0.05
max_edg_val <- 1

# Regime stability
remain_lower <- 0.33
remain_upper <- 0.66

# Estimation parameters
order <- 1
MaxIter <- 50       # Reduced from 200 for quick profiling
verbose <- FALSE

cat("Configuration:\n")
cat(sprintf("  Network: %d nodes, %.0f%% density, %d regimes\n", N, Density*100, M))
cat(sprintf("  Time Series: %d timesteps, %d series\n", T, n_ts))
cat(sprintf("  Estimation: MaxIter=%d\n", MaxIter))
cat("\n")

# =============================================================================
# PROFILE DATA GENERATION
# =============================================================================

cat("Profiling data generation...\n")

profvis({
  Timeseries_data <- generate_timeseries(
    Density = Density,
    min_edg_val = min_edg_val,
    max_edg_val = max_edg_val,
    M = M,
    N = N,
    warmup = warmup,
    T = T,
    totTime = T + warmup,
    n_ts = n_ts,
    remain_lower = remain_lower,
    remain_upper = remain_upper
  )
}, interval = 0.005) -> profile_generation

cat("✓ Generation profiling complete\n\n")

# Save generation profile
htmlwidgets::saveWidget(profile_generation, "profile_generation.html")
cat("Saved: profile_generation.html\n\n")

# =============================================================================
# PROFILE MODEL ESTIMATION
# =============================================================================

cat("Profiling model estimation...\n")

profvis({
  MSAR_models <- estimate_MSAR(
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
}, interval = 0.005) -> profile_estimation

cat("✓ Estimation profiling complete\n\n")

# Save estimation profile
htmlwidgets::saveWidget(profile_estimation, "profile_estimation.html")
cat("Saved: profile_estimation.html\n\n")

# =============================================================================
# SIMPLE TIMING COMPARISON
# =============================================================================

cat("Running timing benchmarks...\n\n")

# Benchmark generation
gen_times <- numeric(3)
for (i in 1:3) {
  start_time <- Sys.time()
  Timeseries_data <- generate_timeseries(
    Density = Density,
    min_edg_val = min_edg_val,
    max_edg_val = max_edg_val,
    M = M,
    N = N,
    warmup = warmup,
    T = T,
    totTime = T + warmup,
    n_ts = n_ts,
    remain_lower = remain_lower,
    remain_upper = remain_upper
  )
  gen_times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
}

cat("Data Generation Timing:\n")
cat(sprintf("  Mean: %.2f seconds\n", mean(gen_times)))
cat(sprintf("  Range: [%.2f, %.2f] seconds\n", min(gen_times), max(gen_times)))
cat("\n")

# Benchmark estimation (single run due to longer time)
cat("Running estimation benchmark (1 run)...\n")
start_time <- Sys.time()
MSAR_models <- estimate_MSAR(
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
est_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("Model Estimation Timing:\n")
cat(sprintf("  Time: %.2f seconds\n", est_time))
cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=================================================================\n")
cat("  Profiling Summary\n")
cat("=================================================================\n")
cat(sprintf("Data Generation:  %.2f sec (%.0f%%)\n",
            mean(gen_times),
            100 * mean(gen_times) / (mean(gen_times) + est_time)))
cat(sprintf("Model Estimation: %.2f sec (%.0f%%)\n",
            est_time,
            100 * est_time / (mean(gen_times) + est_time)))
cat(sprintf("Total:            %.2f sec\n", mean(gen_times) + est_time))
cat("\n")
cat("Profiling reports saved:\n")
cat("  - profile_generation.html\n")
cat("  - profile_estimation.html\n")
cat("\n")
cat("Open these files in a browser to view interactive flame graphs.\n")
cat("=================================================================\n")
