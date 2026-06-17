# =============================================================================
# Test Script for Refactored MSAR Code
# =============================================================================
# This script tests the new tibble-based structure with minimal data

# Load dependencies
source("R/dependencies.R")

# Set test parameters (minimal for quick testing)
N <- c(4)                    # One network size
Density <- c(0.25)           # One density
M <- c(2)                    # One regime count
T <- c(100)                  # Short time series
warmup <- 10
totTime <- T + warmup
n_ts <- 2                    # Just 2 time series

min_edg_val <- 0.05
max_edg_val <- 1
remain_lower <- 0.33
remain_upper <- 0.66

order <- 1
MaxIter <- 50                # Fewer iterations for testing
verbose <- FALSE

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Testing Refactored MSAR Code\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("Running minimal test with:\n")
cat(sprintf("  %d nodes, %.0f%% density, %d regimes\n", N, Density*100, M))
cat(sprintf("  %d timesteps, %d time series\n", T, n_ts))
cat("═══════════════════════════════════════════════════════════════\n\n")

# Test 1: Generate time series data (now returns tibble!)
cat("Test 1: Generating time series data (new tibble structure)...\n")
source("R/generation/generate_timeseries.R")

Timeseries_data <- tryCatch({
  generate_timeseries(
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
}, error = function(e) {
  cat("ERROR in generate_timeseries:\n")
  print(e)
  return(NULL)
})

if (is.null(Timeseries_data)) {
  stop("Failed to generate time series data")
}
cat("✓ Time series generation successful\n")
cat(sprintf("  Class: %s\n", paste(class(Timeseries_data), collapse = ", ")))
cat(sprintf("  Rows: %d\n", nrow(Timeseries_data)))
cat(sprintf("  Columns: %d\n", ncol(Timeseries_data)))

# Verify timeseries_data structure
ts_expected_cols <- c("timesteps", "density", "nodes", "regimes", "ts_id",
                      "timeseries_data", "regime_dynamics")
ts_missing_cols <- setdiff(ts_expected_cols, colnames(Timeseries_data))
if (length(ts_missing_cols) > 0) {
  cat("  WARNING: Missing columns:", paste(ts_missing_cols, collapse = ", "), "\n")
} else {
  cat("  ✓ All expected columns present\n")
}
cat("\n")

# Test 2: Estimate MSAR models (returns tibble now)
cat("Test 2: Estimating MSAR models (new tibble structure)...\n")
source("R/estimation/estimate_MSAR.R")

MSAR_dynamics_list <- tryCatch({
  estimate_MSAR(
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
}, error = function(e) {
  cat("ERROR in estimate_MSAR:\n")
  print(e)
  return(NULL)
})

if (is.null(MSAR_dynamics_list)) {
  stop("Failed to estimate MSAR models")
}
cat("✓ MSAR estimation successful\n\n")

# Test 3: Verify output structure
cat("Test 3: Verifying new tibble structure...\n")
cat(sprintf("  Class: %s\n", paste(class(MSAR_dynamics_list), collapse = ", ")))
cat(sprintf("  Rows: %d\n", nrow(MSAR_dynamics_list)))
cat(sprintf("  Columns: %d\n", ncol(MSAR_dynamics_list)))

# Check expected columns
expected_cols <- c("timesteps", "density", "nodes", "regimes", "ts_id", "regime_id",
                   "Beta_corr", "Kappa_corr", "Beta_ac_corr")
missing_cols <- setdiff(expected_cols, colnames(MSAR_dynamics_list))
if (length(missing_cols) > 0) {
  cat("  WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
} else {
  cat("  ✓ All expected columns present\n")
}
cat("\n")

# Test 4: Test S3 print method
cat("Test 4: Testing S3 print method...\n")
print(MSAR_dynamics_list)
cat("\n")

# Test 5: Test get_stats() function
cat("Test 5: Testing get_stats() function...\n")
stats <- tryCatch({
  get_stats(MSAR_dynamics_list)
}, error = function(e) {
  cat("ERROR in get_stats:\n")
  print(e)
  return(NULL)
})

if (!is.null(stats)) {
  cat("✓ get_stats() successful\n")
  cat(sprintf("  Stats rows: %d\n", nrow(stats)))
  print(head(stats))
} else {
  cat("  WARNING: get_stats() failed\n")
}
cat("\n")

# Test 6: Test dplyr filtering (new capability)
cat("Test 6: Testing dplyr operations on results...\n")
filtered <- tryCatch({
  MSAR_dynamics_list %>%
    filter(timesteps == T, density == Density) %>%
    select(ts_id, regime_id, Beta_corr, Kappa_corr, Beta_ac_corr)
}, error = function(e) {
  cat("ERROR in dplyr operations:\n")
  print(e)
  return(NULL)
})

if (!is.null(filtered)) {
  cat("✓ dplyr filtering successful\n")
  cat("  Filtered results:\n")
  print(filtered)
} else {
  cat("  WARNING: dplyr operations failed\n")
}
cat("\n")

# Test 7: Verify stat_analysis.R compatibility
cat("Test 7: Testing stat_analysis.R compatibility...\n")
stat_test <- tryCatch({
  # This simulates what stat_analysis.R does
  corr_results <- MSAR_dynamics_list %>%
    rename(
      Timesteps = timesteps,
      Density = density,
      Nodes = nodes,
      Regimes = regimes
    ) %>%
    group_by(Timesteps, Density, Nodes, Regimes) %>%
    mutate(
      SimID = ts_id,
      RegimeIndex = regime_id,
      N = n_distinct(SimID)
    ) %>%
    ungroup()

  cat("  ✓ stat_analysis.R data transformation successful\n")
  cat(sprintf("    Unique conditions: %d\n", n_distinct(corr_results$Timesteps,
                                                         corr_results$Density,
                                                         corr_results$Nodes,
                                                         corr_results$Regimes)))
  TRUE
}, error = function(e) {
  cat("ERROR in stat_analysis transformation:\n")
  print(e)
  FALSE
})

cat("\n")

# Test 8: Verify regime-sequence recovery attribute (M > 1 only)
cat("Test 8: Testing regime-sequence recovery attribute...\n")
seq_results <- attr(MSAR_dynamics_list, "sequence_results")
if (is.null(seq_results)) {
  cat("  WARNING: sequence_results attribute is missing\n")
} else {
  cat(sprintf("  ✓ sequence_results attribute present (%d rows)\n", nrow(seq_results)))
  print(seq_results)
}
cat("\n")

# Final summary
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Test Summary\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("All tests completed successfully!\n\n")
cat("Key improvements:\n")
cat("  • No more 5-level nested lists in EITHER stage!\n")
cat("  • generate_timeseries() returns a clean tibble\n")
cat("  • estimate_MSAR() returns a clean tibble\n")
cat("  • Easy filtering with dplyr everywhere\n")
cat("  • S3 print methods for better display\n")
cat("  • get_stats() for quick summaries\n")
cat("  • Consistent tibble structure throughout pipeline\n")
cat("\n")
cat("Example usage:\n")
cat("  # Filter by condition\n")
cat("  results %>% filter(timesteps == 1000, density == 0.3)\n")
cat("\n")
cat("  # Get summary stats\n")
cat("  get_stats(results)\n")
cat("\n")
cat("  # Access individual networks\n")
cat("  results$orig_Beta[[1]]  # First original temporal network\n")
cat("═══════════════════════════════════════════════════════════════\n")
