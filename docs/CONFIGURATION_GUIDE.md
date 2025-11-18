# Configuration Guide

## Overview

All simulation and analysis parameters are configured at the beginning of the main script (`scripts/MSAR_ts_analysis_NHMSAR.R`). This centralized configuration makes it easy to adjust parameters without modifying code throughout the project.

## Configuration Sections

### 1. Network Structure Parameters

```r
N <- c(4)          # Number of nodes in the network
Density <- c(0.25) # Network edge density
```

**N (Number of Nodes)**
- Type: Integer or vector of integers
- Range: Positive integers (typically 3-10)
- Purpose: Defines the size of the network
- Multiple values: Use vector for multiple conditions, e.g., `c(3, 4, 5)`
- Impact: Larger networks require more computational time

**Density (Edge Density)**
- Type: Numeric or vector
- Range: 0 to 1
- Purpose: Proportion of possible edges that exist
- Examples:
  - `0.2` = Sparse network (20% of edges)
  - `0.5` = Medium density
  - `0.8` = Dense network (80% of edges)
- Multiple values: `c(0.2, 0.5, 0.8)` for comparing densities

### 2. Regime Parameters

```r
M <- c(2)              # Number of regimes
remain_lower <- 0.33   # Min probability to stay in regime
remain_upper <- 0.66   # Max probability to stay in regime
```

**M (Number of Regimes)**
- Type: Integer or vector
- Range: Positive integers (typically 2-4)
- Purpose: Number of distinct network states
- Impact: More regimes = more complex dynamics but harder to estimate

**remain_lower & remain_upper**
- Type: Numeric
- Range: 0 to 1, where `remain_lower < remain_upper`
- Purpose: Controls regime stability (how long networks stay in one regime)
- Higher values = More stable regimes
- Examples:
  - `[0.2, 0.4]` = Fast switching between regimes
  - `[0.5, 0.7]` = Moderate stability
  - `[0.8, 0.95]` = Very stable regimes

### 3. Edge Weight Parameters

```r
min_edg_val <- 0.05  # Minimum edge weight
max_edg_val <- 1     # Maximum edge weight
```

**min_edg_val (Edge Threshold)**
- Type: Numeric
- Range: 0 to max_edg_val
- Purpose: Edges below this are set to zero (sparse network enforcement)
- Common values: 0.05, 0.08, 0.1

**max_edg_val**
- Type: Numeric
- Range: Greater than min_edg_val
- Purpose: Maximum absolute edge weight
- Standard: 1 (normalized)

### 4. Time Series Parameters

```r
T <- c(3500)    # Time steps
warmup <- 50    # Warmup period
n_ts <- 30      # Number of time series
```

**T (Time Steps)**
- Type: Integer or vector
- Range: Positive integers (typically 1000-10000)
- Purpose: Length of each time series
- Impact:
  - Longer series = Better estimation but slower
  - Shorter series = Faster but less accurate
- Multiple values: `c(1000, 3500, 5000)` for comparing lengths

**warmup**
- Type: Integer
- Range: Positive integer (typically 50-200)
- Purpose: Initial steps discarded to allow dynamics to stabilize
- Rule of thumb: 1-5% of T

**n_ts (Number of Time Series)**
- Type: Integer
- Range: Positive integer (typically 10-100)
- Purpose: Replications per condition
- Impact: More series = Better statistics but much slower

### 5. Model Estimation Parameters

```r
order <- 1         # AR order
MaxIter <- 200     # Maximum EM iterations
eps <- 1e-5        # Convergence criterion
retry_attempts <- 5
```

**order (Autoregressive Order)**
- Type: Integer
- Range: Positive integer (usually 1)
- Purpose: Lag order for VAR model
- Standard: 1 (first-order VAR)

**MaxIter (Maximum Iterations)**
- Type: Integer
- Range: Positive integer (typically 100-500)
- Purpose: Maximum iterations for EM algorithm
- Impact:
  - Higher = More likely to converge but slower
  - Too low = Risk of non-convergence
- Recommended: 200

**eps (Convergence Epsilon)**
- Type: Numeric
- Range: Small positive number (typically 1e-5 to 1e-7)
- Purpose: Threshold for declaring convergence
- Smaller = Stricter convergence (slower)
- Standard: 1e-5

**retry_attempts**
- Type: Integer
- Range: Positive integer
- Purpose: Number of retries if fitting fails
- Recommended: 3-5

### 6. Output Parameters

```r
verbose <- FALSE      # Progress messages
save_output <- TRUE   # Save results
output_dir <- "output"
```

**verbose**
- Type: Logical
- Values: `TRUE` or `FALSE`
- Purpose: Print detailed progress messages
- Use `TRUE` for debugging

**save_output**
- Type: Logical
- Values: `TRUE` or `FALSE`
- Purpose: Whether to save results to files
- Recommended: `TRUE`

**output_dir**
- Type: String
- Purpose: Directory for saving results
- Default: "output"
- Creates directory if it doesn't exist

### 7. Reproducibility

```r
# set.seed(83742)  # Uncomment for reproducibility
```

**Random Seed**
- Uncomment and set seed for reproducible results
- Use any integer
- Results will be identical across runs with same seed

## Configuration Examples

### Small Quick Test

```r
N <- c(3)
Density <- c(0.3)
M <- c(2)
T <- c(500)
n_ts <- 5
MaxIter <- 100
```

### Standard Analysis

```r
N <- c(4)
Density <- c(0.25)
M <- c(2)
T <- c(3500)
n_ts <- 30
MaxIter <- 200
```

### Comprehensive Factorial Design

```r
N <- c(3, 4, 5)
Density <- c(0.2, 0.5, 0.8)
M <- c(2, 3, 4)
T <- c(1000, 3500, 5000)
n_ts <- 50
MaxIter <- 300
```

⚠️ Warning: This will generate 3 × 3 × 3 × 3 × 50 = 4,050 simulations!

## Tips for Configuration

### Starting Out
1. Use small values for testing (N=3, T=500, n_ts=5)
2. Check output looks reasonable
3. Scale up gradually

### Performance Tuning
- **Bottleneck**: n_ts and T have biggest impact on runtime
- **Memory**: Large N × T × n_ts can exhaust RAM
- **Parallel**: Consider parallelizing across n_ts (future work)

### Estimation Quality
- Minimum recommended: T=1000, n_ts=20
- Good quality: T=3500, n_ts=30
- High quality: T=5000, n_ts=50+

### Troubleshooting
- **Non-convergence**: Increase MaxIter or decrease eps
- **Poor recovery**: Increase T or n_ts
- **Slow runtime**: Reduce n_ts or T
- **Memory issues**: Reduce N or run smaller batches

## Output Files

When `save_output = TRUE`, creates timestamped files:

```
output/
├── Timeseries_data_20250118_120530.rds
├── MSAR_models_20250118_120530.rds
└── CONFIG_20250118_120530.rds
```

**Timeseries_data**: Generated network time series and true dynamics
**MSAR_models**: Estimated models and comparison statistics
**CONFIG**: Configuration used for this run (for reproducibility)

## Best Practices

1. **Document changes**: Keep notes on which parameters worked well
2. **Version control**: Commit configuration before long runs
3. **Save config**: The CONFIG.rds file captures your exact settings
4. **Start small**: Test with reduced parameters before full run
5. **Monitor progress**: Use `verbose = TRUE` for long runs
6. **Check output**: Verify configuration summary at start

## Advanced Usage

### Parameter Sweeps

Run the script multiple times with different configurations:

```r
# Script 1: Vary network size
N <- c(3, 4, 5)
M <- c(2)
# ... run

# Script 2: Vary regimes
N <- c(4)
M <- c(2, 3, 4)
# ... run
```

### Batch Processing

Create wrapper script for multiple runs:

```r
conditions <- expand.grid(
  N = c(3, 4, 5),
  M = c(2, 3)
)

for (i in 1:nrow(conditions)) {
  # Set parameters
  N <- conditions$N[i]
  M <- conditions$M[i]

  # Source main script
  source("scripts/MSAR_ts_analysis_NHMSAR.R")
}
```
