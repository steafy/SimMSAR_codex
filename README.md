# SimMSAR - Markov-Switching Autoregressive Models for Network Dynamics

This project implements simulation and estimation of Markov-Switching Autoregressive (MSAR) models for network dynamics using the NHMSAR package.

## Project Structure

```
SimMSAR_claude/
├── R/                          # R source code
│   ├── estimation/             # MSAR model estimation functions
│   │   ├── estimate_MSAR.R            # Main estimation pipeline
│   │   ├── fit.MSAR_revised_2.R       # Model fitting with EM algorithm
│   │   ├── init.theta.MSAR_revised_2.R # Parameter initialization
│   │   ├── init_and_fit.MSAR_Lasso_2.R # Combined init & fit with retry
│   │   ├── as.thetaMSAR_revised_2.R   # Parameter object conversion
│   │   ├── Mstep.hh.lasso.MSAR_patched_2.R  # M-step with LASSO
│   │   ├── Mstep.hh.reduct.MSAR_patched_2.R # M-step with reduction
│   │   └── EM_converged_patched_2.R   # Convergence checking
│   ├── generation/             # Data generation functions
│   │   ├── generate_timeseries.R      # Main timeseries generation
│   │   ├── generate_netdyn.R          # Network dynamics generation
│   │   ├── generate_Beta.R            # Temporal dynamics (Beta matrix)
│   │   ├── generate_kappa.R           # Precision matrix (kappa)
│   │   ├── generate_transmat.R        # Transition matrix
│   │   ├── generate_random.R          # Random network generation
│   │   └── check_stability.R          # Stability checking
│   ├── utils/                  # Utility functions
│   │   ├── asign_regimes.R            # Regime assignment
│   │   ├── senspec.R                  # Sensitivity/specificity
│   │   ├── summarize_cor.R            # Correlation summary
│   │   └── calculate_MAE.R            # Mean absolute error
│   ├── analysis/               # Statistical analysis
│   │   └── stat_analysis.R            # Post-hoc statistical analysis
│   └── visualization/          # Plotting functions
│       └── Create_3D_surfaceplots_for_means.R
├── scripts/                    # Main execution scripts
│   └── MSAR_ts_analysis_NHMSAR.R      # Main simulation & analysis script
├── Data/                       # Data files (gitignored)
├── archive/                    # Archived/deprecated code
│   ├── old_versions/          # Superseded version 1 files
│   ├── obsolete/              # Broken/non-functional scripts
│   └── experimental/          # Experimental implementations
├── .gitignore
└── README.md
```

## Workflow

### 1. Data Generation
The workflow starts by generating synthetic network dynamics and time-series data:

```r
source("R/generation/generate_timeseries.R")

Timeseries_data <- generate_timeseries(
  Density = c(0.25),      # Network edge density
  N = c(4),               # Number of nodes
  M = c(2),               # Number of regimes
  T = c(3500),            # Time steps
  n_ts = 30,              # Number of time series
  warmup = 50,            # Warmup period
  min_edg_val = 0.05,     # Min edge weight
  max_edg_val = 1,        # Max edge weight
  remain_lower = 0.33,    # Min prob. to stay in regime
  remain_upper = 0.66     # Max prob. to stay in regime
)
```

### 2. Model Estimation
Estimate MSAR models from the generated time-series:

```r
source("R/estimation/estimate_MSAR.R")

MSAR_models <- estimate_MSAR(
  Density = Density,
  N = N,
  M = M,
  T = T,
  n_ts = n_ts,
  order = 1,              # AR order
  MaxIter = 200,          # Max EM iterations
  verbose = FALSE,
  min_edg_val = 0.05,
  Timeseries_data = Timeseries_data
)
```

### 3. Statistical Analysis
Analyze the estimation results:

```r
source("R/analysis/stat_analysis.R")
# Performs PERMANOVA, linear mixed models, and generates plots
```

## Quick Start

Run the main analysis script:

```r
source("scripts/MSAR_ts_analysis_NHMSAR.R")
```

This will:
1. Generate synthetic time-series data with specified network parameters
2. Estimate MSAR models for each time series
3. Compare estimated vs. true network dynamics
4. Save results to `.rds` files

## Dependencies

Required R packages:
- `NHMSAR` - Core MSAR functionality
- `netcontrol` - Network control analysis
- `huge` - Nonparanormal transformation
- `dplyr` - Data manipulation
- `ggplot2` - Visualization
- `lme4`, `lmerTest` - Mixed models
- `lmPerm` - Permutation tests
- `progress` - Progress bars
- `mvtnorm` - Multivariate normal distributions
- `graphicalVAR` - Graphical VAR models
- Additional: `dunn.test`, `sjPlot`, `knitr`, `kableExtra`, `xtable`, `cowplot`

## Output

Results are saved as:
- `Data/Timeseries_data_*.rds` - Generated time-series data
- `Data/MSAR_models_*.rds` - Estimated MSAR models
- `*.rds` files in root - Analysis results
- Plots in designated output directory

## Notes

- All source paths are relative to the project root directory
- The `_2` suffix on estimation files indicates version 2 (current/active version)
- Version 1 files are archived in `archive/old_versions/`
- Large data files (`.rds`) are gitignored to keep repository size manageable

## License

[Add your license information here]

## Contact

[Add your contact information here]
