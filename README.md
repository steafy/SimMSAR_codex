# SimMSAR

Simulation study on parameter and controllability recovery in Markov-Switching
Vector Autoregressive (MS-VAR) models of network dynamics. Synthetic multivariate
time series are generated from known regime-specific networks, Markov-Switching
Autoregressive (MSAR) models are fit back to them, and the recovery of the true
network structure and its average controllability is scored across a factorial
design (nodes x density x regimes x timesteps). This is the companion code to the
master's thesis *Network Controllability in Time-Varying Dynamical Systems*
(Philipps-Universität Marburg, Department of Psychology, 2026).

## Notation

The code uses the same notation as the manuscript.

| Symbol | Meaning |
|--------|---------|
| `A`      | Temporal network — the VAR(1) coefficient matrix (directed, lag-1) |
| `K`      | Contemporaneous network — the precision matrix of the residuals |
| `Sigma`  | Residual covariance; `Sigma = K^{-1}` |
| `AC`     | Average controllability of the temporal network `A` |
| `M`      | Number of regimes (network states) |
| `N`      | Number of nodes |
| `T`      | Number of timesteps per series |

Prefixes `orig_` / `est_` denote the true (generating) and estimated quantities
(e.g. `orig_A`, `est_K`). Recovery outcomes follow the same scheme: `A_corr`,
`K_corr`, `AC_corr_pearson` (correlation of true vs. estimated), `A_sen`/`A_spec`
and `K_sen`/`K_spec` (edge sensitivity/specificity), `NRMSE_A`/`NRMSE_K`.


## Installation

Requires **R >= 4.0.0**.

Most dependencies are on CRAN:

```r
install.packages(c(
  "glmnet", "lars", "mvtnorm", "abind", "dplyr", "tibble", "stringr",
  "future", "future.apply", "progressr", "lme4", "lmerTest", "emmeans",
  "performance", "dunn.test", "ggplot2", "cowplot", "viridisLite",
  "systemfonts", "kableExtra", "xtable", "clue"
))
```

**NHMSAR** is no longer on the current CRAN and must be installed from the CRAN
archive:

```r
install.packages("remotes")
remotes::install_url(
  "https://cran.r-project.org/src/contrib/Archive/NHMSAR/NHMSAR_1.19.tar.gz"
)
```

Sourcing `R/dependencies.R` checks for all required packages and, if any are
missing, prints the exact install command and stops (it never installs
automatically).

Optional: `elementalist` (GitHub-only) enables rounded facet-card styling in the
plots; the plotting code falls back gracefully if it is absent.

```r
remotes::install_github("teunbrand/elementalist")
```

Building the HTML result tables (`kableExtra::save_kable`) additionally needs
**pandoc** on the system path; the LaTeX (`.tex`) exports do not.

## Reproducing the pipeline

Run everything from the repository root (scripts use relative paths).

**(a) Simulation + estimation** — generates the synthetic series, fits the MSAR
models, and writes the raw estimates:

```r
source("scripts/MSAR_ts_analysis_NHMSAR.R")
```

This writes timestamped `output/MSAR_models_<timestamp>.rds`,
`Timeseries_data_<timestamp>.rds` and `CONFIG_<timestamp>.rds`. Expect a long
runtime: the full design is 3 x 3 x 4 x 4 x 150 replications = **21,600 fits**
(several hours on a workstation).

**(b) Analysis** — scores recovery and regenerates every table and plot under
`output/` from a saved estimation run:

```r
MSAR_dynamics_list <- readRDS("output/MSAR_models_final.rds")
source("R/analysis/stat_analysis.R")
```

All scoring (edge thresholding, `K = solve(Sigma)`, correlations, NRMSE,
sensitivity/specificity, average controllability) happens in the analysis step,
so scoring decisions can be revised without re-running the expensive estimation.

## Design factors

| Factor | Levels |
|--------|--------|
| Nodes (`N`)      | 4, 6, 8 |
| Density          | 0.25, 0.50, 0.75 |
| Regimes (`M`)    | 1, 2, 3, 4 |
| Timesteps (`T`)  | 200, 400, 800, 1600 |
| Replications     | 150 per cell |

## Repository structure

```
R/
  dependencies.R              Package checks (fail-fast; no auto-install)
  generation/                 Data generation
    generate_A.R              Temporal network A (stable VAR coefficients)
    generate_K.R              Contemporaneous precision matrix K
    generate_netdyn.R         Assemble one regime's dynamics (A, K, Sigma, AC)
    generate_transmat.R       Markov transition matrix
    generate_random.R         Random signed edge weights
    generate_timeseries.R     Simulate the MSAR series across the design grid
    check_stability.R         VAR stability check
  estimation/                 MSAR fitting (NHMSAR-derived, see attribution)
    fit_msar.R                Core EM fit
    init_theta_msar.R         Parameter initialisation
    as_theta_msar.R           thetaMSAR coercion
    em_converged.R            EM convergence test
    mstep_hh_lasso_msar.R     cv.glmnet LASSO M-step (opt-in engine)
    mstep_hh_lasso_msar_bic.R lars + BIC M-step (default engine)
    mstep_hh_reduct_msar.R    Reduction (support-copy) M-step
    init_and_fit_msar_lasso.R Init + fit wrapper with retries
    estimate_MSAR.R           Parallel estimation driver (raw estimates only)
    match_regimes.R           Match estimated to true regimes (Hungarian)
    reconstruct_regime_sequence.R  Relabel the estimated regime sequence
    stabilize_sigma.R         Opt-in Sigma eigenvalue stabilisation
  utils/
    average_controllability.R Average controllability of A
    senspec.R / undirected_metrics.R  Edge-recovery metrics
    assign_regimes.R          Hungarian assignment helper
    regime_sequence_recovery.R  Hard sequence + Cohen's kappa
  analysis/                   Scoring, modelling, tables, plots
    stat_analysis.R           Analysis entry point (source after loading the .rds)
    data_prep.R transform.R descriptives.R modeling.R sensitivity.R
    sequence_models.R feasibility.R exports.R plots.R reporting.R output_paths.R
scripts/
  MSAR_ts_analysis_NHMSAR.R   Simulation + estimation driver
  migrate_rds_names.R         One-off Beta/Kappa -> A/K rename for old .rds
docs/                         Methodological notes (M-step, Sigma/K degeneracy)
output/                       Final figures (plots/) and tables (results/); the
                              published dataset (.rds) is added via Git LFS
```

## Configuration

The main knobs are set at the top of `scripts/MSAR_ts_analysis_NHMSAR.R`; the
LASSO engine is selected via `options()` (applied per worker via
`lasso_control`). See `docs/MSTEP_LASSO_CV_PENALIZATION.md` for the engine
comparison.

| Setting | Default | Meaning |
|---------|---------|---------|
| `simmsar_lasso_engine`   | `"bic"`   | First M-step engine: `"bic"` (lars + BIC) or `"cvglmnet"` (cross-validated LASSO) |
| `simmsar_lasso_reselect` | `FALSE`   | `cvglmnet` only: re-select the LASSO support every EM iteration instead of freezing iteration 1 |
| `min_edg_val` (`MIN_EDG_VAL`) | `0.05` | Edge-detection threshold; must match the generation intent |
| `AC_HORIZON`             | `25`      | Horizon `T_ac` for average controllability; must equal the generation default |
| `MaxIter`                | `200`     | Maximum EM iterations per fit |
| `workers`                | `5`       | Parallel worker processes (`NULL` = physical cores − 1) |

## Reproducibility

The driver sets `set.seed(20260701)`. Parallel sections use `future.seed = TRUE`,
which derives an L'Ecuyer-CMRG RNG stream from the current seed: the same script,
seed and worker count reproduce results run-to-run.

Note, however, that exact numeric values differ from any prior **serial** run
using the same `set.seed()`. Parallel workers consume random draws in a different
order, and generation additionally draws residuals via a once-per-cell Cholesky
factorisation. This is expected and statistically immaterial — it changes which
particular replicates land where, not the design or the conclusions.

## Known failure modes

- **Non-convergence / estimation failure** is concentrated in the hardest cells
  (short series `T = 200`, high density, `N = 8`, `M = 4`), where the EM step can
  abort or every estimated regime collapses to zero variance. Such fits are
  logged (`summarize_failures()`) and excluded; `R/analysis/sensitivity.R`
  checks that this missingness does not bias the reported estimates.
- **Near-singular `est_Sigma`.** The estimated residual covariance can be
  ill-conditioned; inversion to `K` is gated on its condition number
  (`K_COND_MAX`) and flagged invalid rather than propagated. An opt-in
  eigenvalue stabilisation is available (`simmsar_sigma_stab`, default off). See
  `docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md`.
- **`K_corr` is `NaN`** when the thresholded estimated `K` has no surviving
  off-diagonal edge in any regime (Pearson correlation is undefined for a
  constant vector). This is a scoring artifact of an empty estimated network, not
  an estimation crash.

## Published dataset

`output/MSAR_models_final.rds` (with `output/CONFIG_final.rds`) is the final run
reported in the thesis: 21,591 of the 21,600 expected fits completed
successfully. It is tracked via Git LFS. Load it and source
`R/analysis/stat_analysis.R` to regenerate every table and figure under
`output/`.

## License & attribution

This project is released under the **GNU General Public License v3.0** (see
`LICENSE`).

The estimation core (`R/estimation/fit_msar.R`, `init_theta_msar.R`,
`as_theta_msar.R`, `em_converged.R`, `mstep_hh_lasso_msar.R`,
`mstep_hh_lasso_msar_bic.R`, `mstep_hh_reduct_msar.R`) is adapted from the
**NHMSAR** package by Valerie Monbet (obtained from the CRAN archive; original
license GPL). The adaptations add LASSO / reduction M-step engines and opt-in
Sigma stabilisation; each file carries an attribution header.

If you use this code, please cite the thesis: *Becker, S. (2026). Network
Controllability in Time-Varying Dynamical Systems. Philipps-Universität Marburg, Department of Psychology.*
