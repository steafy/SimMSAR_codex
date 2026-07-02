# =============================================================================
# DIAGNOSTIC (uncommitted): where do degenerate Sigma / Kappa estimates come from?
# =============================================================================
# Companion to docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md. Runs a pilot (n_ts=5)
# across the full T x Density x Nodes x Regimes grid with the default
# cvglmnet + reselect engine, captures per-regime postmix / separation-quality /
# Sigma-conditioning diagnostics (added to estimate_MSAR()), then builds a
# side-table that joins each regime's Sigma/Kappa degeneracy flags to its
# effective sample size, to test the "low postmix drives degeneracy" hypothesis.
#
# NOTHING here is committed. Outputs go to output/diagnostics/ (untracked).
# Run from the project root:
#   Rscript scripts/diagnose_sigma_kappa.R
# =============================================================================

suppressMessages({
  source("R/dependencies.R")
  source("R/generation/generate_timeseries.R")
  source("R/estimation/estimate_MSAR.R")
  source("R/analysis/data_prep.R")            # compute_recovery_metrics + summaries
  source("R/utils/undirected_metrics.R")      # symmetrize / vectorize / senspec_upper_tri
})

set.seed(58396)   # same master seed as the production driver

# ---- pilot grid: same FACTORS as prior pilots ------------------------------
N        <- c(4, 6, 8)
Density  <- c(0.25, 0.5, 0.75)
M        <- c(1, 2, 3, 4)
T        <- c(200, 400, 800, 1600)
n_ts     <- 5
warmup   <- 50
totTime  <- T + warmup
order    <- 1
MaxIter  <- 200
eps      <- 1e-5
min_edg_val <- 0.05
max_edg_val <- 1
remain_lower <- 0.85
remain_upper <- 0.85
workers  <- 5

# Analysis-side gates (must mirror the production analysis defaults) ----------
KAPPA_COND_MAX          <- 1e6
max_plausible_magnitude <- 10

# ---- LASSO engine: production default (cvglmnet + reselect) ----------------
lasso_control <- list(
  simmsar_lasso_engine         = "cvglmnet",
  simmsar_lasso_reselect       = TRUE,
  simmsar_lasso_lambda         = "1se",
  simmsar_lasso_refit          = TRUE,
  simmsar_lasso_adaptive       = FALSE,
  simmsar_lasso_nfolds         = 5,
  simmsar_lasso_nlambda        = 50,
  simmsar_lasso_fixedfolds     = TRUE,
  simmsar_lasso_reselect_iters = Inf,
  simmsar_lasso_reselect_every = 1
)
do.call(options, lasso_control)

out_dir <- "output/diagnostics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("=== DIAGNOSTIC PILOT: data generation ===\n")
t0 <- Sys.time()
Timeseries_data <- generate_timeseries(
  Density = Density, min_edg_val = min_edg_val, max_edg_val = max_edg_val,
  M = M, N = N, warmup = warmup, T = T, totTime = totTime,
  n_ts = n_ts, remain_lower = remain_lower, remain_upper = remain_upper,
  workers = workers
)
cat(sprintf("generation done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("\n=== DIAGNOSTIC PILOT: estimation (cvglmnet + reselect) ===\n")
t1 <- Sys.time()
MSAR <- estimate_MSAR(
  Density = Density, M = M, N = N, T = T, n_ts = n_ts, order = order,
  MaxIter = MaxIter, verbose = FALSE, min_edg_val = min_edg_val,
  Timeseries_data = Timeseries_data, workers = workers,
  lasso_control = lasso_control
)
cat(sprintf("estimation done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t1, units = "mins"))))

# Score with exactly the production analysis gates so the degeneracy flags match
# what the pipeline would report.
MSAR <- compute_recovery_metrics(
  MSAR, min_edg_val = min_edg_val, AC_HORIZON = 25,
  KAPPA_COND_MAX = KAPPA_COND_MAX, max_plausible_magnitude = max_plausible_magnitude
)

saveRDS(MSAR, file.path(out_dir, sprintf("MSAR_diag_%s.rds", stamp)))

# ---------------------------------------------------------------------------
# Build the per-regime degeneracy side-table (self-contained: invert est_Sigma
# here so we characterise degeneracy even for rows the pipeline NA'd out).
# ---------------------------------------------------------------------------
n <- nrow(MSAR)
diag_tbl <- MSAR[, c("timesteps","density","nodes","regimes","ts_id","regime_id",
                     "postmix","postmix_frac","est_Sigma_rcond","est_Sigma_min_eig",
                     "gamma_entropy_mean","iter",
                     "Kappa_corr","Kappa_sen","Kappa_spec","NRMSE_Kappa")]
diag_tbl <- as.data.frame(diag_tbl)

cond_number    <- numeric(n)
solve_failed   <- logical(n)
kappa_offdiag_max <- numeric(n)   # max |offdiag| of inverted (raw) Kappa
kappa_dense    <- logical(n)      # thresholded Kappa has specificity ~0 (all offdiag kept)
kappa_empty    <- logical(n)      # thresholded Kappa has NO offdiag edges (corr undefined)
sigma_trace_frac_lt <- diag_tbl$postmix_frac  # convenience

for (r in seq_len(n)) {
  S <- as.matrix(MSAR$est_Sigma[[r]])
  rc <- tryCatch(rcond(S), error = function(e) 0)
  cond_number[r] <- if (is.finite(rc) && rc > 0) 1 / rc else Inf
  inv <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(inv)) {
    solve_failed[r] <- TRUE
    kappa_offdiag_max[r] <- NA_real_
    kappa_dense[r] <- NA; kappa_empty[r] <- NA
    next
  }
  K <- symmetrize_matrix(inv)
  off <- !diag(TRUE, nrow(K))
  kappa_offdiag_max[r] <- max(abs(K[off]))
  # thresholded off-diagonal (mirror threshold_kappa_offdiag)
  Kt <- K; Kt[off & (abs(K) < min_edg_val)] <- 0
  n_off_edges <- sum(Kt[off] != 0) / 2         # symmetric
  kappa_empty[r] <- (n_off_edges == 0)
  # specificity vs truth (share of true non-edges left at 0); dense if ~0
  ss <- tryCatch(senspec_upper_tri(as.matrix(MSAR$orig_Kappa[[r]]), Kt, diag = FALSE),
                 error = function(e) list(specificity = NA_real_))
  kappa_dense[r] <- isTRUE(ss[["specificity"]] < 0.05)
}

diag_tbl$cond_number       <- cond_number
diag_tbl$solve_failed      <- solve_failed
diag_tbl$kappa_offdiag_max <- kappa_offdiag_max
diag_tbl$kappa_dense       <- kappa_dense
diag_tbl$kappa_empty       <- kappa_empty
diag_tbl$ill_conditioned   <- is.finite(cond_number) & (cond_number > KAPPA_COND_MAX) | !is.finite(cond_number)
diag_tbl$mag_flagged       <- is.finite(kappa_offdiag_max) & (kappa_offdiag_max > max_plausible_magnitude)
diag_tbl$degenerate        <- diag_tbl$solve_failed | diag_tbl$ill_conditioned |
                              diag_tbl$mag_flagged | (diag_tbl$kappa_empty %in% TRUE)

write.csv(diag_tbl, file.path(out_dir, sprintf("sigma_kappa_diag_%s.csv", stamp)),
          row.names = FALSE)
saveRDS(diag_tbl, file.path(out_dir, sprintf("sigma_kappa_diag_%s.rds", stamp)))

# ---------------------------------------------------------------------------
# Headline analysis printed to stdout (captured into the run log).
# ---------------------------------------------------------------------------
cat("\n\n############ DEGENERACY DIAGNOSIS SUMMARY ############\n")
cat(sprintf("Total regime rows: %d  (fits: %d)\n", n,
            nrow(unique(diag_tbl[, c("timesteps","density","nodes","regimes","ts_id")]))))

flag_counts <- c(
  solve_failed    = sum(diag_tbl$solve_failed, na.rm = TRUE),
  ill_conditioned = sum(diag_tbl$ill_conditioned, na.rm = TRUE),
  mag_flagged     = sum(diag_tbl$mag_flagged, na.rm = TRUE),
  kappa_empty     = sum(diag_tbl$kappa_empty %in% TRUE),
  kappa_dense     = sum(diag_tbl$kappa_dense %in% TRUE),
  degenerate_any  = sum(diag_tbl$degenerate, na.rm = TRUE)
)
cat("\n-- flag counts (per regime row) --\n"); print(flag_counts)
cat("\n-- flag rates (%) --\n"); print(round(100 * flag_counts / n, 2))

cat("\n-- postmix / separation by degeneracy status --\n")
agg <- function(mask, lab) {
  d <- diag_tbl[mask & is.finite(diag_tbl$postmix), ]
  cat(sprintf("%-16s n=%4d  postmix: median=%.1f  min=%.1f  |  postmix_frac: median=%.3f  min=%.3f  |  entropy median=%.3f  |  cond median=%.2e\n",
              lab, nrow(d),
              median(d$postmix), min(d$postmix),
              median(d$postmix_frac), min(d$postmix_frac),
              median(d$gamma_entropy_mean),
              median(d$cond_number[is.finite(d$cond_number)])))
}
agg(diag_tbl$degenerate,  "DEGENERATE")
agg(!diag_tbl$degenerate, "healthy")

cat("\n-- degeneracy rate by regimes M --\n")
print(round(100 * tapply(diag_tbl$degenerate, diag_tbl$regimes, mean), 1))
cat("\n-- degeneracy rate by timesteps T --\n")
print(round(100 * tapply(diag_tbl$degenerate, diag_tbl$timesteps, mean), 1))
cat("\n-- degeneracy rate by nodes N --\n")
print(round(100 * tapply(diag_tbl$degenerate, diag_tbl$nodes, mean), 1))
cat("\n-- degeneracy rate by density --\n")
print(round(100 * tapply(diag_tbl$degenerate, diag_tbl$density, mean), 1))

cat("\n-- degeneracy rate by low-postmix quintile (postmix_frac) --\n")
q <- cut(diag_tbl$postmix_frac, breaks = quantile(diag_tbl$postmix_frac, 0:5/5, na.rm = TRUE),
         include.lowest = TRUE)
print(round(100 * tapply(diag_tbl$degenerate, q, mean), 1))

cat("\n-- postmix_frac vs degeneracy: logistic fit --\n")
fit_glm <- tryCatch(
  glm(degenerate ~ postmix_frac + gamma_entropy_mean + factor(nodes) + factor(regimes),
      data = diag_tbl, family = binomial()),
  error = function(e) NULL)
if (!is.null(fit_glm)) print(summary(fit_glm)$coefficients)

# Cross-check: does the "fully empty across ALL regimes of a fit" pattern occur,
# and does it coincide with low postmix?
cat("\n-- fits where EVERY regime's Kappa is empty (all-empty pattern) --\n")
by_fit <- split(diag_tbl, interaction(diag_tbl$timesteps, diag_tbl$density,
                                      diag_tbl$nodes, diag_tbl$regimes, diag_tbl$ts_id, drop = TRUE))
all_empty <- vapply(by_fit, function(f) all(f$kappa_empty %in% TRUE), logical(1))
cat(sprintf("all-empty fits: %d of %d\n", sum(all_empty), length(by_fit)))
if (any(all_empty)) {
  ae <- do.call(rbind, by_fit[all_empty])
  cat("their design cells (T,D,N,M):\n")
  print(unique(ae[, c("timesteps","density","nodes","regimes")]))
  cat(sprintf("their postmix_frac: median=%.3f min=%.3f (vs overall median %.3f)\n",
              median(ae$postmix_frac), min(ae$postmix_frac), median(diag_tbl$postmix_frac)))
}

cat("\nSaved:\n")
cat("  ", file.path(out_dir, sprintf("MSAR_diag_%s.rds", stamp)), "\n")
cat("  ", file.path(out_dir, sprintf("sigma_kappa_diag_%s.csv", stamp)), "\n")
cat("  ", file.path(out_dir, sprintf("sigma_kappa_diag_%s.rds", stamp)), "\n")
cat("\nDONE (stamp ", stamp, ")\n", sep = "")
