#' Estimate MSAR Models from Time Series Data (raw estimates only, parallel)
#'
#' Estimation pipeline that fits Markov-Switching Autoregressive (MSAR) models to
#' each simulated time series and returns the \emph{raw, uninterpreted} estimates
#' together with the regime matching / sequence reconstruction that has to happen
#' while the EM fit object is still in memory. It deliberately does NOT compute
#' any derived outcome metric (correlations, NRMSE, sensitivity/specificity,
#' Kappa inversion, average controllability, regime-sequence accuracy): all of
#' that now lives in the analysis pipeline (\code{R/analysis/}) and is computed
#' from the stored matrices, so every scoring decision is revisable without
#' re-running this expensive step.
#'
#' @param Density Numeric vector. Edge densities used in data generation.
#' @param N Integer vector. Numbers of nodes in networks.
#' @param M Integer vector. Numbers of regimes.
#' @param T Integer vector. Time series lengths.
#' @param n_ts Integer. Number of time series per condition.
#' @param order Integer. AR order for the VAR model (typically 1).
#' @param MaxIter Integer. Maximum EM algorithm iterations.
#' @param verbose Logical. If TRUE, prints detailed progress messages. Note:
#'   with parallel workers, message()/cat() output from individual fits may be
#'   interleaved/delayed; this only affects what you see on screen, not results.
#' @param min_edg_val Numeric. Retained for call-site compatibility but
#'   \strong{no longer used here}: edge thresholding is a scoring decision and
#'   now lives in the analysis pipeline (\code{compute_recovery_metrics()}),
#'   which applies the threshold to the stored raw \code{est_Beta} / inverted
#'   \code{est_Kappa}. Pass the intended threshold there instead.
#' @param Timeseries_data Tibble. Output from \code{generate_timeseries} containing
#'   simulated data and true dynamics in tibble format with list-columns.
#' @param workers Integer. Number of parallel worker processes (via
#'   \code{future::multisession}, which works identically on Windows/Mac/Linux).
#'   Defaults to \code{max(1, parallel::detectCores(logical = FALSE) - 1)}, i.e.
#'   one less than the number of PHYSICAL cores -- this is a conservative
#'   starting point for compute-bound numerical work (EM + LASSO), where
#'   hyperthreading on logical cores typically buys only a modest amount of
#'   extra throughput. Pass \code{workers = 1} to fall back to fully sequential
#'   execution (e.g. for debugging or a sanity-check comparison run).
#' @param lasso_control Optional named list of \code{simmsar_lasso_*} options
#'   controlling the LASSO first-M-step engine (see
#'   \code{docs/MSTEP_LASSO_CV_PENALIZATION.md} and
#'   \code{R/estimation/mstep_hh_lasso_msar.R}). It is applied via
#'   \code{options()} \emph{inside each worker} (options do not propagate to
#'   future workers), so it is the correct way to configure the penalization for
#'   a parallel run. \code{NULL} (default) leaves the library defaults in force
#'   (\code{engine="bic"}, i.e. the legacy lars+BIC step). Build it in the
#'   calling script, e.g. \code{list(simmsar_lasso_engine = "cvglmnet",
#'   simmsar_lasso_reselect = TRUE)}.
#'
#' @return An \code{msar_results} tibble with \strong{one row per regime per
#'   successfully fitted time series} and only the raw quantities:
#'   \describe{
#'     \item{timesteps, density, nodes, regimes, ts_id, regime_id}{Condition / regime identifiers}
#'     \item{orig_Beta, est_Beta}{True / estimated temporal network (raw, NOT thresholded)}
#'     \item{orig_Kappa}{True contemporaneous precision matrix (pass-through from generation)}
#'     \item{orig_Sigma, est_Sigma}{True / estimated residual covariance. \code{est_Sigma}
#'       is stored as-is, however well- or ill-conditioned; it is NOT inverted here.}
#'     \item{orig_Beta_ac}{True Beta average controllability (pass-through from generation)}
#'   }
#'   \code{est_Kappa}, \code{est_Beta_ac} and all derived metrics are produced
#'   downstream in analysis, not here.
#'
#'   A fit-level table (one row per \code{ts_id}/condition, NOT per regime) is
#'   attached via \code{attr(result, "fit_results")} with columns
#'   \code{timesteps, density, nodes, regimes, ts_id, orig_TPM, est_TPM,
#'   orig_regime_sequence, est_regime_sequence}. The estimated TPM and sequence
#'   are relabelled with the same Beta-based matching permutation as the
#'   per-regime table, so regime \code{m} means the same thing across all stored
#'   objects. \code{(density, nodes, regimes, ts_id)} uniquely identifies the
#'   reused generating parameter set, so analysis can build the random-effects
#'   grouping factor from it.
#'
#'   A diagnostic log of discarded fits (reason + design cell) is attached via
#'   \code{attr(result, "failure_log")}; summarise it with
#'   \code{\link{summarize_failures}}.
#'
#' @details
#' Per time series: (1) fit via \code{init_and_fit_msar_lasso} (EM + LASSO, with
#' its own retry logic, unchanged); (2) match estimated regimes to true regimes
#' (Hungarian algorithm on vectorised Beta correlations, \code{\link{match_regimes}});
#' (3) hard-decode the smoothed probabilities into a regime sequence and relabel
#' it with the matching permutation (\code{\link{reconstruct_regime_sequence}});
#' (4) store the raw matrices. No edge thresholding, no covariance inversion, no
#' near-singular-Sigma gate -- those moved to analysis.
#'
#' PARALLELISATION (2026-06): every (T, Density, N, M, replication) combination
#' is a fully independent fit on its own simulated series -- embarrassingly
#' parallel. Implemented via \code{future.apply::future_lapply()} over
#' \code{future::multisession} workers (portable: identical on Windows/Mac/
#' Linux, unlike \code{future::multicore} which silently falls back to
#' sequential on Windows). Three points this design specifically addresses:
#' \itemize{
#'   \item \strong{Memory}: \code{Timeseries_data} is sliced into one minimal
#'     per-task payload BEFORE the parallel section (\code{build_task_data()}),
#'     and the big tibble is then dropped (\code{rm() + gc()}) from the calling
#'     process. \code{future_lapply} ships each worker only the slice of tasks
#'     it actually runs, not the whole input -- avoiding an N-times-duplicated
#'     copy of a potentially multi-GB object across N worker processes.
#'   \item \strong{Reproducibility}: \code{future.seed = TRUE} derives an
#'     L'Ecuyer-CMRG parallel RNG stream from the current \code{.Random.seed},
#'     giving fully reproducible results across repeated runs of the SAME
#'     script -- but the actual numeric values will differ from any prior
#'     SERIAL run with the same \code{set.seed()}, since draws are consumed in
#'     a different order. This is expected and statistically immaterial.
#'   \item \strong{Load balancing}: \code{future.scheduling = Inf} dispatches
#'     one task at a time rather than pre-chunking, because fit cost is highly
#'     heterogeneous across the design grid (e.g. N=8/T=4000/M=4 is vastly more
#'     expensive than N=4/T=250/M=1); per-task dispatch overhead is negligible
#'     next to multi-second fit times, so this trades a little overhead for
#'     much better balance across workers.
#' }
#' Packages are NOT part of \code{future}'s automatic global-export (only plain
#' R objects/closures are), so each worker re-attaches them via
#' \code{ensure_worker_packages()}, guarded to run only once per worker process
#' (multisession workers are a persistent pool reused across tasks, not
#' respawned per call).
#'
#' @note Dependencies are loaded centrally via R/dependencies.R.
#'
#' @seealso
#' \code{\link{generate_timeseries}}, \code{\link{init_and_fit_msar_lasso}},
#' \code{\link{match_regimes}}, \code{\link{reconstruct_regime_sequence}},
#' \code{\link{summarize_failures}}
#'
#' @export
# Dependencies are loaded centrally via R/dependencies.R

# Load functions from NHMSAR
source("R/estimation/fit_msar.R")
source("R/estimation/init_theta_msar.R")
source("R/estimation/init_and_fit_msar_lasso.R")

# Load regime matching / sequence reconstruction (deterministic; metrics live in analysis)
source("R/utils/assign_regimes.R")
source("R/utils/regime_sequence_recovery.R")
source("R/estimation/match_regimes.R")
source("R/estimation/reconstruct_regime_sequence.R")


# -----------------------------------------------------------------------------
# Re-attach packages inside a worker process. Plain R objects/functions defined
# via source() are auto-exported by future's globals-detection, but attached
# PACKAGES are process-level state and are never auto-exported -- every worker
# needs its own library() calls. Guarded via an option flag so a persistent
# multisession worker only pays this cost once, not on every one of thousands
# of tasks it processes over the life of the future::plan().
ensure_worker_packages <- function() {
  if (!isTRUE(getOption("simmsar_worker_pkgs_loaded"))) {
    if (exists("load_packages", mode = "function")) {
      load_packages(quietly = TRUE)
    } else {
      source("R/dependencies.R")
    }
    options(simmsar_worker_pkgs_loaded = TRUE)
  }
}

# -----------------------------------------------------------------------------
# Typed empty schemas, shared by the per-cell combine (fit_one_cell) and the
# final top-level combine (estimate_MSAR) -- both can hit the "nothing
# succeeded" edge case (a whole cell where every replicate failed; or, at the
# top level, every cell). Centralised here so both stay in sync.
empty_regime_rows_schema <- function() {
  tibble::tibble(
    timesteps = integer(0), density = numeric(0), nodes = integer(0),
    regimes = integer(0), ts_id = integer(0), regime_id = integer(0),
    orig_Beta = list(), est_Beta = list(), orig_Kappa = list(),
    orig_Sigma = list(), est_Sigma = list(), orig_Beta_ac = list(),
    # --- Sigma/Kappa degeneracy diagnostics (uncommitted; see
    #     docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md). Per (matched) regime:
    postmix           = numeric(0),  # effective obs count for this regime
    postmix_frac      = numeric(0),  # share of total effective obs (sum-to-1 over regimes of a fit)
    est_Sigma_rcond   = numeric(0),  # rcond of the stored est_Sigma (1/condition number)
    est_Sigma_min_eig = numeric(0),  # smallest eigenvalue of est_Sigma
    # --- fit-level, repeated on every regime row of the same fit:
    gamma_entropy_mean = numeric(0), # mean per-timestep posterior entropy (nats), 0..log(M)
    iter               = integer(0)  # EM iterations to convergence
  )
}

# -----------------------------------------------------------------------------
# Regime-separation summaries from the smoothed posterior (model_fit$smoothedprob).
# For the single-series pipeline (N.samples = 1, M > 1) smoothedprob is a
# (T-1) x M matrix; for M == 1 it may be absent. Returns postmix (effective obs
# count per regime, summed over time) and the mean per-timestep posterior entropy
# (a separation-quality measure independent of the raw postmix sum: high entropy =
# the model is chronically unsure which regime it is in). Estimated-regime order.
fit_separation_summary <- function(smoothedprob, M, T_minus_1) {
  if (M == 1 || is.null(smoothedprob)) {
    return(list(postmix = T_minus_1, entropy_mean = 0))
  }
  # Reshape to a (rows x M) matrix of per-timestep posteriors. For N.samples = 1
  # smoothedprob is already (T-1) x M; the 3D [N.samples, time, M] case is folded
  # so the last (regime) dim is preserved as columns.
  if (length(dim(smoothedprob)) == 3) {
    d3 <- dim(smoothedprob)
    P  <- matrix(aperm(smoothedprob, c(2, 1, 3)), ncol = d3[3])
  } else {
    P <- as.matrix(smoothedprob)
  }
  postmix <- colSums(P, na.rm = TRUE)
  # per-row entropy H_t = -sum_j p_tj log p_tj (0 for degenerate rows)
  logP <- ifelse(P > 0, log(P), 0)
  H    <- -rowSums(P * logP, na.rm = TRUE)
  list(postmix = postmix, entropy_mean = mean(H, na.rm = TRUE))
}
empty_fit_row_schema <- function() {
  tibble::tibble(
    timesteps = integer(0), density = numeric(0), nodes = integer(0),
    regimes = integer(0), ts_id = integer(0),
    orig_TPM = list(), est_TPM = list(),
    orig_regime_sequence = list(), est_regime_sequence = list()
  )
}
empty_failure_schema <- function() {
  tibble::tibble(
    timesteps = integer(0), density = numeric(0), nodes = integer(0),
    regimes = integer(0), ts_id = integer(0),
    stage = character(0), error = character(0)
  )
}
bind_or_empty <- function(piece_list, empty_fun) {
  piece_list <- piece_list[!vapply(piece_list, is.null, logical(1))]
  if (length(piece_list) > 0) dplyr::bind_rows(piece_list) else empty_fun()
}

# -----------------------------------------------------------------------------
# Build one payload per DESIGN CELL (T, Density, N, M combination), each
# bundling all n_ts replicates. Dispatching whole cells rather than
# individual replicates as separate future_lapply tasks is the key lever for
# wall-clock time: every dispatched task pays a fixed serialisation/IPC cost
# under multisession (sending data to a worker, returning a result), and with
# n_ts in the tens to hundreds, dispatching per-replicate makes that fixed
# cost a large fraction of total time for what are often only few-second
# fits. Bundling by cell divides the number of dispatches by n_ts, while
# leaving cross-CELL load-balancing untouched -- the real cost heterogeneity
# in this design (N=8/T=4000/M=4 vs. N=4/T=250/M=1) is between cells, not
# between replicates of the same cell, so nothing is lost there.
build_cell_data <- function(Timeseries_data, T, Density, N, M, n_ts) {
  
  cell_grid <- expand.grid(
    k = seq_along(M), j = seq_along(N), i = seq_along(Density), t = seq_along(T),
    KEEP.OUT.ATTRS = FALSE
  )
  n_cells <- nrow(cell_grid)
  cell_data <- vector("list", n_cells)
  
  for (c in seq_len(n_cells)) {
    t <- cell_grid$t[c]; i <- cell_grid$i[c]; j <- cell_grid$j[c]; k <- cell_grid$k[c]
    
    replicates <- vector("list", n_ts)
    for (l in seq_len(n_ts)) {
      current_row <- Timeseries_data %>%
        dplyr::filter(timesteps == T[t], density == Density[i],
                      nodes == N[j], regimes == M[k], ts_id == l)
      
      replicates[[l]] <- list(
        l = l,
        timeseries_data = if (nrow(current_row) > 0) current_row$timeseries_data[[1]] else NULL,
        regime_dynamics = if (nrow(current_row) > 0) current_row$regime_dynamics[[1]] else NULL,
        regime_sequence = if (nrow(current_row) > 0) current_row$regime_sequence[[1]] else NULL,
        transmat        = if (nrow(current_row) > 0) current_row$transmat[[1]]        else NULL
      )
    }
    
    cell_data[[c]] <- list(
      t = t, i = i, j = j, k = k,
      timesteps_val = T[t], density_val = Density[i],
      nodes_val     = N[j], regimes_val = M[k],
      replicates    = replicates
    )
  }
  
  cell_data
}

# -----------------------------------------------------------------------------
# Fit ONE replicate within a cell. Pure function: takes its own minimal data
# slice, returns its own result pieces -- no shared mutable state.
fit_one_replicate <- function(cell, rep, order, MaxIter, verbose) {
  
  loc <- paste0("T=", cell$timesteps_val, ", D=", cell$density_val,
                ", N=", cell$nodes_val, ", M=", cell$regimes_val, ", ts=", rep$l)
  
  make_failure <- function(stage, error_msg) {
    tibble::tibble(
      timesteps = cell$timesteps_val, density = cell$density_val,
      nodes     = cell$nodes_val,     regimes = cell$regimes_val, ts_id = rep$l,
      stage     = stage,
      error     = if (is.null(error_msg) || length(error_msg) == 0) NA_character_
      else paste(trimws(error_msg), collapse = " ")
    )
  }
  
  if (is.null(rep$timeseries_data)) {
    message("No timeseries data found for condition: ", loc)
    return(list(regime_rows = NULL, fit_row = NULL,
                failure = make_failure("no_data", "no timeseries data found for condition")))
  }
  
  if (verbose) {
    print("", quote = FALSE)
    print(paste("Timeseries:", loc), quote = FALSE)
  }
  
  current_ts      <- rep$timeseries_data
  current_regimes <- rep$regime_dynamics
  current_ts_norm <- current_ts  # Normalisation hook (currently identity)
  
  timesteps_t <- cell$timesteps_val
  d           <- ncol(current_ts_norm)
  current_ts_norm_array <- array(data = current_ts_norm, dim = c(timesteps_t, 1, d))
  
  # Fit MSAR model (EM + LASSO + retry; unchanged)
  result    <- init_and_fit_msar_lasso(
    data    = current_ts_norm_array,
    M       = cell$regimes_val,
    order   = order,
    MaxIter = MaxIter,
    retry   = 5,
    verbose = verbose
  )
  model_fit <- result[["fit"]]
  error     <- result[["error"]]
  
  if (is.null(model_fit)) {
    message("ignore fit! ", "\n\tTimeseries_data: ", loc, "\n\terror: ", error)
    return(list(regime_rows = NULL, fit_row = NULL,
                failure = make_failure("fit_null", error)))
  }
  
  # --- Raw estimates -----------------------------------------------------
  est_Betas  <- model_fit[["theta"]][["A"]]      # list per regime, each $A1
  est_sigmas <- model_fit[["theta"]][["sigma"]]  # list per regime
  est_TPM    <- model_fit[["theta"]][["transmat"]]
  
  est_Betas_mat <- setNames(
    lapply(seq_len(cell$regimes_val), function(m) est_Betas[[m]][["A1"]]),
    paste0("Regime", seq_len(cell$regimes_val))
  )
  orig_Betas <- setNames(
    lapply(seq_len(cell$regimes_val), function(m) current_regimes[[paste0("Regime", m)]][["Beta"]]),
    paste0("Regime", seq_len(cell$regimes_val))
  )
  
  # --- Regime matching (Hungarian on raw vectorised Beta) ----------------
  assigned_regimes <- match_regimes(orig_Betas, est_Betas_mat)
  if (is.null(assigned_regimes)) {
    return(list(regime_rows = NULL, fit_row = NULL,
                failure = make_failure("zero_var_beta",
                                       "estimated Beta has zero variance in at least one regime")))
  }
  
  # --- Sigma/Kappa degeneracy diagnostics (uncommitted; see doc) ----------
  # Regime separation from the smoothed posterior (estimated-regime order), and
  # per-regime conditioning of the STORED est_Sigma. postmix[est_idx] is the
  # effective observation count that produced est_sigmas[[est_idx]]; the leading
  # hypothesis is that low postmix drives the near-singular / magnitude-degenerate
  # Kappa. rcond/min-eig are computed on exactly the matrix that analysis inverts.
  sep <- fit_separation_summary(model_fit[["smoothedprob"]],
                                cell$regimes_val, cell$timesteps_val - 1L)
  postmix_vec  <- sep$postmix
  postmix_frac <- postmix_vec / sum(postmix_vec)
  entropy_mean <- sep$entropy_mean
  em_iters     <- if (!is.null(model_fit[["Iter"]])) as.integer(model_fit[["Iter"]]) else NA_integer_
  sig_rcond <- function(S) { S <- as.matrix(S); tryCatch(rcond(S), error = function(e) NA_real_) }
  sig_mineig <- function(S) {
    S <- as.matrix(S)
    tryCatch(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values),
             error = function(e) NA_real_)
  }

  # --- One raw row per regime (true-regime order 1..M) --------------------
  regime_rows <- vector("list", nrow(assigned_regimes))
  for (m in seq_len(nrow(assigned_regimes))) {
    orig_idx <- assigned_regimes[[m, 1]]   # true regime label (== m)
    est_idx  <- assigned_regimes[[m, 2]]   # matched estimated regime label
    orig_reg <- current_regimes[[paste0("Regime", orig_idx)]]
    est_S    <- est_sigmas[[est_idx]]

    regime_rows[[m]] <- tibble::tibble(
      timesteps = cell$timesteps_val, density = cell$density_val,
      nodes     = cell$nodes_val,     regimes = cell$regimes_val,
      ts_id     = rep$l,              regime_id = orig_idx,
      orig_Beta    = list(orig_reg[["Beta"]]),
      est_Beta     = list(est_Betas_mat[[est_idx]]),
      orig_Kappa   = list(orig_reg[["kappa"]]),
      orig_Sigma   = list(orig_reg[["sigma"]]),
      est_Sigma    = list(est_S),
      orig_Beta_ac = list(orig_reg[["Beta_ac"]]),
      postmix            = postmix_vec[est_idx],
      postmix_frac       = postmix_frac[est_idx],
      est_Sigma_rcond    = sig_rcond(est_S),
      est_Sigma_min_eig  = sig_mineig(est_S),
      gamma_entropy_mean = entropy_mean,
      iter               = em_iters
    )
  }
  regime_rows <- dplyr::bind_rows(regime_rows)
  
  # --- Fit-level row: TPM + regime sequences, relabelled to true order ----
  perm <- assigned_regimes[, 2]
  est_TPM_relabelled <- est_TPM[perm, perm, drop = FALSE]
  dimnames(est_TPM_relabelled) <- list(paste0("Regime", seq_len(cell$regimes_val)),
                                       paste0("Regime", seq_len(cell$regimes_val)))
  
  # True sequence alignment: regime_sequence has length (totTime - 1); the
  # trimmed series keeps the last T rows, and the order-1 fit drops the first
  # of those as the initial AR lag -> smoothed probs have T - 1 rows. Take the
  # last T true entries, drop the first to match.
  true_seq_full <- rep$regime_sequence
  true_seq      <- utils::tail(true_seq_full, cell$timesteps_val)[-1]
  
  # M == 1 has no latent switching: the decoded sequence is trivially
  # all-regime-1, and smoothedprob may be absent -- handle without touching
  # the EM object. Sequence recovery (RQ4) is analysed for regimes >= 2 only.
  if (cell$regimes_val > 1) {
    seq_mapped <- reconstruct_regime_sequence(
      model_fit[["smoothedprob"]], assigned_regimes
    )$mapped
  } else {
    seq_mapped <- rep(1L, length(true_seq))
  }
  
  fit_row     <- NULL
  seq_failure <- NULL  # only set below if the sequence row has to be skipped
  if (length(true_seq) == length(seq_mapped)) {
    fit_row <- tibble::tibble(
      timesteps = cell$timesteps_val, density = cell$density_val,
      nodes     = cell$nodes_val,     regimes = cell$regimes_val, ts_id = rep$l,
      orig_TPM = list(rep$transmat),
      est_TPM  = list(est_TPM_relabelled),
      orig_regime_sequence = list(true_seq),
      est_regime_sequence  = list(seq_mapped)
    )
  } else {
    # NOTE: this is NOT a failure for Beta_corr/Kappa_corr/Beta_ac_corr_pearson --
    # regime_rows above is already built and is returned regardless of what
    # happens here. Only the RQ4 sequence row (fit_row) is unavailable for this
    # replicate. Logged via the same make_failure() mechanism as fit_null/
    # zero_var_beta/no_data (previously this branch only printed a message()
    # and returned failure = NULL, so it was invisible to summarize_failures()
    # and to every downstream feasibility/sensitivity analysis -- this is what
    # produced the ~21% "extra", untracked missingness in RQ4 relative to the
    # corr-outcome failure rate).
    #
    # Diagnostic fields (parseable via e.g. stringr::str_match(error,
    # "smoothedprob_rows=(\\d+)")) let summarize_failures()'s existing by_cell
    # breakdown show WHERE (which T/N/M/Density) the mismatch concentrates, and
    # the length comparison shows WHICH side is off: true_seq vs. true_seq_full
    # (a generation/trimming-side issue) or smoothedprob_rows vs. the T-1
    # expectation documented in reconstruct_regime_sequence() (an estimation-
    # side issue, e.g. EM convergence affecting how many rows smoothedprob has).
    n_smoothed <- if (cell$regimes_val > 1) nrow(model_fit[["smoothedprob"]]) else NA_integer_
    diag_msg <- sprintf(
      paste("seq length mismatch: true_seq_full=%d true_seq=%d",
            "smoothedprob_rows=%s est_seq=%d expected_T_minus_1=%d order=%d"),
      length(true_seq_full), length(true_seq),
      if (is.na(n_smoothed)) "NA" else as.character(n_smoothed),
      length(seq_mapped), cell$timesteps_val - 1L, order
    )
    message("Skipping fit-level sequence row for ", loc, ": ", diag_msg)
    seq_failure <- make_failure("seq_length_mismatch", diag_msg)
  }
  
  list(regime_rows = regime_rows, fit_row = fit_row, failure = seq_failure)
}

# -----------------------------------------------------------------------------
# Fit an entire CELL (all n_ts replicates). This -- not fit_one_replicate --
# is the unit dispatched to future_lapply. Packages are attached once here
# (not per replicate): cheap after the first call on a given persistent
# worker, so doing it at the cell level rather than the replicate level
# saves nothing functionally but keeps the call site singular and obvious.
fit_one_cell <- function(cell, order, MaxIter, verbose, progress_fun = NULL,
                         lasso_control = NULL) {

  ensure_worker_packages()

  # Apply the LASSO / penalization options INSIDE the worker process. options()
  # are process-level state and are NOT exported to future workers, so setting
  # them in the calling (main) process would have no effect on the parallel
  # fits -- they must be (re)applied here, in the worker, before any fit_msar
  # call. lasso_control is a named list of simmsar_lasso_* options built by the
  # caller (see estimate_MSAR() / MSAR_ts_analysis_NHMSAR.R).
  if (!is.null(lasso_control) && length(lasso_control) > 0) {
    do.call(options, lasso_control)
  }

  regime_list  <- vector("list", length(cell$replicates))
  fit_list     <- vector("list", length(cell$replicates))
  failure_list <- vector("list", length(cell$replicates))
  
  for (r in seq_along(cell$replicates)) {
    out <- fit_one_replicate(cell, cell$replicates[[r]], order, MaxIter, verbose)
    regime_list[[r]]  <- out$regime_rows
    fit_list[[r]]     <- out$fit_row
    failure_list[[r]] <- out$failure
    if (!is.null(progress_fun)) progress_fun()
  }
  
  list(
    regime_rows = bind_or_empty(regime_list,  empty_regime_rows_schema),
    fit_row     = bind_or_empty(fit_list,     empty_fit_row_schema),
    failure     = bind_or_empty(failure_list, empty_failure_schema)
  )
}


estimate_MSAR <- function(Density, N, M, T, n_ts, order, MaxIter, verbose, min_edg_val,
                          Timeseries_data, workers = NULL, lasso_control = NULL) {
  
  if (is.null(workers)) {
    workers <- max(1, parallel::detectCores(logical = FALSE) - 1)
  }
  
  print(sprintf("Estimating MSAR models from timeseries (%d parallel worker%s)",
                workers, if (workers == 1) "" else "s"), quote = FALSE)
  
  # future::multisession works identically on Windows/Mac/Linux (unlike
  # future::multicore, which silently falls back to sequential on Windows).
  # workers = 1 collapses to fully sequential execution (also useful as a
  # debugging / sanity-check fallback).
  old_plan <- future::plan()
  if (workers <= 1) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Build the per-cell payload BEFORE going parallel, then drop the (possibly
  # multi-GB) full Timeseries_data tibble from the calling process: each
  # worker will only receive the slice of cell_data it actually runs, via
  # future_lapply's own chunking of its X argument -- NOT a closure-captured
  # copy of the whole object. Dispatching whole CELLS (all n_ts replicates
  # bundled) rather than individual replicates is what keeps per-dispatch
  # IPC/serialisation overhead from dominating wall-clock time -- see the
  # comment on build_cell_data() above.
  cat("Slicing time series data into per-cell payloads...\n")
  cell_data <- build_cell_data(Timeseries_data, T, Density, N, M, n_ts)
  n_cells   <- length(cell_data)
  n_total   <- n_cells * n_ts
  rm(Timeseries_data); invisible(gc(verbose = FALSE))
  
  # Live cross-process progress bar (the old `progress` package's progress_bar
  # only updates a LOCAL copy inside whichever worker calls tick() -- it never
  # reaches the main console under multisession. progressr is the standard
  # future-aware replacement: each worker's p() call is relayed back here).
  # Progress is still reported per REPLICATE (n_total steps), even though the
  # dispatch unit is the cell -- fit_one_cell() calls progress_fun() after
  # each of its n_ts replicates, so the bar doesn't go quiet for the full
  # duration of a large cell.
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers(progressr::handler_progress(
      format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining"
    ))
    results <- progressr::with_progress({
      p <- progressr::progressor(steps = n_total)
      # Pass fit_one_cell BARE (with order/MaxIter/verbose/progress_fun as named
      # ... args), exactly like the no-progressr branch below. Do NOT wrap it in
      # an anonymous `function(cell) fit_one_cell(...)` defined here: that closure
      # would close over THIS (estimate_MSAR) evaluation frame, which still holds
      # the multi-GB `cell_data`. future serialises an exported function together
      # with its environment, so the wrapper would drag a full copy of cell_data
      # into every future as the global 'FUN' -- tripping
      # future.globals.maxSize (the "FUN is 900 MiB of class function" error) and
      # duplicating the whole dataset to each worker. fit_one_cell lives in the
      # global env (sourced), so passing it bare exports only the small function,
      # while cell_data still ships correctly -- and only sliced -- as X.
      future.apply::future_lapply(
        cell_data,
        fit_one_cell,
        order = order, MaxIter = MaxIter, verbose = verbose, progress_fun = p,
        lasso_control = lasso_control,
        future.seed = TRUE,
        future.scheduling = Inf  # one CELL per dispatch: cost is highly
        # heterogeneous across cells (e.g. N=8/
        # T=4000/M=4 vs. N=4/T=250/M=1), so
        # fine-grained scheduling balances workers
        # much better than pre-chunking. Bundling
        # all n_ts replicates per cell already keeps
        # the number of dispatches small (n_cells,
        # not n_cells*n_ts), so this no longer
        # incurs per-replicate dispatch overhead.
      )
    })
  } else {
    message("Package 'progressr' not installed -- running without a live progress bar.\n",
            "Install it (install.packages(\"progressr\")) for progress reporting across workers.")
    results <- future.apply::future_lapply(
      cell_data,
      fit_one_cell, order = order, MaxIter = MaxIter, verbose = verbose,
      lasso_control = lasso_control,
      future.seed = TRUE, future.scheduling = Inf
    )
  }
  
  # ---------------------------------------------------------------------------
  # Combine the per-cell results returned by all workers (each is already
  # internally bound across that cell's n_ts replicates, see fit_one_cell()).
  # ---------------------------------------------------------------------------
  regime_list  <- lapply(results, `[[`, "regime_rows")
  fit_list     <- lapply(results, `[[`, "fit_row")
  failure_list <- lapply(results, `[[`, "failure")
  
  msar_results <- bind_or_empty(regime_list,  empty_regime_rows_schema)
  fit_results  <- bind_or_empty(fit_list,     empty_fit_row_schema)
  failure_log  <- bind_or_empty(failure_list, empty_failure_schema)
  
  class(msar_results) <- c("msar_results", class(msar_results))
  attr(msar_results, "fit_results") <- fit_results
  attr(msar_results, "failure_log") <- failure_log
  # Persist the exact LASSO/penalization config that produced this object, so a
  # diagnostic run's engine settings never have to be reconstructed after the
  # fact (see docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md). NULL means library
  # defaults (engine="bic") were in force.
  attr(msar_results, "lasso_control") <- lasso_control

  msar_results
}


#' Summarise Discarded-Fit Reasons from estimate_MSAR()
#'
#' Tabulates the diagnostic failure log attached as
#' \code{attr(result, "failure_log")} by coarse error type, estimation stage,
#' and design cell. Use this to diagnose the missingness mechanism: which kind
#' of numerical failure dominates which condition.
#'
#' @param msar_results An msar_results object from \code{estimate_MSAR()}.
#' @param by Character vector of design columns to break failures down by.
#'   Default groups by all four design factors.
#'
#' @return Invisibly, a list with \code{by_type}, \code{by_cell}, \code{raw}.
#'   Also prints \code{by_type} for a quick overview.
#'
#' @details
#' \code{stage} distinguishes \emph{where} the fit was dropped:
#' \itemize{
#'   \item \code{fit_null}: EM aborted on error/warning across all retries.
#'   \item \code{zero_var_beta}: an estimated Beta had zero variance, so regime
#'     matching was undefined and the fit was discarded.
#'   \item \code{no_data}: no simulated series matched the condition (lookup miss).
#'   \item \code{seq_length_mismatch}: the EM fit and regime matching succeeded
#'     (so this replicate's regime_rows -- and hence Beta_corr/Kappa_corr/
#'     Beta_ac_corr_pearson -- are unaffected and already included), but the
#'     true and reconstructed regime sequences had different lengths, so the
#'     fit-level sequence row needed for RQ4 could not be built. RQ4-specific
#'     missingness only; does not reduce the sample for RQ1--RQ3.
#' }
#' The near-singular \code{est_Sigma} case is NOT logged here any more -- it is a
#' per-regime validity flag handled in analysis (\code{sigma_validity_log}).
#'
#' @seealso \code{\link{estimate_MSAR}}, \code{summarize_sigma_validity}
#' @export
summarize_failures <- function(msar_results,
                               by = c("timesteps", "density", "nodes", "regimes")) {
  
  fl <- attr(msar_results, "failure_log")
  if (is.null(fl) || nrow(fl) == 0) {
    message("No failures logged.")
    return(invisible(NULL))
  }
  
  classify <- function(x) {
    x <- tolower(ifelse(is.na(x), "", x))
    dplyr::case_when(
      grepl("seq length mismatch", x) ~ "RQ4 sequence-length mismatch (not a fit failure)",
      grepl("singular|rcond|cxx|det\\(|positive.?definite|chol", x) ~ "singular/ill-conditioned covariance",
      grepl("zero variance|sd", x) ~ "degenerate (zero-variance) estimate",
      grepl("nan|non-finite|infinite|\\binf\\b|missing value|na/nan", x) ~ "non-finite / NaN in likelihood",
      grepl("lars|glmnet|lasso|lambda|penal", x) ~ "LASSO / penalised-path issue",
      grepl("no timeseries data", x) ~ "no data (lookup miss)",
      x == "" ~ "unclassified (empty)",
      TRUE ~ "other"
    )
  }
  
  fl <- fl %>% dplyr::mutate(error_type = classify(error))
  
  by_type <- fl %>% dplyr::count(stage, error_type, name = "n", sort = TRUE)
  by_cell <- fl %>%
    dplyr::count(dplyr::across(dplyr::all_of(by)), error_type, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  
  cat("=== DISCARDED FITS: reason breakdown ===\n")
  cat(sprintf("Total logged failures: %d\n\n", nrow(fl)))
  print(by_type, n = Inf)
  cat("\nUse $by_cell for the per-condition breakdown, $raw for the full log.\n")
  
  invisible(list(by_type = by_type, by_cell = by_cell, raw = fl))
}


#' Print Method for msar_results
#'
#' @param x An msar_results object
#' @param ... Additional arguments (unused)
#'
#' @export
print.msar_results <- function(x, ...) {
  cat("MSAR Results (raw estimates)\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Total observations (regime rows): %d\n", nrow(x)))
  cat("Conditions tested:\n")
  cat(sprintf("  Timesteps: %s\n", paste(unique(x$timesteps), collapse = ", ")))
  cat(sprintf("  Density: %s\n", paste(unique(x$density), collapse = ", ")))
  cat(sprintf("  Nodes: %s\n", paste(unique(x$nodes), collapse = ", ")))
  cat(sprintf("  Regimes: %s\n", paste(unique(x$regimes), collapse = ", ")))
  cat("════════════════════════════════════════════════════════════════\n")
  cat("Raw estimates only; derived metrics are computed in the analysis pipeline\n")
  cat("(compute_recovery_metrics()). Fit-level TPM/sequences in attr 'fit_results'.\n")
  NextMethod()
}