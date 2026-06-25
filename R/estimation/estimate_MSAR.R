#' Estimate MSAR Models from Time Series Data (raw estimates only)
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
#' @param verbose Logical. If TRUE, prints detailed progress messages.
#' @param min_edg_val Numeric. Retained for call-site compatibility but
#'   \strong{no longer used here}: edge thresholding is a scoring decision and
#'   now lives in the analysis pipeline (\code{compute_recovery_metrics()}),
#'   which applies the threshold to the stored raw \code{est_Beta} / inverted
#'   \code{est_Kappa}. Pass the intended threshold there instead.
#' @param Timeseries_data Tibble. Output from \code{generate_timeseries} containing
#'   simulated data and true dynamics in tibble format with list-columns.
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


estimate_MSAR <- function(Density, N, M, T, n_ts, order, MaxIter, verbose, min_edg_val,
                          Timeseries_data) {

  # Setup progressbar
  print("Estimating MSAR models from timeseries", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(T) * length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )

  # PERFORMANCE: Pre-allocate the per-regime tibble to its maximum possible size
  # (avoids costly row-binding in the loop). Raw quantities only.
  max_rows <- length(T) * length(Density) * length(N) * length(M) * n_ts * max(M)

  msar_results <- tibble::tibble(
    timesteps = integer(max_rows),
    density   = numeric(max_rows),
    nodes     = integer(max_rows),
    regimes   = integer(max_rows),
    ts_id     = integer(max_rows),
    regime_id = integer(max_rows),

    # Raw network matrices (list-columns). est_Beta / est_Sigma stored as-is.
    orig_Beta    = vector("list", max_rows),
    est_Beta     = vector("list", max_rows),
    orig_Kappa   = vector("list", max_rows),
    orig_Sigma   = vector("list", max_rows),
    est_Sigma    = vector("list", max_rows),
    orig_Beta_ac = vector("list", max_rows)
  )
  row_idx <- 0

  # Pre-allocate the fit-level table (one row per fit; TPM + sequences are
  # properties of the whole fit, not of an individual regime).
  max_fit_rows <- length(T) * length(Density) * length(N) * length(M) * n_ts

  fit_results <- tibble::tibble(
    timesteps = integer(max_fit_rows),
    density   = numeric(max_fit_rows),
    nodes     = integer(max_fit_rows),
    regimes   = integer(max_fit_rows),
    ts_id     = integer(max_fit_rows),
    orig_TPM  = vector("list", max_fit_rows),
    est_TPM   = vector("list", max_fit_rows),
    orig_regime_sequence = vector("list", max_fit_rows),
    est_regime_sequence  = vector("list", max_fit_rows)
  )
  fit_idx <- 0

  # Diagnostic failure log: records WHY and WHERE each fit was discarded, so the
  # missingness mechanism is inspectable. Attached as attr "failure_log";
  # summarise with summarize_failures(). NOTE: the old whole-fit drop on a
  # near-singular est_Sigma is gone -- that is now a per-regime NA flag applied
  # in analysis (sigma_validity_log), a distinct mechanism from this log.
  failure_records <- list()
  log_failure <- function(stage, error_msg) {
    failure_records[[length(failure_records) + 1L]] <<- tibble::tibble(
      timesteps = T[t],
      density   = Density[i],
      nodes     = N[j],
      regimes   = M[k],
      ts_id     = l,
      stage     = stage,
      error     = if (is.null(error_msg) || length(error_msg) == 0)
                    NA_character_
                  else paste(trimws(error_msg), collapse = " ")
    )
  }

  for (t in seq_along(T)) {
    for (i in seq_along(Density)) {
      for (j in seq_along(N)) {
        for (k in seq_along(M)) {
          for (l in 1:n_ts) {

            loc <- paste0("t=", t, ", i=", i, ", j=", j, ", k=", k, ", l=", l)

            if (verbose) {
              print("", quote = FALSE)
              print(paste("Timeseries:", loc), quote = FALSE)
            }

            # Get complete data for current timeseries from tibble
            current_row <- Timeseries_data %>%
              filter(timesteps == T[t], density == Density[i],
                     nodes == N[j], regimes == M[k], ts_id == l)

            if (nrow(current_row) == 0) {
              message("No timeseries data found for condition: ",
                      "T=", T[t], " D=", Density[i], " N=", N[j], " M=", M[k], " ts=", l)
              log_failure("no_data", "no timeseries data found for condition")
              pb$tick()
              next
            }

            current_ts      <- current_row$timeseries_data[[1]]
            current_regimes <- current_row$regime_dynamics[[1]]

            # Normalisation hook (currently identity; see generate_timeseries.R)
            current_ts_norm <- current_ts

            # Make 3D array (time, samples = 1, nodes) expected by the fitter
            timesteps <- T[t]
            d         <- ncol(current_ts_norm)
            current_ts_norm_array <- array(data = current_ts_norm,
                                           dim = c(timesteps, 1, d))

            # Fit MSAR model (EM + LASSO + retry; unchanged)
            result    <- init_and_fit_msar_lasso(
              data    = current_ts_norm_array,
              M       = M[k],
              order   = order,
              MaxIter = MaxIter,
              retry   = 5,
              verbose = verbose
            )
            model_fit <- result[["fit"]]
            error     <- result[["error"]]

            if (is.null(model_fit)) {
              message("ignore fit! ", "\n\tTimeseries_data: ", loc, "\n\terror: ", error)
              log_failure("fit_null", error)
              pb$tick()
              next
            }

            # --- Raw estimates ---------------------------------------------------
            est_Betas  <- model_fit[["theta"]][["A"]]      # list per regime, each $A1
            est_sigmas <- model_fit[["theta"]][["sigma"]]  # list per regime
            est_TPM    <- model_fit[["theta"]][["transmat"]]

            # Named lists (Regime1..M) of raw matrices, in regime order.
            est_Betas_mat <- setNames(
              lapply(seq_len(M[k]), function(m) est_Betas[[m]][["A1"]]),
              paste0("Regime", seq_len(M[k]))
            )
            orig_Betas <- setNames(
              lapply(seq_len(M[k]), function(m) current_regimes[[paste0("Regime", m)]][["Beta"]]),
              paste0("Regime", seq_len(M[k]))
            )

            # --- Regime matching (Hungarian on raw vectorised Beta) --------------
            assigned_regimes <- match_regimes(orig_Betas, est_Betas_mat)
            if (is.null(assigned_regimes)) {
              # Degenerate (zero-variance) estimated Beta -> cannot match.
              log_failure("zero_var_beta", "estimated Beta has zero variance in at least one regime")
              pb$tick()
              next
            }

            # --- Store one raw row per regime (true-regime order 1..M) -----------
            for (m in seq_len(nrow(assigned_regimes))) {
              orig_idx <- assigned_regimes[[m, 1]]   # true regime label (== m)
              est_idx  <- assigned_regimes[[m, 2]]   # matched estimated regime label
              row_idx  <- row_idx + 1

              orig_reg <- current_regimes[[paste0("Regime", orig_idx)]]

              msar_results$timesteps[row_idx] <- T[t]
              msar_results$density[row_idx]   <- Density[i]
              msar_results$nodes[row_idx]     <- N[j]
              msar_results$regimes[row_idx]   <- M[k]
              msar_results$ts_id[row_idx]     <- l
              msar_results$regime_id[row_idx] <- orig_idx

              msar_results$orig_Beta[[row_idx]]    <- orig_reg[["Beta"]]
              msar_results$est_Beta[[row_idx]]     <- est_Betas_mat[[est_idx]]
              msar_results$orig_Kappa[[row_idx]]   <- orig_reg[["kappa"]]
              msar_results$orig_Sigma[[row_idx]]   <- orig_reg[["sigma"]]
              msar_results$est_Sigma[[row_idx]]    <- est_sigmas[[est_idx]]
              msar_results$orig_Beta_ac[[row_idx]] <- orig_reg[["Beta_ac"]]
            }

            # --- Fit-level row: TPM + regime sequences, relabelled to true order -
            # Permutation: for true regime t1, the matching estimated label.
            perm <- assigned_regimes[, 2]
            est_TPM_relabelled <- est_TPM[perm, perm, drop = FALSE]
            dimnames(est_TPM_relabelled) <- list(paste0("Regime", seq_len(M[k])),
                                                 paste0("Regime", seq_len(M[k])))

            # True sequence alignment: regime_sequence has length (totTime - 1);
            # the trimmed series keeps the last T rows, and the order-1 fit drops
            # the first of those as the initial AR lag -> smoothed probs have
            # T - 1 rows. Take the last T true entries, drop the first to match.
            true_seq_full <- current_row$regime_sequence[[1]]
            true_seq      <- utils::tail(true_seq_full, T[t])[-1]

            # M == 1 has no latent switching: the decoded sequence is trivially
            # all-regime-1, and smoothedprob may be absent -- handle without
            # touching the EM object. Sequence recovery (RQ4) is analysed for
            # regimes >= 2 only anyway.
            if (M[k] > 1) {
              seq_mapped <- reconstruct_regime_sequence(
                model_fit[["smoothedprob"]], assigned_regimes
              )$mapped
            } else {
              seq_mapped <- rep(1L, length(true_seq))
            }

            # Guard against a length mismatch (e.g. unexpected smoothedprob
            # dimensions): store the fit-level row only when the sequences align,
            # so downstream accuracy/Cohen's kappa is well-defined.
            if (length(true_seq) == length(seq_mapped)) {
              fit_idx <- fit_idx + 1
              fit_results$timesteps[fit_idx] <- T[t]
              fit_results$density[fit_idx]   <- Density[i]
              fit_results$nodes[fit_idx]     <- N[j]
              fit_results$regimes[fit_idx]   <- M[k]
              fit_results$ts_id[fit_idx]     <- l
              fit_results$orig_TPM[[fit_idx]] <- current_row$transmat[[1]]
              fit_results$est_TPM[[fit_idx]]  <- est_TPM_relabelled
              fit_results$orig_regime_sequence[[fit_idx]] <- true_seq
              fit_results$est_regime_sequence[[fit_idx]]  <- seq_mapped
            } else {
              message("Skipping fit-level sequence row for ", loc,
                      ": length mismatch (true = ", length(true_seq),
                      ", est = ", length(seq_mapped), ")")
            }

            pb$tick()
          }
        }
      }
    }
  }

  # Trim to actual size. seq_len() (not 1:n) so n == 0 selects nothing.
  msar_results <- msar_results[seq_len(row_idx), ]
  fit_results  <- fit_results[seq_len(fit_idx), ]

  class(msar_results) <- c("msar_results", class(msar_results))
  attr(msar_results, "fit_results") <- fit_results

  failure_log <- if (length(failure_records) > 0) {
    dplyr::bind_rows(failure_records)
  } else {
    tibble::tibble(
      timesteps = integer(0), density = numeric(0), nodes = integer(0),
      regimes = integer(0), ts_id = integer(0),
      stage = character(0), error = character(0)
    )
  }
  attr(msar_results, "failure_log") <- failure_log

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
