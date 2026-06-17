#' Estimate MSAR Models from Time Series Data
#'
#' Main estimation pipeline that fits Markov-Switching Autoregressive (MSAR) models
#' to time series data and compares estimated network dynamics to true dynamics.
#' Processes multiple conditions in parallel and computes recovery statistics.
#'
#' @param Density Numeric vector. Edge densities used in data generation.
#' @param N Integer vector. Numbers of nodes in networks.
#' @param M Integer vector. Numbers of regimes.
#' @param T Integer vector. Time series lengths.
#' @param n_ts Integer. Number of time series per condition.
#' @param order Integer. AR order for the VAR model (typically 1).
#' @param MaxIter Integer. Maximum EM algorithm iterations.
#' @param verbose Logical. If TRUE, prints detailed progress messages.
#' @param min_edg_val Numeric. Minimum edge threshold - values below are set to zero.
#'   Only used when use_bootstrap = FALSE. Never applied to the diagonal of Kappa.
#' @param Timeseries_data Tibble. Output from \code{generate_timeseries} containing
#'   simulated data and true dynamics in tibble format with list-columns.
#' @param use_bootstrap Logical. If TRUE, use bootstrap stability selection instead
#'   of simple thresholding. Default: FALSE. Note: This is computationally expensive.
#' @param n_bootstrap Integer. Number of bootstrap samples when use_bootstrap = TRUE.
#'   Default: 50. Higher values give more stable selection but take longer.
#' @param bootstrap_threshold Numeric. Selection threshold (0-1) for bootstrap.
#'   Edges selected in > threshold proportion of bootstraps are kept. Default: 0.6.
#'
#' @return An msar_results tibble (one row per regime per successfully fitted
#'   time series) with the following columns:
#'   \describe{
#'     \item{timesteps, density, nodes, regimes, ts_id, regime_id}{Condition identifiers}
#'     \item{orig_Beta, est_Beta}{True / estimated temporal network (list-columns)}
#'     \item{orig_Kappa, est_Kappa}{True / estimated contemporaneous network (list-columns)}
#'     \item{orig_Beta_ac, est_Beta_ac}{True / estimated Beta average controllability (list-columns)}
#'     \item{Beta_corr}{Correlation between true and estimated Beta (full matrix)}
#'     \item{Beta_sen, Beta_spec}{Sensitivity/specificity for Beta edge detection (full matrix)}
#'     \item{MAE_Beta}{Mean absolute error for Beta}
#'     \item{Beta_ac_corr}{Correlation of average controllability for Beta}
#'     \item{Kappa_corr}{Correlation between true and estimated Kappa (off-diagonal upper triangle)}
#'     \item{Kappa_sen, Kappa_spec}{Sensitivity/specificity for Kappa edge detection (off-diagonal)}
#'     \item{MAE_Kappa}{Mean absolute error for Kappa (off-diagonal)}
#'   }
#'
#'   Regime-sequence recovery statistics (accuracy, Cohen's kappa; M > 1 only)
#'   are attached as a separate tibble via \code{attr(result, "sequence_results")},
#'   so that the primary return value remains a plain msar_results tibble for
#'   backward compatibility with the rest of the pipeline (\code{stat_analysis.R}
#'   etc. expect \code{inherits(MSAR_dynamics_list, "msar_results")} to hold and
#'   treat it directly as a data frame).
#'
#' @details
#' The function processes each time series through the following pipeline:
#'
#' 1. **Normalization**: Applies nonparanormal transformation via \code{huge.npn}
#'
#' 2. **Model Estimation**: Fits MSAR model using \code{init_and_fit_msar_lasso}
#'    with EM algorithm and LASSO regularization
#'
#' 3. **Network Recovery**:
#'    \itemize{
#'      \item Extracts Beta (temporal) and Sigma (covariance) from fitted model
#'      \item Computes Kappa = Sigma^{-1} (precision matrix), keeping its full diagonal
#'      \item Applies min_edg_val threshold to estimated Beta (full matrix) and
#'        to estimated Kappa (off-diagonal only)
#'    }
#'
#' 4. **Regime Assignment**: For multi-regime models, matches estimated to true
#'    regimes based on maximum Beta correlations using \code{assign_regimes}.
#'    This single Beta-based mapping is applied consistently to Beta recovery,
#'    Kappa recovery, Beta-AC recovery, and the regime-sequence recovery.
#'
#' 5. **Comparison Metrics**: Computes for each regime pair:
#'    \itemize{
#'      \item Correlation: Pearson correlation between vectorized networks
#'      \item Sensitivity/Specificity: Edge detection accuracy
#'      \item MAE: Mean absolute error of edge weights
#'      \item Average controllability correlation (Beta only)
#'    }
#'
#' 6. **Regime-Sequence Recovery** (M > 1 only): The hard (argmax) estimated
#'    regime sequence is derived from the smoothed probabilities, relabeled
#'    using the Beta-based regime mapping, and compared to the true regime
#'    sequence via accuracy and Cohen's kappa.
#'
#' Failed estimations (NULL fits or zero-variance networks) are skipped with warning messages.
#'
#' @note
#' Dependencies are loaded centrally via R/dependencies.R
#' Required packages: huge, NHMSAR, dplyr, tibble, progress
#'
#' The function displays progress bars for model estimation.
#'
#' @seealso
#' \code{\link{generate_timeseries}} for generating input data
#' \code{\link{init_and_fit_msar_lasso}} for model fitting
#' \code{\link{assign_regimes}} for regime matching
#' \code{\link{senspec}} for sensitivity/specificity
#' \code{\link{calculate_MAE}} for mean absolute error
#' \code{\link{get_stats}} for computing summary statistics
#'
#' @examples
#' \dontrun{
#' # Generate data
#' ts_data <- generate_timeseries(
#'   Density = 0.3, N = 4, M = 2, T = 1000,
#'   n_ts = 10, warmup = 50, totTime = 1050,
#'   mean_rep = 10, sd_rep = 2,
#'   min_edg_val = 0.05, max_edg_val = 1,
#'   remain_lower = 0.33, remain_upper = 0.66
#' )
#'
#' # Estimate models
#' results <- estimate_MSAR(
#'   Density = 0.3, N = 4, M = 2, T = 1000,
#'   n_ts = 10, order = 1, MaxIter = 200,
#'   verbose = FALSE, min_edg_val = 0.05,
#'   Timeseries_data = ts_data
#' )
#'
#' # Access results (a tibble!)
#' results %>% filter(timesteps == 1000, density == 0.3)
#'
#' # Get summary statistics
#' get_stats(results)
#'
#' # Regime-sequence recovery (M > 1 only)
#' attr(results, "sequence_results")
#' }
#'
#' @export
# Dependencies are loaded centrally via R/dependencies.R
# Required packages: huge, NHMSAR, dplyr, tibble, progress

# Load functions from NHMSAR
source("R/estimation/fit_msar.R")
source("R/estimation/init_theta_msar.R")
source("R/estimation/init_and_fit_msar_lasso.R")

# Load functions to assign and compare regime dynamics
source("R/utils/assign_regimes.R")
source("R/utils/senspec.R")
source("R/utils/summarize_cor.R")
source("R/utils/calculate_MAE.R")
source("R/utils/undirected_metrics.R")
source("R/utils/average_controllability.R")
source("R/utils/regime_sequence_recovery.R")

# Load bootstrap stability selection
source("R/estimation/bootstrap_stability.R")



estimate_MSAR <- function(Density, N, M, T, n_ts, order, MaxIter, verbose, min_edg_val, Timeseries_data,
                          use_bootstrap = FALSE, n_bootstrap = 50, bootstrap_threshold = 0.6) {

# Setup progressbar
print("Estimating MSAR models from timeseries", quote = FALSE)
pb <- progress_bar$new(
  format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
  total = length(T) * length(Density) * length(N) * length(M) * n_ts,
  clear = FALSE,
  width = 80
)

# PERFORMANCE: Pre-allocate tibble to maximum possible size
# This avoids costly row-binding in loops
max_rows <- length(T) * length(Density) * length(N) * length(M) * n_ts * max(M)

msar_results <- tibble::tibble(
  timesteps = integer(max_rows),
  density = numeric(max_rows),
  nodes = integer(max_rows),
  regimes = integer(max_rows),
  ts_id = integer(max_rows),
  regime_id = integer(max_rows),

  # List-columns for network matrices
  orig_Beta = vector("list", max_rows),
  est_Beta = vector("list", max_rows),
  orig_Kappa = vector("list", max_rows),
  est_Kappa = vector("list", max_rows),
  orig_Beta_ac = vector("list", max_rows),
  est_Beta_ac = vector("list", max_rows),

  # Scalar metrics
  Beta_corr = numeric(max_rows),
  Beta_sen = numeric(max_rows),
  Beta_spec = numeric(max_rows),
  MAE_Beta = numeric(max_rows),
  Beta_ac_corr = numeric(max_rows),
  Kappa_corr = numeric(max_rows),
  Kappa_sen = numeric(max_rows),
  Kappa_spec = numeric(max_rows),
  MAE_Kappa = numeric(max_rows)
)

row_idx <- 0  # Index counter for tibble filling

# Pre-allocate tibble for regime-sequence recovery (M > 1 only)
max_seq_rows <- length(T) * length(Density) * length(N) * length(M) * n_ts

sequence_results <- tibble::tibble(
  timesteps = integer(max_seq_rows),
  density = numeric(max_seq_rows),
  nodes = integer(max_seq_rows),
  regimes = integer(max_seq_rows),
  ts_id = integer(max_seq_rows),
  accuracy = numeric(max_seq_rows),
  cohens_kappa = numeric(max_seq_rows)
)
seq_idx <- 0  # Index counter for sequence_results

# markov regime-switching autoregression models with expectation maximization method
for (t in seq_along(T)) {
  for (i in seq_along(Density)) {
    for (j in seq_along(N)) {
      for (k in seq_along(M)) {
        for (l in 1:n_ts) {

          loc <-
            paste("i=", i, ", j=", j, ", k=", k, ", l=", l, sep = "")

          if (verbose)
          {
            print("", quote = FALSE)
            print("", quote = FALSE)
            print(paste("Timeseries:", loc), quote = FALSE)
            print("", quote = FALSE)
          }

          # Get complete data for current timeseries from tibble
          current_row <- Timeseries_data %>%
            filter(timesteps == T[t], density == Density[i],
                   nodes == N[j], regimes == M[k], ts_id == l)

          # Check if we found the data
          if (nrow(current_row) == 0) {
            message("No timeseries data found for condition: ",
                    "T=", T[t], " D=", Density[i], " N=", N[j], " M=", M[k], " ts=", l)
            pb$tick()
            next
          }

          # Extract timeseries data from list-column
          current_ts <- current_row$timeseries_data[[1]]

          # Extract regimes data from list-column
          current_regimes <- current_row$regime_dynamics[[1]]

          # Normalize timeseries data
          current_ts_norm <- huge.npn(current_ts, verbose = verbose)

          # Make array from timeseries data
          timesteps <- T[t]                   # Set no. of timesteps
          d <- ncol(current_ts_norm)          # Set no. of nodes
          N.samples <- 1                      # Set no. of time series
          current_ts_norm_array <- array(data = current_ts_norm,
                                         dim = c(timesteps, N.samples, d))

          # Estimate MSAR models
          result <- init_and_fit_msar_lasso(
            data = current_ts_norm_array,
            M = M[k],
            order = order,
            MaxIter = MaxIter,
            retry = 5,
            verbose = verbose
          )

          model_fit <- result[["fit"]]
          error <- result[["error"]]

          # browser()

          if (is.null(model_fit)) {
            message("ignore fit! ", "\n\tTimeseries_data: ", loc, "\n\terror: ", error)
            pb$tick()
            next
          }

         # Sys.sleep(1)

          # Extract estimated Betas (temporal network) for all regimes from MSAR model
          est_Betas <- model_fit[["theta"]][["A"]]

          # Extract estimated sigma for all regimes
          est_sigmas <- model_fit[["theta"]][["sigma"]]

          # Make Kappa from estimated sigma for all regimes (contemporaneous network)
          est_Kappas <- list()
          kappa_failed <- FALSE
          for(m in 1:M[k]) {
            est_kappa <- tryCatch(
              solve(est_sigmas[[m]]),
              error = function(e) NULL
            )
            if (is.null(est_kappa)) {
              kappa_failed <- TRUE
              break
            }
            # Kappa is undirected; symmetrize to reduce numeric asymmetry
            est_Kappas[[paste0("Regime", m)]] <- symmetrize_matrix(est_kappa)
          }
          if (kappa_failed) {
            message("ignore fit! ", "\n\tTimeseries_data: ", loc,
                    "\n\terror: singular covariance matrix in at least one regime")
            pb$tick()
            next
          }

          # Collect raw estimated Beta matrices per regime
          est_Betas_mat <- list()
          for(m in 1:M[k]) {
            est_Betas_mat[[paste0("Regime", m)]] <- est_Betas[[m]][["A1"]]
          }

          # Remove spurious edges using bootstrap stability selection or simple threshold.
          # The diagonal of Kappa is never thresholded (structurally dominant, required
          # for positive definiteness).
          if (use_bootstrap) {
            # Bootstrap stability selection
            if (verbose) {
              message("Running bootstrap stability selection...")
            }

            stability_result <- bootstrap_stability_selection(
              data_norm = current_ts_norm,
              M = M[k],
              order = order,
              MaxIter = MaxIter,
              n_bootstrap = n_bootstrap,
              block_size = min(20, floor(T[t] / 10)),
              threshold = bootstrap_threshold,
              verbose = verbose
            )

            if (!is.null(stability_result)) {
              # Apply bootstrap selection
              filtered <- apply_bootstrap_selection(est_Betas_mat, est_Kappas, stability_result)
              est_Betas_mat <- filtered$est_Betas
              est_Kappas <- filtered$est_Kappas

              if (verbose) {
                for (m in 1:M[k]) {
                  n_beta <- sum(stability_result$Beta_selected[[m]])
                  n_kappa <- sum(stability_result$Kappa_selected[[m]])
                  message("  Regime ", m, ": ", n_beta, " Beta edges, ", n_kappa, " Kappa edges selected")
                }
              }
            } else {
              # Bootstrap failed, fall back to simple threshold
              if (verbose) {
                message("Bootstrap failed, using simple threshold")
              }
              for(m in 1:M[k]) {
                for(n in 1:N[j]) {
                  for(o in 1:N[j]) {
                    if(abs(est_Betas_mat[[m]][n, o]) < min_edg_val){
                      est_Betas_mat[[m]][n, o] <- 0
                    }
                    if(n != o && abs(est_Kappas[[m]][n, o]) < min_edg_val){
                      est_Kappas[[m]][n, o] <- 0
                    }
                  }
                }
              }
            }
          } else {
            # Simple threshold (original approach)
            for(m in 1:M[k]) {
              for(n in 1:N[j]) {
                for(o in 1:N[j]) {
                  if(abs(est_Betas_mat[[m]][n, o]) < min_edg_val){
                    est_Betas_mat[[m]][n, o] <- 0
                  }
                  if(n != o && abs(est_Kappas[[m]][n, o]) < min_edg_val){
                    est_Kappas[[m]][n, o] <- 0
                  }
                }
              }
            }
          }

          # Extract original Beta
          get_org_Beta <- function(current_regimes) {
            vectors <- list()
            for (regime in names(current_regimes)) {
              vec_name <- paste(regime, "Beta", sep = "_")
              vectors[[vec_name]] <- as.vector(current_regimes[[regime]]$Beta)
            }
            return(vectors)
          }

          # Extract estimated Beta
          get_est_Beta <- function(est_Betas_mat) {
            vectors <- list()
            for (regime in names(est_Betas_mat)) {
              vec_name <- paste(regime, "Beta", sep = "_")
              vectors[[vec_name]] <- as.vector(est_Betas_mat[[regime]])
            }
            return(vectors)
          }

          # Make matrix of vectors from original and estimated Beta
          vecs_org_Beta <- get_org_Beta(current_regimes)
          vecs_est_Beta <- get_est_Beta(est_Betas_mat)

          # Check if estimated Beta has sd == 0
          if (any(sapply(vecs_est_Beta, function(vec) sd(vec) == 0))) {
            pb$tick()
            next
          }

          # Calculate correlations for each pair of original and estimated Beta
          # (this Beta-based correlation matrix is also used for regime assignment)
          cor_results <- outer(names(vecs_org_Beta),
                               names(vecs_est_Beta),
                               Vectorize(function(n1, n2) {
                                 cor(vecs_org_Beta[[n1]],
                                     vecs_est_Beta[[n2]])
                               }
                               )
          )


          # Assign original and estimated regimes based on highest correlations of Beta
          if (M[k] > 1) {
            assigned_regimes <- assign_regimes(cor_results)

            # Process each regime pair
            for (m in 1:nrow(assigned_regimes)) {
              row_idx <- row_idx + 1

              orig_Beta <- current_regimes[[assigned_regimes[[m, 1]]]][["Beta"]]
              est_Beta <- est_Betas_mat[[assigned_regimes[[m, 2]]]]
              Beta_senspec <- senspec(orig_Beta, est_Beta)
              MAE_Beta <- calculate_MAE(orig_Beta, est_Beta)

              est_Beta_ac <- average_controllability(est_Beta)
              orig_Beta_ac <- current_regimes[[assigned_regimes[[m, 1]]]][["Beta_ac"]]
              cor_Beta_ac <- cor(orig_Beta_ac, est_Beta_ac, method = "pearson")

              orig_Kappa <- current_regimes[[assigned_regimes[[m, 1]]]][["kappa"]]
              est_Kappa <- est_Kappas[[assigned_regimes[[m, 2]]]]
              Kappa_senspec <- senspec_upper_tri(orig_Kappa, est_Kappa, diag = FALSE)
              cor_Kappa <- cor(
                vectorize_upper_tri(orig_Kappa, diag = FALSE),
                vectorize_upper_tri(est_Kappa, diag = FALSE),
                method = "pearson"
              )
              MAE_Kappa <- mae_upper_tri(orig_Kappa, est_Kappa, diag = FALSE)

              # Store in tibble
              msar_results$timesteps[row_idx] <- T[t]
              msar_results$density[row_idx] <- Density[i]
              msar_results$nodes[row_idx] <- N[j]
              msar_results$regimes[row_idx] <- M[k]
              msar_results$ts_id[row_idx] <- l
              msar_results$regime_id[row_idx] <- m

              msar_results$orig_Beta[[row_idx]] <- orig_Beta
              msar_results$est_Beta[[row_idx]] <- est_Beta
              msar_results$orig_Kappa[[row_idx]] <- orig_Kappa
              msar_results$est_Kappa[[row_idx]] <- est_Kappa
              msar_results$orig_Beta_ac[[row_idx]] <- orig_Beta_ac
              msar_results$est_Beta_ac[[row_idx]] <- est_Beta_ac

              msar_results$Beta_corr[row_idx] <- assigned_regimes[[m, 3]]
              msar_results$Beta_sen[row_idx] <- Beta_senspec[["sensitivity"]]
              msar_results$Beta_spec[row_idx] <- Beta_senspec[["specificity"]]
              msar_results$MAE_Beta[row_idx] <- MAE_Beta
              msar_results$Beta_ac_corr[row_idx] <- cor_Beta_ac
              msar_results$Kappa_corr[row_idx] <- cor_Kappa
              msar_results$Kappa_sen[row_idx] <- Kappa_senspec[["sensitivity"]]
              msar_results$Kappa_spec[row_idx] <- Kappa_senspec[["specificity"]]
              msar_results$MAE_Kappa[row_idx] <- MAE_Kappa
            }

            # ---- Regime-sequence recovery (M > 1 only) ----
            # True sequence: regime_sequence has length (totTime - 1); the
            # trimmed timeseries data keeps only the last T rows (after warmup
            # removal in generate_timeseries.R). The order-1 MSAR fit further
            # drops the very first of those T rows (used only as the initial
            # AR lag, no regime probability estimated for it), so the smoothed
            # probabilities have T - 1 rows. We align by taking the last T
            # entries of the full true sequence, then dropping the first of
            # those to match.
            true_seq_full <- current_row$regime_sequence[[1]]
            true_seq <- utils::tail(true_seq_full, T[t])[-1]

            # Estimated hard sequence from smoothed probabilities, relabeled
            # using the Beta-based regime mapping (assigned_regimes: col 1 =
            # true label, col 2 = estimated label that maps to it)
            est_seq_raw <- get_hard_regime_sequence(model_fit[["smoothedprob"]])
            label_map <- setNames(assigned_regimes[, 1], assigned_regimes[, 2])
            est_seq_mapped <- as.integer(label_map[as.character(est_seq_raw)])

            if (length(true_seq) == length(est_seq_mapped)) {
              seq_idx <- seq_idx + 1
              sequence_results$timesteps[seq_idx] <- T[t]
              sequence_results$density[seq_idx] <- Density[i]
              sequence_results$nodes[seq_idx] <- N[j]
              sequence_results$regimes[seq_idx] <- M[k]
              sequence_results$ts_id[seq_idx] <- l
              sequence_results$accuracy[seq_idx] <- mean(true_seq == est_seq_mapped)
              sequence_results$cohens_kappa[seq_idx] <- cohens_kappa_manual(true_seq, est_seq_mapped)
            } else {
              message("Skipping sequence recovery for ", loc,
                      ": length mismatch (true = ", length(true_seq),
                      ", est = ", length(est_seq_mapped), ")")
            }

          } else {
            # Single regime case (no sequence recovery: trivial with M = 1)
            row_idx <- row_idx + 1

            orig_Beta <- current_regimes[["Regime1"]][["Beta"]]
            est_Beta <- est_Betas_mat[["Regime1"]]
            Beta_senspec <- senspec(orig_Beta, est_Beta)
            cor_Beta <- cor(as.vector(orig_Beta), as.vector(est_Beta), method = "pearson")
            MAE_Beta <- calculate_MAE(orig_Beta, est_Beta)

            est_Beta_ac <- average_controllability(est_Beta)
            orig_Beta_ac <- current_regimes[["Regime1"]][["Beta_ac"]]
            cor_Beta_ac <- cor(orig_Beta_ac, est_Beta_ac, method = "pearson")

            orig_Kappa <- current_regimes[["Regime1"]][["kappa"]]
            est_Kappa <- est_Kappas[["Regime1"]]
            cor_Kappa <- cor(
              vectorize_upper_tri(orig_Kappa, diag = FALSE),
              vectorize_upper_tri(est_Kappa, diag = FALSE),
              method = "pearson"
            )
            Kappa_senspec <- senspec_upper_tri(orig_Kappa, est_Kappa, diag = FALSE)
            MAE_Kappa <- mae_upper_tri(orig_Kappa, est_Kappa, diag = FALSE)

            # Store in tibble
            msar_results$timesteps[row_idx] <- T[t]
            msar_results$density[row_idx] <- Density[i]
            msar_results$nodes[row_idx] <- N[j]
            msar_results$regimes[row_idx] <- M[k]
            msar_results$ts_id[row_idx] <- l
            msar_results$regime_id[row_idx] <- 1

            msar_results$orig_Beta[[row_idx]] <- orig_Beta
            msar_results$est_Beta[[row_idx]] <- est_Beta
            msar_results$orig_Kappa[[row_idx]] <- orig_Kappa
            msar_results$est_Kappa[[row_idx]] <- est_Kappa
            msar_results$orig_Beta_ac[[row_idx]] <- orig_Beta_ac
            msar_results$est_Beta_ac[[row_idx]] <- est_Beta_ac

            msar_results$Beta_corr[row_idx] <- cor_Beta
            msar_results$Beta_sen[row_idx] <- Beta_senspec[["sensitivity"]]
            msar_results$Beta_spec[row_idx] <- Beta_senspec[["specificity"]]
            msar_results$MAE_Beta[row_idx] <- MAE_Beta
            msar_results$Beta_ac_corr[row_idx] <- cor_Beta_ac
            msar_results$Kappa_corr[row_idx] <- cor_Kappa
            msar_results$Kappa_sen[row_idx] <- Kappa_senspec[["sensitivity"]]
            msar_results$Kappa_spec[row_idx] <- Kappa_senspec[["specificity"]]
            msar_results$MAE_Kappa[row_idx] <- MAE_Kappa
          }

          # Update progressbar
          pb$tick()
        }
      }
    }
  }
}

# Trim to actual size (remove pre-allocated empty rows)
msar_results <- msar_results[1:row_idx, ]
sequence_results <- sequence_results[1:seq_idx, ]

# Add S3 class
class(msar_results) <- c("msar_results", class(msar_results))

# Attach regime-sequence recovery as an attribute (see @return for rationale)
attr(msar_results, "sequence_results") <- sequence_results

return(msar_results)

}


#' Get Summary Statistics from MSAR Results
#'
#' Computes summary statistics (mean, sd, median, min, max, N) for each metric
#' grouped by experimental conditions.
#'
#' @param msar_results An msar_results object from \code{estimate_MSAR}
#' @param group_by Character vector of grouping variables.
#'   Default: c("timesteps", "density", "nodes", "regimes")
#'
#' @return A tibble with summary statistics for each metric
#'
#' @export
get_stats <- function(msar_results,
                      group_by = c("timesteps", "density", "nodes", "regimes")) {

  if (!inherits(msar_results, "msar_results")) {
    warning("Input is not an msar_results object. Treating as tibble.")
  }

  msar_results %>%
    dplyr::group_by(across(all_of(group_by))) %>%
    dplyr::summarise(
      # Beta correlations
      Beta_corr_mean = mean(Beta_corr, na.rm = TRUE),
      Beta_corr_sd = sd(Beta_corr, na.rm = TRUE),
      Beta_corr_median = median(Beta_corr, na.rm = TRUE),
      Beta_corr_min = min(Beta_corr, na.rm = TRUE),
      Beta_corr_max = max(Beta_corr, na.rm = TRUE),

      # Beta MAE
      Beta_MAE_mean = mean(MAE_Beta, na.rm = TRUE),
      Beta_MAE_sd = sd(MAE_Beta, na.rm = TRUE),
      Beta_MAE_median = median(MAE_Beta, na.rm = TRUE),
      Beta_MAE_min = min(MAE_Beta, na.rm = TRUE),
      Beta_MAE_max = max(MAE_Beta, na.rm = TRUE),

      # Beta sensitivity
      Beta_sen_mean = mean(Beta_sen, na.rm = TRUE),
      Beta_sen_sd = sd(Beta_sen, na.rm = TRUE),
      Beta_sen_median = median(Beta_sen, na.rm = TRUE),
      Beta_sen_min = min(Beta_sen, na.rm = TRUE),
      Beta_sen_max = max(Beta_sen, na.rm = TRUE),

      # Beta specificity
      Beta_spec_mean = mean(Beta_spec, na.rm = TRUE),
      Beta_spec_sd = sd(Beta_spec, na.rm = TRUE),
      Beta_spec_median = median(Beta_spec, na.rm = TRUE),
      Beta_spec_min = min(Beta_spec, na.rm = TRUE),
      Beta_spec_max = max(Beta_spec, na.rm = TRUE),

      # Beta AC correlations
      Beta_ac_corr_mean = mean(Beta_ac_corr, na.rm = TRUE),
      Beta_ac_corr_sd = sd(Beta_ac_corr, na.rm = TRUE),
      Beta_ac_corr_median = median(Beta_ac_corr, na.rm = TRUE),
      Beta_ac_corr_min = min(Beta_ac_corr, na.rm = TRUE),
      Beta_ac_corr_max = max(Beta_ac_corr, na.rm = TRUE),

      # Kappa correlations
      Kappa_corr_mean = mean(Kappa_corr, na.rm = TRUE),
      Kappa_corr_sd = sd(Kappa_corr, na.rm = TRUE),
      Kappa_corr_median = median(Kappa_corr, na.rm = TRUE),
      Kappa_corr_min = min(Kappa_corr, na.rm = TRUE),
      Kappa_corr_max = max(Kappa_corr, na.rm = TRUE),

      # Kappa MAE
      Kappa_MAE_mean = mean(MAE_Kappa, na.rm = TRUE),
      Kappa_MAE_sd = sd(MAE_Kappa, na.rm = TRUE),
      Kappa_MAE_median = median(MAE_Kappa, na.rm = TRUE),
      Kappa_MAE_min = min(MAE_Kappa, na.rm = TRUE),
      Kappa_MAE_max = max(MAE_Kappa, na.rm = TRUE),

      # Kappa sensitivity
      Kappa_sen_mean = mean(Kappa_sen, na.rm = TRUE),
      Kappa_sen_sd = sd(Kappa_sen, na.rm = TRUE),
      Kappa_sen_median = median(Kappa_sen, na.rm = TRUE),
      Kappa_sen_min = min(Kappa_sen, na.rm = TRUE),
      Kappa_sen_max = max(Kappa_sen, na.rm = TRUE),

      # Kappa specificity
      Kappa_spec_mean = mean(Kappa_spec, na.rm = TRUE),
      Kappa_spec_sd = sd(Kappa_spec, na.rm = TRUE),
      Kappa_spec_median = median(Kappa_spec, na.rm = TRUE),
      Kappa_spec_min = min(Kappa_spec, na.rm = TRUE),
      Kappa_spec_max = max(Kappa_spec, na.rm = TRUE),

      # Sample size
      N = n(),

      .groups = "drop"
    )
}


#' Print Method for msar_results
#'
#' @param x An msar_results object
#' @param ... Additional arguments (unused)
#'
#' @export
print.msar_results <- function(x, ...) {
  cat("MSAR Results\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Total observations: %d\n", nrow(x)))
  cat(sprintf("Conditions tested:\n"))
  cat(sprintf("  Timesteps: %s\n", paste(unique(x$timesteps), collapse = ", ")))
  cat(sprintf("  Density: %s\n", paste(unique(x$density), collapse = ", ")))
  cat(sprintf("  Nodes: %s\n", paste(unique(x$nodes), collapse = ", ")))
  cat(sprintf("  Regimes: %s\n", paste(unique(x$regimes), collapse = ", ")))
  cat(sprintf("  Time series per condition: %d\n", max(x$ts_id)))
  cat("════════════════════════════════════════════════════════════════\n")
  cat("\nSummary of Beta Correlations:\n")
  print(summary(x$Beta_corr))
  cat("\nUse get_stats() for detailed summary statistics.\n")
  cat("Use dplyr::filter() to subset by conditions.\n")
  NextMethod()
}
