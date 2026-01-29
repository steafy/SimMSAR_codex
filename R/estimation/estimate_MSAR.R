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
#'   Only used when use_bootstrap = FALSE.
#' @param Timeseries_data Tibble. Output from \code{generate_timeseries} containing
#'   simulated data and true dynamics in tibble format with list-columns.
#' @param use_bootstrap Logical. If TRUE, use bootstrap stability selection instead
#'   of simple thresholding. Default: FALSE. Note: This is computationally expensive.
#' @param n_bootstrap Integer. Number of bootstrap samples when use_bootstrap = TRUE.
#'   Default: 50. Higher values give more stable selection but take longer.
#' @param bootstrap_threshold Numeric. Selection threshold (0-1) for bootstrap.
#'   Edges selected in > threshold proportion of bootstraps are kept. Default: 0.6.
#'
#' @return A tibble (msar_results object) with the following columns:
#'   \describe{
#'     \item{timesteps}{Number of time steps in the series}
#'     \item{density}{Edge density of the network}
#'     \item{nodes}{Number of nodes in the network}
#'     \item{regimes}{Number of regimes in the model}
#'     \item{ts_id}{Time series ID (1 to n_ts)}
#'     \item{regime_id}{Regime ID within the time series}
#'     \item{orig_Wtemp}{True temporal network (list-column)}
#'     \item{est_Wtemp}{Estimated temporal network (list-column)}
#'     \item{orig_Wcont}{True contemporaneous network (list-column)}
#'     \item{est_Wcont}{Estimated contemporaneous network (list-column)}
#'     \item{orig_Wtemp_ac}{True temporal network average controllability (list-column)}
#'     \item{est_Wtemp_ac}{Estimated temporal network average controllability (list-column)}
#'     \item{orig_Wcont_ac}{True contemporaneous network average controllability (list-column)}
#'     \item{est_Wcont_ac}{Estimated contemporaneous network average controllability (list-column)}
#'     \item{Wtemp_corr}{Correlation between true and estimated Wtemp}
#'     \item{Wtemp_sen}{Sensitivity for Wtemp edge detection}
#'     \item{Wtemp_spec}{Specificity for Wtemp edge detection}
#'     \item{MAE_Wtemp}{Mean absolute error for Wtemp}
#'     \item{Wtemp_ac_corr}{Correlation of average controllability for Wtemp}
#'     \item{Wcont_corr}{Correlation between true and estimated Wcont}
#'     \item{Wcont_sen}{Sensitivity for Wcont edge detection}
#'     \item{Wcont_spec}{Specificity for Wcont edge detection}
#'     \item{MAE_Wcont}{Mean absolute error for Wcont}
#'     \item{Wcont_ac_corr}{Correlation of average controllability for Wcont}
#'   }
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
#'      \item Computes Kappa = Sigma^{-1} (precision matrix)
#'      \item Derives Wcont from Kappa: Wcont[i,j] = -Kappa[i,j] / sqrt(Kappa[i,i] * Kappa[j,j])
#'      \item Derives Wtemp from Beta, Sigma, Kappa
#'      \item Applies min_edg_val threshold to estimated networks
#'    }
#'
#' 4. **Regime Assignment**: For multi-regime models, matches estimated to true
#'    regimes based on maximum Wtemp correlations using \code{assign_regimes}
#'
#' 5. **Comparison Metrics**: Computes for each regime pair:
#'    \itemize{
#'      \item Correlation: Pearson correlation between vectorized networks
#'      \item Sensitivity/Specificity: Edge detection accuracy
#'      \item MAE: Mean absolute error of edge weights
#'      \item Controllability correlation: For network control metrics
#'    }
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
#' # Access results (now a tibble!)
#' results %>% filter(timesteps == 1000, density == 0.3)
#'
#' # Get summary statistics
#' get_stats(results)
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
  orig_Wtemp = vector("list", max_rows),
  est_Wtemp = vector("list", max_rows),
  orig_Wcont = vector("list", max_rows),
  est_Wcont = vector("list", max_rows),
  orig_Wtemp_ac = vector("list", max_rows),
  est_Wtemp_ac = vector("list", max_rows),
  orig_Wcont_ac = vector("list", max_rows),
  est_Wcont_ac = vector("list", max_rows),

  # Scalar metrics
  Wtemp_corr = numeric(max_rows),
  Wtemp_sen = numeric(max_rows),
  Wtemp_spec = numeric(max_rows),
  MAE_Wtemp = numeric(max_rows),
  Wtemp_ac_corr = numeric(max_rows),
  Wcont_corr = numeric(max_rows),
  Wcont_sen = numeric(max_rows),
  Wcont_spec = numeric(max_rows),
  MAE_Wcont = numeric(max_rows),
  Wcont_ac_corr = numeric(max_rows)
)

row_idx <- 0  # Index counter for tibble filling

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
            print(paste("fit parms:", app_parms), quote = FALSE)
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

          # Extract estimated Betas for all regimes from MSAR model
          est_Betas <- model_fit[["theta"]][["A"]]

          # Extract estimated sigma for all regimes
          est_sigmas <- model_fit[["theta"]][["sigma"]]

          # Make kappa from estimated sigma for all regimes
          est_kappas <- list()
          for(m in 1:M[k]) {
            est_kappas[[paste0("Regime", m)]] <- solve(est_sigmas[[m]])
          }


          # Make estimated Wcont from estimated kappa
          est_Wconts <- list()
          for(m in 1:M[k]) {
            est_Wcont <- matrix(0, nrow = N[j], ncol = N[j])
            est_kappa <- est_kappas[[m]]
            for(n in 1:N[j]) {
              for(o in 1:N[j]) {
                if(n != o) {
                  est_Wcont[n, o] <- -est_kappa[n, o] / sqrt(est_kappa[n, n] * est_kappa[o, o])
                  est_Wconts[[paste0("Regime", m)]] <- est_Wcont
                }
              }
            }
          }

          # Make estimated Wtemps from estimates of Beta, kappa and sigma
          est_Wtemps <- list()
          for(m in 1:M[k]) {
            est_Wtemp <- matrix(0, nrow = N[j], ncol = N[j])
            est_Beta <- est_Betas[[m]][["A1"]]
            est_kappa <- est_kappas[[m]]
            est_sigma <- est_sigmas[[m]]
            for(n in 1:N[j]) {
              for(o in 1:N[j]) { # raus: if(n != o) {
                # if(n != o) {
                #est_Wtemp[n, o] <- est_Beta[n, o] / sqrt(est_sigma[i, i] %*% est_kappa[j ,j] + est_Beta[i, j]^2)
                est_Wtemp[n, o] <- est_Beta[n, o] / sqrt(est_sigma[n, n] %*% est_kappa[o, o] + est_Beta[n, o]^2)
                est_Wtemps[[paste0("Regime", m)]] <- est_Wtemp
                # }
              }
            }
          }

          # Remove spurious edges using bootstrap stability selection or simple threshold
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
              filtered <- apply_bootstrap_selection(est_Wtemps, est_Wconts, stability_result)
              est_Wtemps <- filtered$est_Wtemps
              est_Wconts <- filtered$est_Wconts

              if (verbose) {
                for (m in 1:M[k]) {
                  n_wtemp <- sum(stability_result$Wtemp_selected[[m]])
                  n_wcont <- sum(stability_result$Wcont_selected[[m]])
                  message("  Regime ", m, ": ", n_wtemp, " Wtemp edges, ", n_wcont, " Wcont edges selected")
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
                    if(abs(est_Wtemps[[m]][n, o]) < min_edg_val){
                      est_Wtemps[[m]][n, o] <- 0
                    }
                    if(abs(est_Wconts[[m]][n, o]) < min_edg_val){
                      est_Wconts[[m]][n, o] <- 0
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
                  if(abs(est_Wtemps[[m]][n, o]) < min_edg_val){
                    est_Wtemps[[m]][n, o] <- 0
                  }
                  if(abs(est_Wconts[[m]][n, o]) < min_edg_val){
                    est_Wconts[[m]][n, o] <- 0
                  }
                }
              }
            }
          }

          # Extract original Wtemp
          get_org_Wtemp <- function(current_regimes) {
            vectors <- list()
            for (regime in names(current_regimes)) {
              vec_name <- paste(regime, "Wtemp", sep = "_")
              vectors[[vec_name]] <- as.vector(current_regimes[[regime]]$Wtemp)
            }
            return(vectors)
          }

          # Extract estimated Wtemp
          get_est_Wtemp  <- function(model_fit) {  # vorher function(model_fit)
            vectors <- list()
            for (regime in names(est_Wtemps)) {
              vec_name <- paste(regime, "Wtemp", sep = "_")
              vectors[[vec_name]] <- as.vector(est_Wtemps[[regime]])
            }
            return(vectors)
          }

          # Make matrix of vectors from original and estimated Wtemp
          vecs_org_Wtemp <- get_org_Wtemp(current_regimes)
          vecs_est_Wtemp <- get_est_Wtemp(model_fit)

          # Check if estimated Wtemp has sd == 0
          if (any(sapply(vecs_est_Wtemp,function(vec) sd(vec) == 0))) {
            pb$tick()
            next
          }

          # Calculate correlations for each pair of original and estimated Wtemp
          cor_results <- outer(names(vecs_org_Wtemp),
                               names(vecs_est_Wtemp),
                               Vectorize(function(n1, n2) {
                                 cor(vecs_org_Wtemp[[n1]],
                                     vecs_est_Wtemp[[n2]])
                               }
                               )
          )


          # Assign original and estimated regimes based on highest correlations of Wtemp
          if (M[k] > 1) {
            assigned_regimes <- assign_regimes(cor_results)

            # Process each regime pair
            for (m in 1:nrow(assigned_regimes)) {
              row_idx <- row_idx + 1

              orig_Wtemp <- current_regimes[[assigned_regimes[[m, 1]]]][["Wtemp"]]
              est_Wtemp <- est_Wtemps[[assigned_regimes[[m, 2]]]]
              Wtemp_senspec <- senspec(orig_Wtemp, est_Wtemp)
              MAE_Wtemp <- calculate_MAE(orig_Wtemp, est_Wtemp)

              est_Wtemp_ac <- ave_control_centrality(est_Wtemp)
              orig_Wtemp_ac <- current_regimes[[assigned_regimes[[m, 1]]]][["Wtemp_ac"]]
              cor_Wtemp_ac <- cor(orig_Wtemp_ac, est_Wtemp_ac, method = "pearson")

              orig_Wcont <- current_regimes[[assigned_regimes[[m, 1]]]][["Wcont"]]
              est_Wcont <- est_Wconts[[assigned_regimes[[m, 2]]]]
              Wcont_senspec <- senspec(orig_Wcont, est_Wcont)
              cor_Wcont <- cor(as.vector(orig_Wcont), as.vector(est_Wcont), method = "pearson")
              MAE_Wcont <- calculate_MAE(orig_Wcont, est_Wcont)

              est_Wcont_ac <- ave_control_centrality(est_Wcont)
              orig_Wcont_ac <- current_regimes[[assigned_regimes[[m, 1]]]][["Wcont_ac"]]
              cor_Wcont_ac <- cor(orig_Wcont_ac, est_Wcont_ac, method = "pearson")

              # Store in tibble
              msar_results$timesteps[row_idx] <- T[t]
              msar_results$density[row_idx] <- Density[i]
              msar_results$nodes[row_idx] <- N[j]
              msar_results$regimes[row_idx] <- M[k]
              msar_results$ts_id[row_idx] <- l
              msar_results$regime_id[row_idx] <- m

              msar_results$orig_Wtemp[[row_idx]] <- orig_Wtemp
              msar_results$est_Wtemp[[row_idx]] <- est_Wtemp
              msar_results$orig_Wcont[[row_idx]] <- orig_Wcont
              msar_results$est_Wcont[[row_idx]] <- est_Wcont
              msar_results$orig_Wtemp_ac[[row_idx]] <- orig_Wtemp_ac
              msar_results$est_Wtemp_ac[[row_idx]] <- est_Wtemp_ac
              msar_results$orig_Wcont_ac[[row_idx]] <- orig_Wcont_ac
              msar_results$est_Wcont_ac[[row_idx]] <- est_Wcont_ac

              msar_results$Wtemp_corr[row_idx] <- assigned_regimes[[m, 3]]
              msar_results$Wtemp_sen[row_idx] <- Wtemp_senspec[["sensitivity"]]
              msar_results$Wtemp_spec[row_idx] <- Wtemp_senspec[["specificity"]]
              msar_results$MAE_Wtemp[row_idx] <- MAE_Wtemp
              msar_results$Wtemp_ac_corr[row_idx] <- cor_Wtemp_ac
              msar_results$Wcont_corr[row_idx] <- cor_Wcont
              msar_results$Wcont_sen[row_idx] <- Wcont_senspec[["sensitivity"]]
              msar_results$Wcont_spec[row_idx] <- Wcont_senspec[["specificity"]]
              msar_results$MAE_Wcont[row_idx] <- MAE_Wcont
              msar_results$Wcont_ac_corr[row_idx] <- cor_Wcont_ac
            }
          } else {
            # Single regime case
            row_idx <- row_idx + 1

            orig_Wtemp <- current_regimes[["Regime1"]][["Wtemp"]]
            est_Wtemp <- est_Wtemps[["Regime1"]]
            Wtemp_senspec <- senspec(orig_Wtemp, est_Wtemp)
            cor_Wtemp <- cor(as.vector(orig_Wtemp), as.vector(est_Wtemp), method = "pearson")
            MAE_Wtemp <- calculate_MAE(orig_Wtemp, est_Wtemp)

            est_Wtemp_ac <- ave_control_centrality(est_Wtemp)
            orig_Wtemp_ac <- current_regimes[["Regime1"]][["Wtemp_ac"]]
            cor_Wtemp_ac <- cor(orig_Wtemp_ac, est_Wtemp_ac, method = "pearson")

            orig_Wcont <- current_regimes[["Regime1"]][["Wcont"]]
            est_Wcont <- est_Wconts[["Regime1"]]
            cor_Wcont <- cor(as.vector(orig_Wcont), as.vector(est_Wcont), method = "pearson")
            Wcont_senspec <- senspec(orig_Wcont, est_Wcont)

            est_Wcont_ac <- ave_control_centrality(est_Wcont)
            orig_Wcont_ac <- current_regimes[["Regime1"]][["Wcont_ac"]]
            cor_Wcont_ac <- cor(orig_Wcont_ac, est_Wcont_ac, method = "pearson")
            MAE_Wcont <- calculate_MAE(orig_Wcont, est_Wcont)

            # Store in tibble
            msar_results$timesteps[row_idx] <- T[t]
            msar_results$density[row_idx] <- Density[i]
            msar_results$nodes[row_idx] <- N[j]
            msar_results$regimes[row_idx] <- M[k]
            msar_results$ts_id[row_idx] <- l
            msar_results$regime_id[row_idx] <- 1

            msar_results$orig_Wtemp[[row_idx]] <- orig_Wtemp
            msar_results$est_Wtemp[[row_idx]] <- est_Wtemp
            msar_results$orig_Wcont[[row_idx]] <- orig_Wcont
            msar_results$est_Wcont[[row_idx]] <- est_Wcont
            msar_results$orig_Wtemp_ac[[row_idx]] <- orig_Wtemp_ac
            msar_results$est_Wtemp_ac[[row_idx]] <- est_Wtemp_ac
            msar_results$orig_Wcont_ac[[row_idx]] <- orig_Wcont_ac
            msar_results$est_Wcont_ac[[row_idx]] <- est_Wcont_ac

            msar_results$Wtemp_corr[row_idx] <- cor_Wtemp
            msar_results$Wtemp_sen[row_idx] <- Wtemp_senspec[["sensitivity"]]
            msar_results$Wtemp_spec[row_idx] <- Wtemp_senspec[["specificity"]]
            msar_results$MAE_Wtemp[row_idx] <- MAE_Wtemp
            msar_results$Wtemp_ac_corr[row_idx] <- cor_Wtemp_ac
            msar_results$Wcont_corr[row_idx] <- cor_Wcont
            msar_results$Wcont_sen[row_idx] <- Wcont_senspec[["sensitivity"]]
            msar_results$Wcont_spec[row_idx] <- Wcont_senspec[["specificity"]]
            msar_results$MAE_Wcont[row_idx] <- MAE_Wcont
            msar_results$Wcont_ac_corr[row_idx] <- cor_Wcont_ac
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

# Add S3 class
class(msar_results) <- c("msar_results", class(msar_results))

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
      # Wtemp correlations
      Wtemp_corr_mean = mean(Wtemp_corr, na.rm = TRUE),
      Wtemp_corr_sd = sd(Wtemp_corr, na.rm = TRUE),
      Wtemp_corr_median = median(Wtemp_corr, na.rm = TRUE),
      Wtemp_corr_min = min(Wtemp_corr, na.rm = TRUE),
      Wtemp_corr_max = max(Wtemp_corr, na.rm = TRUE),

      # Wtemp MAE
      Wtemp_MAE_mean = mean(MAE_Wtemp, na.rm = TRUE),
      Wtemp_MAE_sd = sd(MAE_Wtemp, na.rm = TRUE),
      Wtemp_MAE_median = median(MAE_Wtemp, na.rm = TRUE),
      Wtemp_MAE_min = min(MAE_Wtemp, na.rm = TRUE),
      Wtemp_MAE_max = max(MAE_Wtemp, na.rm = TRUE),

      # Wtemp sensitivity
      Wtemp_sen_mean = mean(Wtemp_sen, na.rm = TRUE),
      Wtemp_sen_sd = sd(Wtemp_sen, na.rm = TRUE),
      Wtemp_sen_median = median(Wtemp_sen, na.rm = TRUE),
      Wtemp_sen_min = min(Wtemp_sen, na.rm = TRUE),
      Wtemp_sen_max = max(Wtemp_sen, na.rm = TRUE),

      # Wtemp specificity
      Wtemp_spec_mean = mean(Wtemp_spec, na.rm = TRUE),
      Wtemp_spec_sd = sd(Wtemp_spec, na.rm = TRUE),
      Wtemp_spec_median = median(Wtemp_spec, na.rm = TRUE),
      Wtemp_spec_min = min(Wtemp_spec, na.rm = TRUE),
      Wtemp_spec_max = max(Wtemp_spec, na.rm = TRUE),

      # Wtemp AC correlations
      Wtemp_ac_corr_mean = mean(Wtemp_ac_corr, na.rm = TRUE),
      Wtemp_ac_corr_sd = sd(Wtemp_ac_corr, na.rm = TRUE),
      Wtemp_ac_corr_median = median(Wtemp_ac_corr, na.rm = TRUE),
      Wtemp_ac_corr_min = min(Wtemp_ac_corr, na.rm = TRUE),
      Wtemp_ac_corr_max = max(Wtemp_ac_corr, na.rm = TRUE),

      # Wcont correlations
      Wcont_corr_mean = mean(Wcont_corr, na.rm = TRUE),
      Wcont_corr_sd = sd(Wcont_corr, na.rm = TRUE),
      Wcont_corr_median = median(Wcont_corr, na.rm = TRUE),
      Wcont_corr_min = min(Wcont_corr, na.rm = TRUE),
      Wcont_corr_max = max(Wcont_corr, na.rm = TRUE),

      # Wcont MAE
      Wcont_MAE_mean = mean(MAE_Wcont, na.rm = TRUE),
      Wcont_MAE_sd = sd(MAE_Wcont, na.rm = TRUE),
      Wcont_MAE_median = median(MAE_Wcont, na.rm = TRUE),
      Wcont_MAE_min = min(MAE_Wcont, na.rm = TRUE),
      Wcont_MAE_max = max(MAE_Wcont, na.rm = TRUE),

      # Wcont sensitivity
      Wcont_sen_mean = mean(Wcont_sen, na.rm = TRUE),
      Wcont_sen_sd = sd(Wcont_sen, na.rm = TRUE),
      Wcont_sen_median = median(Wcont_sen, na.rm = TRUE),
      Wcont_sen_min = min(Wcont_sen, na.rm = TRUE),
      Wcont_sen_max = max(Wcont_sen, na.rm = TRUE),

      # Wcont specificity
      Wcont_spec_mean = mean(Wcont_spec, na.rm = TRUE),
      Wcont_spec_sd = sd(Wcont_spec, na.rm = TRUE),
      Wcont_spec_median = median(Wcont_spec, na.rm = TRUE),
      Wcont_spec_min = min(Wcont_spec, na.rm = TRUE),
      Wcont_spec_max = max(Wcont_spec, na.rm = TRUE),

      # Wcont AC correlations
      Wcont_ac_corr_mean = mean(Wcont_ac_corr, na.rm = TRUE),
      Wcont_ac_corr_sd = sd(Wcont_ac_corr, na.rm = TRUE),
      Wcont_ac_corr_median = median(Wcont_ac_corr, na.rm = TRUE),
      Wcont_ac_corr_min = min(Wcont_ac_corr, na.rm = TRUE),
      Wcont_ac_corr_max = max(Wcont_ac_corr, na.rm = TRUE),

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
  cat("\nSummary of Wtemp Correlations:\n")
  print(summary(x$Wtemp_corr))
  cat("\nUse get_stats() for detailed summary statistics.\n")
  cat("Use dplyr::filter() to subset by conditions.\n")
  NextMethod()
}

