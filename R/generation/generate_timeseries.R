#' Generate MSAR Time Series Data
#'
#' Generates multivariate time series data from Markov-Switching Autoregressive (MSAR)
#' models with network dynamics. Creates random network dynamics for each regime and
#' simulates time series using regime-switching behavior controlled by a transition matrix.
#'
#' @param Density Numeric vector of network edge densities (0 to 1).
#'   Defines the proportion of possible edges that exist in the network.
#' @param N Integer vector of node counts. Number of variables/nodes in the network.
#' @param M Integer vector of regime counts. Number of distinct network states.
#' @param T Integer vector of time steps. Length of each generated time series (excluding warmup).
#' @param n_ts Integer. Number of time series to generate per condition.
#' @param warmup Integer. Number of initial time steps to discard (allows dynamics to stabilize).
#' @param totTime Integer vector. Total time steps including warmup (typically T + warmup).
#' @param mean_rep Numeric. Mean length of regime repetitions in sequence (currently unused).
#' @param sd_rep Numeric. Standard deviation for regime repetition lengths (currently unused).
#' @param min_edg_val Numeric. Minimum absolute edge weight. Edges below this are set to zero.
#' @param max_edg_val Numeric. Maximum absolute edge weight.
#' @param remain_lower Numeric. Lower bound for probability of remaining in same regime (0 to 1).
#' @param remain_upper Numeric. Upper bound for probability of remaining in same regime (0 to 1).
#'
#' @return A nested list structure containing:
#'   \describe{
#'     \item{Timeseries}{Matrix of simulated time series data (T × N)}
#'     \item{Regime_sequence}{Vector of regime indices for each time step}
#'     \item{Regime_dynamics}{List of true network dynamics for each regime, including:
#'       \itemize{
#'         \item mu: Mean vector
#'         \item W_temp: Temporal network (lag-1 effects)
#'         \item Beta: Autoregressive coefficient matrix
#'         \item sigma: Residual covariance matrix
#'         \item kappa: Precision matrix
#'         \item W_cont: Contemporaneous network (partial correlations)
#'       }
#'     }
#'     \item{Transmat}{Transition probability matrix between regimes}
#'   }
#'
#'   The list is nested by: Timesteps > Density > Nodes > Regimes > Timeseries_i
#'
#' @details
#' The function operates in two stages:
#'
#' 1. **Generate Network Dynamics**: For each combination of parameters, generates
#'    network dynamics (temporal and contemporaneous networks, covariance matrices)
#'    using \code{generate_netdyn}. Repeats generation until minimum edge value
#'    constraints are satisfied.
#'
#' 2. **Generate Time Series**: Simulates VAR(1) time series data with regime switching:
#'    \itemize{
#'      \item Initial state sampled from uniform(0, 5)
#'      \item Regime transitions follow a Markov process via transition matrix
#'      \item State evolution: X[t] = mu + Beta(X[t-1] - mu) + epsilon
#'      \item epsilon ~ MVN(0, sigma) with regime-specific parameters
#'    }
#'
#' Progress bars display generation status for both stages.
#'
#' @note
#' Dependencies are loaded centrally via R/dependencies.R
#' Required packages: netcontrol, graphicalVAR, abind, mvtnorm, progress
#'
#' @seealso
#' \code{\link{generate_netdyn}} for network dynamics generation
#' \code{\link{generate_transmat}} for transition matrix generation
#'
#' @examples
#' \dontrun{
#' # Generate simple 2-regime time series
#' ts_data <- generate_timeseries(
#'   Density = 0.3,
#'   N = 4,
#'   M = 2,
#'   T = 1000,
#'   n_ts = 10,
#'   warmup = 50,
#'   totTime = 1050,
#'   mean_rep = 10,
#'   sd_rep = 2,
#'   min_edg_val = 0.05,
#'   max_edg_val = 1,
#'   remain_lower = 0.33,
#'   remain_upper = 0.66
#' )
#' }
#'
#' @export
# Dependencies are loaded centrally via R/dependencies.R
# Required packages: netcontrol, graphicalVAR, abind, mvtnorm, progress

source("R/generation/generate_netdyn.R")
source("R/generation/generate_transmat.R")

generate_timeseries <- function(Density,
                                N,
                                M,
                                T,
                                n_ts,
                                warmup,
                                totTime,
                                mean_rep,
                                sd_rep,
                                min_edg_val,
                                max_edg_val,
                                remain_lower,
                                remain_upper
                                ) {
  
  # Setup progressbar
  print("Generating networks", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )
  
  
  # Generate set of dynamics for each regime
  #set.seed(123)
  DynamicsMatrices_list <- list()
  for (i in seq_along(Density)) {
    dyn_level2 <- list()
    for (j in seq_along(N)) {
      dyn_level3 <- list()
      for (k in 1:n_ts) {
        dyn_level4 <- list()
        for (l in seq_along(M)) {
          dyn_level5 <- list()
          for (m in 1:M[l]) {
            repeat {
              W <- generate_netdyn(N[j], Density[i], min_edg_val, max_edg_val)
              if (all(abs(W[["Wtemp"]][W[["Wtemp"]] != 0]) >= min_edg_val) &&
                  all(abs(W[["Wcont"]][W[["Wcont"]] != 0]) >= min_edg_val) &&
                  any(W[["Wcont"]] != 0)) {
                break
                }
             }
            dyn_level5[[paste0("Regime", m)]] <- W
          }
          pb$tick()
          
          dyn_level4[[paste0(M[l], "_Regimes")]] <- dyn_level5[1:M[l]]
        }
        dyn_level3[[paste0("Set_no.", k)]] <- dyn_level4
      }
      dyn_level2[[paste0(N[j], "_Nodes")]] <- dyn_level3
    }
    DynamicsMatrices_list[[paste0("Density_", Density[i])]] <- dyn_level2
  }
 
  
  # Setup progressbar
  print("Generating timeseries", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(T) * length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )
  
  # Generate timeseries data
  #set.seed(456)
  Timeseries_data <- list()
  for (t in seq_along(totTime)) {
    ts_level1 <- list()
    for (i in seq_along(Density)) {
      ts_level2 <- list()
      for (j in seq_along(N)) {
        ts_level3 <- list()
        for (k in seq_along(M)) {
          ts_level4 <- list()
          all_Wtemp_sd <- numeric()
          for (l in 1:n_ts) {

            # Generate initial vector for timeseries
            init <- runif(N[j], min = 0, max = 5)

            # Get dynamics for current timeseries
            dynamics <- DynamicsMatrices_list[[i]][[j]][[l]][[k]]

            ## Generate TS from Transmat
            # Generate transmat
            transmat <- generate_transmat(M[k], remain_lower, remain_upper)

            # PERFORMANCE: Pre-allocate regime sequence vector
            # This avoids O(n²) vector copying when growing with c()
            Rseq <- integer(totTime[t] - 1)

            # Initialize state matrix
            X <- matrix(init, nrow = totTime[t], ncol = N[j])
            for (m in 2:totTime[t]) {

              # Set starting regime index
              if (m == 2) {
                reg_index <- sample.int(M[k], 1, replace = TRUE)
                Rseq[1] <- reg_index
              }

              # Determine regime in next timestep from current timestep and transmat
              reg_index <- sample.int(M[k], 1, prob = transmat[reg_index, ])
              Rseq[m - 1] <- reg_index

              # Get dynamics for current regime
              curreg_mu <- dynamics[[reg_index]][["mu"]]
              curreg_W_temp <- dynamics[[reg_index]][["W_temp"]]
              curreg_Beta <- dynamics[[reg_index]][["Beta"]]
              curreg_sigma <- dynamics[[reg_index]][["sigma"]]
        
        
              # Generate next state vector from previous x, Beta and residuals from sigma
              X[m, ] <- curreg_mu + 
                t(curreg_Beta %*% (X[m - 1, ] - curreg_mu)) + 
                mvtnorm::rmvnorm(1, rep(0, N[j]), curreg_sigma)
              
            }
            pb$tick()
            
            ts_level4[[paste0("Timeseries_", l)]] <- list("Timeseries" = X[-seq_len(warmup), , drop = FALSE],
                                                          "Regime_sequence" = Rseq,
                                                          "Regime_dynamics" = dynamics,
                                                          "Transmat" = transmat)
          }
          ts_level3[[paste0(M[k], "_Regimes")]] <- ts_level4
        }
        ts_level2[[paste0(N[j], "_Nodes")]] <- ts_level3
      }
      ts_level1[[paste0("Density_", Density[i]*100, "%")]] <- ts_level2
    }
    Timeseries_data[[paste0(T[t], "_Timesteps")]] <- ts_level1
  }
  return(Timeseries_data = Timeseries_data)
}

