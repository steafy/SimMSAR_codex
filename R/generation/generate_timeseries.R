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
#' @return A tibble (timeseries_data object) with the following columns:
#'   \describe{
#'     \item{timesteps}{Number of time steps in the series (excluding warmup)}
#'     \item{density}{Edge density of the network}
#'     \item{nodes}{Number of nodes in the network}
#'     \item{regimes}{Number of regimes in the model}
#'     \item{ts_id}{Time series ID (1 to n_ts)}
#'     \item{timeseries_data}{List-column: Matrix of time series data (T × N)}
#'     \item{regime_sequence}{List-column: Vector of regime indices for each time step}
#'     \item{regime_dynamics}{List-column: List of network dynamics for each regime}
#'     \item{transmat}{List-column: Transition probability matrix between regimes}
#'   }
#'
#'   Each regime_dynamics list contains:
#'   \itemize{
#'     \item mu: Mean vector
#'     \item Wtemp: Temporal network (lag-1 effects)
#'     \item Beta: Autoregressive coefficient matrix
#'     \item sigma: Residual covariance matrix
#'     \item kappa: Precision matrix
#'     \item Wcont: Contemporaneous network (partial correlations)
#'     \item Wtemp_ac: Average controllability for Wtemp
#'     \item Wcont_ac: Average controllability for Wcont
#'   }
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
#' Required packages: netcontrol, graphicalVAR, abind, mvtnorm, progress, tibble, dplyr
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
#'
#' # Access specific time series
#' ts_data %>% filter(timesteps == 1000, density == 0.3, ts_id == 1)
#' }
#'
#' @export
# Dependencies are loaded centrally via R/dependencies.R
# Required packages: netcontrol, graphicalVAR, abind, mvtnorm, progress, tibble, dplyr

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
                                remain_upper,
                                max_attempts = 1000
                                ) {

  # Setup progressbar for network generation
  print("Generating networks", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )


  # Generate set of dynamics for each regime
  # Store in a tibble instead of nested lists
  dynamics_rows <- length(Density) * length(N) * length(M) * n_ts

  dynamics_tibble <- tibble::tibble(
    density = numeric(dynamics_rows),
    nodes = integer(dynamics_rows),
    regimes = integer(dynamics_rows),
    set_id = integer(dynamics_rows),
    dynamics = vector("list", dynamics_rows)
  )

  row_idx <- 0

  for (i in seq_along(Density)) {
    for (j in seq_along(N)) {
      for (k in 1:n_ts) {
        for (l in seq_along(M)) {
          row_idx <- row_idx + 1

          # Generate dynamics for all regimes in this condition
          regime_dynamics <- list()
          for (m in 1:M[l]) {
            attempts <- 0
            repeat {
              attempts <- attempts + 1
              W <- generate_netdyn(N[j], Density[i], min_edg_val, max_edg_val)
              if (all(abs(W[["Wtemp"]][W[["Wtemp"]] != 0]) >= min_edg_val) &&
                  all(abs(W[["Wcont"]][W[["Wcont"]] != 0]) >= min_edg_val) &&
                  any(W[["Wcont"]] != 0)) {
                break
              }
              if (attempts >= max_attempts) {
                stop(
                  "Failed to generate valid network dynamics after ",
                  max_attempts,
                  " attempts. Check Density/min_edg_val/max_edg_val settings."
                )
              }
            }
            regime_dynamics[[paste0("Regime", m)]] <- W
          }

          pb$tick()

          # Store in tibble
          dynamics_tibble$density[row_idx] <- Density[i]
          dynamics_tibble$nodes[row_idx] <- N[j]
          dynamics_tibble$regimes[row_idx] <- M[l]
          dynamics_tibble$set_id[row_idx] <- k
          dynamics_tibble$dynamics[[row_idx]] <- regime_dynamics
        }
      }
    }
  }


  # Setup progressbar for timeseries generation
  print("Generating timeseries", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(T) * length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )

  # Generate timeseries data
  # Pre-allocate tibble
  ts_rows <- length(T) * length(Density) * length(N) * length(M) * n_ts

  timeseries_tibble <- tibble::tibble(
    timesteps = integer(ts_rows),
    density = numeric(ts_rows),
    nodes = integer(ts_rows),
    regimes = integer(ts_rows),
    ts_id = integer(ts_rows),
    timeseries_data = vector("list", ts_rows),
    regime_sequence = vector("list", ts_rows),
    regime_dynamics = vector("list", ts_rows),
    transmat = vector("list", ts_rows)
  )

  ts_idx <- 0

  for (t in seq_along(totTime)) {
    for (i in seq_along(Density)) {
      for (j in seq_along(N)) {
        for (k in seq_along(M)) {
          for (l in 1:n_ts) {
            ts_idx <- ts_idx + 1

            # Generate initial vector for timeseries
            init <- runif(N[j], min = 0, max = 5)

            # Get dynamics for current timeseries from dynamics_tibble
            dynamics <- dynamics_tibble %>%
              filter(density == Density[i], nodes == N[j],
                     regimes == M[k], set_id == l) %>%
              pull(dynamics) %>%
              .[[1]]

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
              curreg_Beta <- dynamics[[reg_index]][["Beta"]]
              curreg_sigma <- dynamics[[reg_index]][["sigma"]]


              # Generate next state vector from previous x, Beta and residuals from sigma
              X[m, ] <- curreg_mu +
                t(curreg_Beta %*% (X[m - 1, ] - curreg_mu)) +
                mvtnorm::rmvnorm(1, rep(0, N[j]), curreg_sigma)

            }
            pb$tick()

            # Store in tibble
            timeseries_tibble$timesteps[ts_idx] <- T[t]
            timeseries_tibble$density[ts_idx] <- Density[i]
            timeseries_tibble$nodes[ts_idx] <- N[j]
            timeseries_tibble$regimes[ts_idx] <- M[k]
            timeseries_tibble$ts_id[ts_idx] <- l
            timeseries_tibble$timeseries_data[[ts_idx]] <- X[-seq_len(warmup), , drop = FALSE]
            timeseries_tibble$regime_sequence[[ts_idx]] <- Rseq
            timeseries_tibble$regime_dynamics[[ts_idx]] <- dynamics
            timeseries_tibble$transmat[[ts_idx]] <- transmat
          }
        }
      }
    }
  }

  # Add S3 class
  class(timeseries_tibble) <- c("timeseries_data", class(timeseries_tibble))

  return(timeseries_tibble)
}


#' Print Method for timeseries_data
#'
#' @param x A timeseries_data object
#' @param ... Additional arguments (unused)
#'
#' @export
print.timeseries_data <- function(x, ...) {
  cat("MSAR Timeseries Data\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Total time series: %d\n", nrow(x)))
  cat(sprintf("Conditions:\n"))
  cat(sprintf("  Timesteps: %s\n", paste(unique(x$timesteps), collapse = ", ")))
  cat(sprintf("  Density: %s\n", paste(unique(x$density), collapse = ", ")))
  cat(sprintf("  Nodes: %s\n", paste(unique(x$nodes), collapse = ", ")))
  cat(sprintf("  Regimes: %s\n", paste(unique(x$regimes), collapse = ", ")))
  cat(sprintf("  Time series per condition: %d\n", max(x$ts_id)))
  cat("════════════════════════════════════════════════════════════════\n")
  cat("\nUse dplyr::filter() to subset by conditions.\n")
  cat("Access data: $timeseries_data, $regime_dynamics, etc.\n")
  NextMethod()
}
