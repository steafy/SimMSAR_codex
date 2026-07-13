#' Generate Network Dynamics for a Single Regime
#'
#' Generates complete network dynamics for one regime, including the temporal
#' (lag-1) and contemporaneous networks, along with associated statistics and
#' controllability measures.
#'
#' @param N Integer. Number of nodes in the network.
#' @param Density Numeric (0 to 1). Target edge density for the network.
#' @param min_edg_val Numeric. Minimum absolute edge weight (edges below this are zero).
#' @param max_edg_val Numeric. Maximum absolute edge weight.
#'
#' @return A list containing network dynamics and statistics:
#'   \describe{
#'     \item{mu}{Mean vector (length N)}
#'     \item{A}{Lag-1 autoregressive coefficient matrix (N × N); the temporal network}
#'     \item{A_sd}{Standard deviation of A values}
#'     \item{A_strength}{Node strength for A network}
#'     \item{A_mean_strength}{Mean node strength for A}
#'     \item{A_density}{Edge density of A}
#'     \item{A_weighted_density}{Weighted density of A}
#'     \item{A_stability}{Stability indicator for A (TRUE/FALSE)}
#'     \item{AC}{Average controllability of A (node-wise vector)}
#'     \item{sigma}{Residual covariance matrix (N × N)}
#'     \item{K}{Precision matrix (N × N); the contemporaneous network. Diagonal
#'       is dominant by construction (see \code{generate_K}) and is kept.}
#'     \item{K_pos.definit}{"Yes"/"No" indicating if K is positive definite}
#'     \item{K_thresh}{Minimum absolute off-diagonal non-zero value in K}
#'     \item{K_density}{Off-diagonal edge density of K}
#'     \item{K_weighted_density}{Off-diagonal weighted density of K}
#'   }
#'
#' @details
#' The function generates network dynamics in several steps:
#'
#' 1. **Mean vector**: Sampled from uniform(0, 5)
#'
#' 2. **Temporal network**: A matrix generated via \code{generate_A};
#'    used directly (raw, unstandardized) as the temporal network.
#'
#' 3. **Contemporaneous network**: K (precision matrix) generated via
#'    \code{generate_K}; Sigma = K^{-1}. K is used directly
#'    (including its dominant diagonal) as the contemporaneous network.
#'
#' 4. **Statistics**: Network density, strength, and average controllability
#'    (for A only) are computed.
#'
#' @note
#' Dependencies are loaded centrally via R/dependencies.R
#' Required packages: Matrix
#'
#' @seealso
#' \code{\link{generate_A}} for temporal network generation
#' \code{\link{generate_K}} for precision matrix generation
#' \code{\link{check_stability}} for stability checking
#' \code{\link{average_controllability}} for the controllability metric
#'
#' @examples
#' \dontrun{
#' # Generate dynamics for 4-node network with 30% density
#' dynamics <- generate_netdyn(N = 4, Density = 0.3,
#'                             min_edg_val = 0.05, max_edg_val = 1)
#' }
#'
#' @export
# Dependencies are loaded centrally via R/dependencies.R
# Required packages: Matrix

# Load functions
source("R/generation/generate_random.R")
source("R/generation/generate_A.R")
source("R/generation/generate_K.R")
source("R/generation/check_stability.R")
source("R/utils/average_controllability.R")


generate_netdyn <- function (N, Density, min_edg_val, max_edg_val){

  # Generate means for each regime
  mu <- runif(N, min = 0, max = 5)

  ## Generate lag-1 regression matrix (A) -- this is the temporal network
  A_list <- generate_A(N, Density, min_edg_val, max_edg_val)
  A <- A_list$A

  # Calculate average controllability for A (temporal network only)
  AC <- average_controllability(A)


  ## Generate contemporaneous network dynamics
  # Make precision matrix K -- this is the contemporaneous network
  K_list <- generate_K(N, Density, min_edg_val, max_edg_val)

  # Calculate sigma from K
  K <- K_list$K
  sigma <- solve(K)

  # Check if all off-diagonal non-zero elements of K are above threshold
  # (K is symmetric, so the lower triangle represents all off-diagonal edges)
  K_offdiag <- K[lower.tri(K, diag = FALSE)]
  K_thresh <- if (any(K_offdiag != 0)) min(abs(K_offdiag[K_offdiag != 0])) else NA_real_

  # Calculate empirical off-diagonal density of K
  K_dens <- mean(K[lower.tri(K, diag = FALSE)] != 0)
  K_weigh_dens <- mean(abs(K[lower.tri(K, diag = FALSE)]))


  return(list(mu = mu,
              A = A_list$A,
              A_sd = A_list$sd_A,
              A_strength = A_list$A_str,
              A_mean_strength = A_list$A_mean_str,
              A_density = A_list$A_dens,
              A_weighted_density = A_list$A_weigh_dens,
              A_stability = A_list$A_stability,
              AC = AC,
              sigma = sigma,
              K = K_list$K,
              K_pos.definit = K_list$K_posdef,
              K_thresh = K_thresh,
              K_density = K_dens,
              K_weighted_density = K_weigh_dens
              ))
  }
