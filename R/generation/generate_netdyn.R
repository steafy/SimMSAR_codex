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
#'     \item{Beta}{Lag-1 autoregressive coefficient matrix (N × N); the temporal network}
#'     \item{Beta_sd}{Standard deviation of Beta values}
#'     \item{Beta_strength}{Node strength for Beta network}
#'     \item{Beta_mean_strength}{Mean node strength for Beta}
#'     \item{Beta_density}{Edge density of Beta}
#'     \item{Beta_weighted_density}{Weighted density of Beta}
#'     \item{Beta_stability}{Stability indicator for Beta (TRUE/FALSE)}
#'     \item{Beta_ac}{Average controllability of Beta (node-wise vector)}
#'     \item{sigma}{Residual covariance matrix (N × N)}
#'     \item{kappa}{Precision matrix (N × N); the contemporaneous network. Diagonal
#'       is dominant by construction (see \code{generate_kappa}) and is kept.}
#'     \item{kappa_pos.definit}{"Yes"/"No" indicating if kappa is positive definite}
#'     \item{kappa_thresh}{Minimum absolute off-diagonal non-zero value in kappa}
#'     \item{kappa_density}{Off-diagonal edge density of kappa}
#'     \item{kappa_weighted_density}{Off-diagonal weighted density of kappa}
#'   }
#'
#' @details
#' The function generates network dynamics in several steps:
#'
#' 1. **Mean vector**: Sampled from uniform(0, 5)
#'
#' 2. **Temporal network**: Beta matrix generated via \code{generate_Beta};
#'    used directly (raw, unstandardized) as the temporal network.
#'
#' 3. **Contemporaneous network**: Kappa (precision matrix) generated via
#'    \code{generate_kappa}; Sigma = kappa^{-1}. Kappa is used directly
#'    (including its dominant diagonal) as the contemporaneous network.
#'
#' 4. **Statistics**: Network density, strength, and average controllability
#'    (for Beta only) are computed.
#'
#' @note
#' Dependencies are loaded centrally via R/dependencies.R
#' Required packages: Matrix
#'
#' @seealso
#' \code{\link{generate_Beta}} for temporal network generation
#' \code{\link{generate_kappa}} for precision matrix generation
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
source("R/generation/generate_Beta.R")
source("R/generation/generate_kappa.R")
source("R/generation/check_stability.R")
source("R/utils/average_controllability.R")


generate_netdyn <- function (N, Density, min_edg_val, max_edg_val){

  # Generate means for each regime
  mu <- runif(N, min = 0, max = 5)

  ## Generate lag-1 regression matrix (Beta) -- this is the temporal network
  Beta_list <- generate_Beta(N, Density, min_edg_val, max_edg_val)
  Beta <- Beta_list$Beta

  # Calculate average controllability for Beta (temporal network only)
  Beta_ac <- average_controllability(Beta)


  ## Generate contemporaneous network dynamics
  # Make precision matrix kappa -- this is the contemporaneous network
  kappa_list <- generate_kappa(N, Density, min_edg_val, max_edg_val)

  # Calculate sigma from kappa
  kappa <- kappa_list$kappa
  sigma <- solve(kappa)

  # Check if all off-diagonal non-zero elements of kappa are above threshold
  # (kappa is symmetric, so the lower triangle represents all off-diagonal edges)
  kappa_offdiag <- kappa[lower.tri(kappa, diag = FALSE)]
  kappa_thresh <- if (any(kappa_offdiag != 0)) min(abs(kappa_offdiag[kappa_offdiag != 0])) else NA_real_

  # Calculate empirical off-diagonal density of kappa
  kappa_dens <- mean(kappa[lower.tri(kappa, diag = FALSE)] != 0)
  kappa_weigh_dens <- mean(abs(kappa[lower.tri(kappa, diag = FALSE)]))


  return(list(mu = mu,
              Beta = Beta_list$Beta,
              Beta_sd = Beta_list$sd_Beta,
              Beta_strength = Beta_list$Beta_str,
              Beta_mean_strength = Beta_list$Beta_mean_str,
              Beta_density = Beta_list$Beta_dens,
              Beta_weighted_density = Beta_list$Beta_weigh_dens,
              Beta_stability = Beta_list$Beta_stability,
              Beta_ac = Beta_ac,
              sigma = sigma,
              kappa = kappa_list$kappa,
              kappa_pos.definit = kappa_list$kappa_posdef,
              kappa_thresh = kappa_thresh,
              kappa_density = kappa_dens,
              kappa_weighted_density = kappa_weigh_dens
              ))
  }
