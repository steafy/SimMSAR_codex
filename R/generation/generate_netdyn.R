#' Generate Network Dynamics for a Single Regime
#'
#' Generates complete network dynamics for one regime, including temporal (lag-1) and
#' contemporaneous networks, along with associated statistics and controllability measures.
#'
#' @param N Integer. Number of nodes in the network.
#' @param Density Numeric (0 to 1). Target edge density for the network.
#' @param min_edg_val Numeric. Minimum absolute edge weight (edges below this are zero).
#' @param max_edg_val Numeric. Maximum absolute edge weight.
#'
#' @return A list containing network dynamics and statistics:
#'   \describe{
#'     \item{mu}{Mean vector (length N)}
#'     \item{Beta}{Lag-1 autoregressive coefficient matrix (N × N)}
#'     \item{Beta_sd}{Standard deviation of Beta values}
#'     \item{Beta_strength}{Node strength for Beta network}
#'     \item{Beta_mean_strength}{Mean node strength for Beta}
#'     \item{Beta_density}{Edge density of Beta}
#'     \item{Beta_weighted_density}{Weighted density of Beta}
#'     \item{Beta_stability}{Stability indicator for Beta (TRUE/FALSE)}
#'     \item{Wtemp}{Temporal partial correlation matrix (N × N)}
#'     \item{Wtemp_thresh}{Minimum non-zero value in Wtemp}
#'     \item{Wtemp_density}{Edge density of Wtemp}
#'     \item{Wtemp_weighted_density}{Weighted density of Wtemp}
#'     \item{Wtemp_ac}{Average controllability for Wtemp}
#'     \item{sigma}{Residual covariance matrix (N × N)}
#'     \item{kappa}{Precision matrix (N × N)}
#'     \item{Wcont}{Contemporaneous partial correlation matrix (N × N)}
#'     \item{Wcont_thresh}{Minimum non-zero value in Wcont}
#'     \item{Wcont_strength}{Node strength for Wcont}
#'     \item{Wcont_mean_strength}{Mean node strength for Wcont}
#'     \item{Wcont_density}{Edge density of Wcont}
#'     \item{Wcont_weighted_density}{Weighted density of Wcont}
#'     \item{kappa_pos.definit}{"Yes"/"No" indicating if kappa is positive definite}
#'     \item{Wcont_stability}{Stability indicator for Wcont (TRUE/FALSE)}
#'     \item{Wcont_ac}{Average controllability for Wcont}
#'   }
#'
#' @details
#' The function generates network dynamics in several steps:
#'
#' 1. **Mean vector**: Sampled from uniform(0, 5)
#'
#' 2. **Temporal network**:
#'    \itemize{
#'      \item Beta matrix generated via \code{generate_Beta}
#'      \item Wtemp calculated as standardized Beta: Wtemp[i,j] = Beta[i,j] / sqrt(sigma[i,i] * kappa[j,j] + Beta[i,j]^2)
#'    }
#'
#' 3. **Contemporaneous network**:
#'    \itemize{
#'      \item Kappa (precision matrix) generated via \code{generate_kappa}
#'      \item Sigma = kappa^{-1}
#'      \item Wcont calculated from kappa: Wcont[i,j] = -kappa[i,j] / sqrt(kappa[i,i] * kappa[j,j])
#'    }
#'
#' 4. **Statistics**: Network density, strength, controllability computed for both networks
#'
#' @note
#' Dependencies are loaded centrally via R/dependencies.R
#' Required packages: netcontrol, Matrix
#'
#' @seealso
#' \code{\link{generate_Beta}} for temporal network generation
#' \code{\link{generate_kappa}} for precision matrix generation
#' \code{\link{check_stability}} for stability checking
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
# Required packages: netcontrol, Matrix

# Load functions
source("R/generation/generate_random.R")
source("R/generation/generate_Beta.R")
source("R/generation/generate_kappa.R")
source("R/generation/check_stability.R")


generate_netdyn <- function (N, Density, min_edg_val, max_edg_val){

  # Generate means for each regime
  mu <- runif(N, min = 0, max = 5)
 
  ## Generate lag-1 regression matrix (Beta)
  Beta_list <- generate_Beta(N, Density, min_edg_val, max_edg_val)
  Beta <- Beta_list$Beta
  
  
  ## Generate contemporaneous network dynamics
  # Make precision matrix kappa
  kappa_list <- generate_kappa(N, Density, min_edg_val, max_edg_val)

  # Calculate sigma from kappa
  kappa <- kappa_list$kappa
  sigma <- solve(kappa)

  # Calculate matrix of partial correlations (Wcont) from kappa
  Wcont <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        Wcont[i, j] <- -kappa[i, j] / sqrt(kappa[i, i] * kappa[j, j])
      }
    }
  }

  # Check if all correlations in Wconr are above 0.1
  Wcont_thresh <- min(abs(Wcont[Wcont != 0]))
  
  # Check stablity of Wcont
  Wcont_stability <- check_stability(Wcont)

  # Calculate node strength
  Wcont_str <- rowSums(Wcont)
  Wcont_mean_str <- mean(Wcont_str)

  # Calculates empirical density of Wcont
  # Wcont_dens <- sum(Wcont[lower.tri(Wcont, diag = FALSE)] != 0) / ((N*(N-1))/2)
  Wcont_dens <- mean(Wcont[lower.tri(Wcont, diag = FALSE)] != 0)
  Wcont_weigh_dens <- mean(abs(Wcont[lower.tri(Wcont, diag = FALSE)]))


  # Calculate average controllability for Wcont
  Wcont_ac <- ave_control_centrality(Wcont)

  # Calculate Wtemp (parial correlation matrix) from Beta, kappa and sigma
  Wtemp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {   # if raus if (i != j) {
      #if (i != j) {
        Wtemp[i, j] <- Beta[i, j] / sqrt(sigma[i, i] %*% kappa[j ,j] + Beta[i, j]^2)
     # }
    }
  }

  # Calculate empirical density of Wtemp
  Wtemp_dens <- mean(Wtemp != 0)
  Wtemp_weigh_dens <- mean(abs(Wtemp))

  # Check if all correlations in Wtemp are above 0.1
  Wtemp_thresh <- min(abs(Wtemp[Wtemp != 0]))
  
  # Calculate average controllability for Wtemp
  Wtemp_ac <- ave_control_centrality(Wtemp)


  return(list(mu = mu,
              Beta = Beta_list$Beta,
              Beta_sd = Beta_list$sd_Beta,
              Beta_strength = Beta_list$Beta_str,
              Beta_mean_strength = Beta_list$Beta_mean_str,
              Beta_density = Beta_list$Beta_dens,
              Beta_weighted_density = Beta_list$Beta_weigh_dens,
              Beta_stability = Beta_list$Beta_stability,
              Wtemp = Wtemp,
              Wtemp_thresh = Wtemp_thresh,
              Wtemp_density = Wtemp_dens,
              Wtemp_weighted_density = Wtemp_weigh_dens,
              Wtemp_ac = Wtemp_ac,
              sigma = sigma,
              kappa = kappa_list$kappa,
              Wcont = Wcont,
              Wcont_thresh = Wcont_thresh,
              Wcont_strength = Wcont_str,
              Wcont_mean_strength = Wcont_mean_str,
              Wcont_density = Wcont_dens,
              Wcont_weighted_density = Wcont_weigh_dens,
              kappa_pos.definit = kappa_list$kappa_posdef,
              Wcont_stability = Wcont_stability,
              Wcont_ac = Wcont_ac
              ))
  }
