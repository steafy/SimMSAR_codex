#' Generate Temporal Network Dynamics (Beta Matrix)
#'
#' Generates a stable autoregressive coefficient matrix (Beta) for VAR(1) models with
#' specified edge density and weight constraints. Ensures stability by scaling to keep
#' spectral radius below 1.
#'
#' @param N Integer. Number of nodes in the network.
#' @param Density Numeric (0 to 1). Target proportion of non-zero edges in Beta.
#' @param min_edg_val Numeric. Minimum absolute value for non-zero edges.
#' @param max_edg_val Numeric. Maximum absolute value for non-zero edges.
#'
#' @return A list containing:
#'   \describe{
#'     \item{Beta}{Autoregressive coefficient matrix (N × N)}
#'     \item{Beta_stability}{Logical. TRUE if stable (spectral radius < 1)}
#'     \item{sd_Beta}{Standard deviation of all Beta values}
#'     \item{Beta_str}{Node strength vector (row sums)}
#'     \item{Beta_mean_str}{Mean node strength}
#'     \item{Beta_dens}{Empirical edge density (proportion of non-zero elements)}
#'     \item{Beta_weigh_dens}{Weighted density (mean absolute value)}
#'   }
#'
#' @details
#' The function follows these steps:
#'
#' 1. **Initial generation**: Creates sparse Beta matrix with random edge weights
#'    sampled uniformly from [min_edg_val, max_edg_val] with random signs
#'
#' 2. **Stability enforcement**: Scales Beta to ensure spectral radius < 1:
#'    \itemize{
#'      \item Computes eigenvalues and spectral radius
#'      \item If spectral radius >= 1, scales by 0.99/spectral_radius
#'    }
#'
#' 3. **Edge adjustment**: Re-adjusts non-zero edges to be within [min_edg_val, max_edg_val]
#'    while preserving signs and stability
#'
#' 4. **Statistics**: Computes density, strength, and stability metrics
#'
#' @note
#' The stability constraint ensures the VAR process is stationary. This may result in
#' actual edge densities slightly different from the target Density parameter.
#'
#' @seealso
#' \code{\link{check_stability}} for stability verification
#' \code{\link{generate_netdyn}} which uses this function
#'
#' @examples
#' \dontrun{
#' # Generate stable Beta matrix for 5 nodes with 30% density
#' beta_result <- generate_Beta(N = 5, Density = 0.3,
#'                              min_edg_val = 0.05, max_edg_val = 1)
#' print(beta_result$Beta_stability)  # Should be TRUE
#' }
#'
#' @export
generate_Beta <- function(N, Density, min_edg_val, max_edg_val) {
  
  # Generate Beta with non-zero elements between min_edg_val & max_edg_val
  Beta <- matrix(0, nrow = N, ncol = N)
  num_nonzero <- round(Density * N * N)
  indices <- sample(N*N, num_nonzero, replace = FALSE)
  signs <- sample(c(-1, 1), num_nonzero, replace = TRUE)
  values <- runif(num_nonzero, min_edg_val, max_edg_val) * signs
  Beta[indices] <- values
  
  # Scale Beta to ensure stability
  eigvals <- eigen(Beta)$values
  spectral_radius <- max(abs(eigvals))
  scaling_factor <- ifelse(spectral_radius >= 1, 0.99 / spectral_radius, 1)
  Beta <- Beta * scaling_factor
  
  # Adjust non-zero elements to between min_edg_val & max_edg_val
  adjust_non_zero_elements <- function(Beta, min_edg_val, max_edg_val) {
    non_zero_elements <- which(Beta != 0, arr.ind = TRUE)
    for (idx in 1:nrow(non_zero_elements)) {
      i <- non_zero_elements[idx, 1]
      j <- non_zero_elements[idx, 2]
      Beta[i, j] <- sign(Beta[i, j]) * max(min(abs(Beta[i, j]), max_edg_val), min_edg_val)
    }
    return(Beta)
  }
  
  Beta <- adjust_non_zero_elements(Beta, min_edg_val, max_edg_val)
  
  # Check stability of Beta
  Beta_stability <- check_stability(Beta)
  
  # Calculate sd of Beta
  sd_Beta <- sd(as.vector(Beta))
  
  # Calculate node strength of Beta
  Beta_str <- rowSums(Beta)
  Beta_mean_str <- mean(Beta_str)
 
  # Calculate empirical density of Beta
  Beta_dens <- mean(Beta != 0)
  Beta_weigh_dens <- mean(abs(Beta))
  
 
  return(list(
    Beta = Beta,
    Beta_stability = Beta_stability,
    sd_Beta = sd_Beta,
    Beta_str = Beta_str,
    Beta_mean_str = Beta_mean_str,
    Beta_dens = Beta_dens,
    Beta_weigh_dens = Beta_weigh_dens
  ))
  
}

