#' Check Stability of a Matrix
#'
#' Determines if a matrix is stable by checking whether all eigenvalues have
#' absolute value less than 1. Used to verify stability of autoregressive
#' coefficient matrices in VAR models.
#'
#' @param Matrix Numeric matrix to check for stability.
#'
#' @return Logical:
#'   \describe{
#'     \item{TRUE}{All eigenvalues have |lambda| < 1 (stable)}
#'     \item{FALSE}{At least one eigenvalue has |lambda| >= 1 (unstable)}
#'   }
#'
#' @details
#' A matrix is considered stable if its spectral radius (maximum absolute eigenvalue)
#' is less than 1. For autoregressive models, this ensures:
#' \itemize{
#'   \item The process is stationary
#'   \item Impulse responses decay over time
#'   \item The variance is finite
#' }
#'
#' The function computes all eigenvalues and checks if max(|eigenvalues|) < 1.
#'
#' @note
#' This stability criterion applies to VAR(1) models. For higher-order VAR models,
#' the companion matrix should be checked instead.
#'
#' @seealso
#' \code{\link{generate_Beta}} which uses this function
#' \code{\link{generate_netdyn}} which reports stability
#'
#' @examples
#' \dontrun{
#' # Check stability of a simple matrix
#' M <- matrix(c(0.3, 0.1, 0.2, 0.4), 2, 2)
#' check_stability(M)  # Should return "Yes"
#'
#' # Unstable matrix
#' M2 <- matrix(c(0.8, 0.5, 0.6, 0.7), 2, 2)
#' check_stability(M2)  # Likely returns "No"
#' }
#'
#' @export
check_stability <- function(Matrix) {
  eigen <- eigen(Matrix)
  stability <- all(abs(eigen$values) < 1)
  
  return(stability)
}
