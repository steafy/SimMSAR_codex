#' Calculate Mean Absolute Error Between Two Matrices
#'
#' Computes the mean absolute error (MAE) between two matrices of the same dimensions.
#' Measures average magnitude of errors in network edge weight estimation.
#'
#' @param A Matrix. True/original network matrix.
#' @param B Matrix. Estimated network matrix (same dimensions as A).
#'
#' @return Numeric. Mean absolute error averaged across all matrix elements.
#'
#' @details
#' The MAE is calculated as:
#'
#' \deqn{MAE = \frac{1}{n \times m} \sum_{i=1}^{n} \sum_{j=1}^{m} |A_{ij} - B_{ij}|}
#'
#' where n × m is the total number of matrix elements.
#'
#' **Interpretation**:
#' \itemize{
#'   \item MAE = 0: Perfect match
#'   \item Smaller MAE: Better estimation accuracy
#'   \item MAE measures both edge presence errors AND edge weight errors
#'   \item Scale depends on edge weight range (e.g., [-1, 1] vs [0, 5])
#' }
#'
#' Unlike correlation (which measures pattern similarity), MAE measures
#' absolute accuracy of edge weight recovery.
#'
#' @note
#' \itemize{
#'   \item Assumes A and B have identical dimensions
#'   \item All elements contribute equally (no weighting by importance)
#'   \item For sparse networks, many elements are zero, which reduces MAE
#'   \item Sensitive to outliers in edge weight differences
#' }
#'
#' @seealso
#' \code{\link{estimate_MSAR}} which uses this to evaluate recovery accuracy
#' \code{\link{senspec}} for edge detection accuracy (binary)
#'
#' @examples
#' \dontrun{
#' # True network
#' true_net <- matrix(c(0, 0.5, 0.3, 0.7), 2, 2)
#'
#' # Estimated network
#' est_net <- matrix(c(0, 0.6, 0.2, 0.8), 2, 2)
#'
#' mae <- calculate_MAE(true_net, est_net)
#' # MAE = (|0-0| + |0.5-0.6| + |0.3-0.2| + |0.7-0.8|) / 4
#' #     = (0 + 0.1 + 0.1 + 0.1) / 4 = 0.075
#' }
#'
#' @export
calculate_MAE <- function(A, B) {
 
  # Number of elements from matrix A
  n_elements <- prod(dim(A))
  
  # Absolute error for all nodes between matrices
  abs_errors <- abs(A - B)
  
  # Calculate MAE
  MAE <- sum(abs_errors) / n_elements
  
  return(MAE)
}