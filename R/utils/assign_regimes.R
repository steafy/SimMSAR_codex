#' Assign Estimated Regimes to True Regimes
#'
#' Maps estimated regimes to original (true) regimes based on maximum correlation
#' matching. Solves the regime label-switching problem by finding the globally
#' optimal one-to-one assignment that maximizes total correlation.
#'
#' @param cor_results Matrix of correlations between original and estimated regimes.
#'   Rows represent original regimes, columns represent estimated regimes.
#'   Typically correlations of vectorized network matrices (A or K).
#'
#' @return Matrix with n rows (number of regimes) and 3 columns:
#'   \describe{
#'     \item{orig_Reg_No}{Original regime index (1 to n)}
#'     \item{est_Reg_No}{Matched estimated regime index}
#'     \item{Value}{Correlation value for this pairing}
#'   }
#'
#' @details
#' The function solves the assignment problem using the Hungarian algorithm
#' (Kuhn-Munkres algorithm) via the \code{clue} package. This guarantees the
#' globally optimal one-to-one matching that maximizes total correlation.
#'
#' The Hungarian algorithm is superior to greedy approaches because it considers
#' all possible assignments simultaneously, avoiding suboptimal local decisions.
#'
#' For example, with correlation matrix:
#' \preformatted{
#'         Est1   Est2
#' True1   0.90   0.85
#' True2   0.80   0.10
#' }
#'
#' A greedy algorithm would assign True1→Est1 (0.90), forcing True2→Est2 (0.10),
#' with total 1.00. The Hungarian algorithm finds the optimal: True1→Est2 (0.85),
#' True2→Est1 (0.80), with total 1.65.
#'
#' @note
#' This is critical for MSAR models because regime labels are arbitrary.
#' The EM algorithm may converge to a solution where estimated regime 1
#' corresponds to true regime 2, etc. This function corrects that optimally.
#'
#' Requires the \code{clue} package for \code{solve_LSAP}.
#'
#' @examples
#' \dontrun{
#' # Correlation matrix: rows = true regimes, cols = estimated regimes
#' cor_matrix <- matrix(c(0.90, 0.80, 0.85, 0.10), 2, 2)
#' assignments <- asign_regimes(cor_matrix)
#' # Optimal result: true regime 1 -> est regime 2 (0.85)
#' #                 true regime 2 -> est regime 1 (0.80)
#' # Total correlation: 1.65 (vs 1.00 with greedy)
#' }
#'
#' @seealso \code{\link[clue]{solve_LSAP}} for the Hungarian algorithm implementation
#'
#' @export
assign_regimes <- function(cor_results) {

  n <- nrow(cor_results)

  # Handle single regime case

if (n == 1) {
    assign_regs <- matrix(c(1, 1, cor_results[1, 1]), nrow = 1, ncol = 3,
                         dimnames = list(NULL, c("orig_Reg_No", "est_Reg_No", "Value")))
    return(assign_regs)
  }

  # Use Hungarian algorithm to find optimal assignment

  # solve_LSAP minimizes cost, so we convert correlations to costs
  # by using (1 - correlation) or max - correlation
  cost_matrix <- max(cor_results) - cor_results

  # Solve the Linear Sum Assignment Problem
  assignment <- clue::solve_LSAP(cost_matrix, maximum = FALSE)

  # Build result matrix
  assign_regs <- matrix(nrow = n, ncol = 3,
                       dimnames = list(NULL, c("orig_Reg_No", "est_Reg_No", "Value")))

  for (i in 1:n) {
    est_idx <- as.integer(assignment[i])
    assign_regs[i, ] <- c(i, est_idx, cor_results[i, est_idx])
  }

  return(assign_regs)
}


