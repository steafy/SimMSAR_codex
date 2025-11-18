#' Assign Estimated Regimes to True Regimes
#'
#' Maps estimated regimes to original (true) regimes based on maximum correlation
#' matching. Solves the regime label-switching problem by finding the optimal
#' one-to-one assignment that maximizes total correlation.
#'
#' @param cor_results Matrix of correlations between original and estimated regimes.
#'   Rows represent original regimes, columns represent estimated regimes.
#'   Typically correlations of vectorized network matrices (Wtemp or Wcont).
#'
#' @return Matrix with n rows (number of regimes) and 3 columns:
#'   \describe{
#'     \item{orig_Reg_No}{Original regime index (1 to n)}
#'     \item{est_Reg_No}{Matched estimated regime index}
#'     \item{Value}{Correlation value for this pairing}
#'   }
#'
#' @details
#' The function solves the assignment problem using a greedy algorithm:
#'
#' 1. **Sort**: Orders rows by their maximum correlation (descending)
#'
#' 2. **Sequential Assignment**: For each row (starting with highest max correlation):
#'    \itemize{
#'      \item Finds the column with highest correlation among unused columns
#'      \item Assigns that pairing
#'      \item Marks the column as used
#'    }
#'
#' This greedy approach ensures:
#' \itemize{
#'   \item Each estimated regime is assigned to at most one true regime
#'   \item Regimes with clearer matches (higher correlations) are assigned first
#'   \item The label-switching problem in mixture models is resolved
#' }
#'
#' @note
#' This is critical for MSAR models because regime labels are arbitrary.
#' The EM algorithm may converge to a solution where estimated regime 1
#' corresponds to true regime 2, etc. This function corrects that.
#'
#' @examples
#' \dontrun{
#' # Correlation matrix: rows = true regimes, cols = estimated regimes
#' cor_matrix <- matrix(c(0.85, 0.35, 0.40, 0.80), 2, 2)
#' assignments <- asign_regimes(cor_matrix)
#' # Result: true regime 1 -> est regime 1 (0.85)
#' #         true regime 2 -> est regime 2 (0.80)
#' }
#'
#' @export
asign_regimes <- function(cor_results) {
  
  # Add column with indices for the original regimes to correlation results
  cor_results <- cbind(cor_results, "orig. Regime No." = 1:nrow(cor_results))
  
  # Determine highest correlation for each row
  max_vals <- numeric(nrow(cor_results))
  for (i in 1:nrow(cor_results)) {
    max_vals[i] <- max(cor_results[i, -ncol(cor_results)])
    }
  
  # Sort matrix with respect to highest correlation in descending order
  cor_matrix_sorted <- cor_results[order(max_vals, decreasing = TRUE), ]
  
  # Select highest values, without using a column more than once
  n <- nrow(cor_matrix_sorted)
  asign_regs <- matrix(nrow = n, ncol = 3, dimnames = list(NULL, c("orig_Reg_No", "est_Reg_No", "Value")))
  used_cols <- integer(0)
  
  # Get correlation values of each row excluding columns with previous higher values 
  for (i in 1:n) {
    cur_row_vals <- cor_matrix_sorted[i, -ncol(cor_matrix_sorted)]
    avail_cols <- setdiff(1:length(cur_row_vals), used_cols)
    
    if (length(avail_cols) > 0) {
      max_col <- which.max(cur_row_vals[avail_cols])
      asign_regs[i, ] <- c(cor_matrix_sorted[i, ncol(cor_matrix_sorted)], avail_cols[max_col], cur_row_vals[avail_cols[max_col]])
      used_cols <- c(used_cols, avail_cols[max_col])
    } else {
      asign_regs[i, ] <- c(cor_matrix_sorted[i, ncol(cor_matrix_sorted)], NA, NA)
    }
  }
  return(asign_regs)
}


