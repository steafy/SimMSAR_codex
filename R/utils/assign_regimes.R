# Optimally assign estimated regimes to true regimes from a true-vs-estimated
# correlation matrix, via the Hungarian algorithm (clue::solve_LSAP on
# max(cor) - cor). Returns a matrix with columns orig_Reg_No, est_Reg_No, Value.
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


