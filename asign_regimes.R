# Function to map estimated regimes on original regimes based on correlations

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


