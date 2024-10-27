# Function to calculate true positive rate (sensitivity) and true negative rate (specificity)
# of original and estimated networks

senspec <- function(orig_mat, est_mat) {
  true_pos <- (orig_mat != 0) & (est_mat != 0)
  true_neg <- (orig_mat == 0) & (est_mat == 0)
  
  sensitivity <- sum(true_pos) / sum(orig_mat != 0)
  specificity <- sum(true_neg) / sum(orig_mat == 0)

  return(list("sensitivity" = sensitivity,
              "specificity" = specificity
              )
         ) 
}

