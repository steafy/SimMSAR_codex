# Edge-recovery sensitivity/specificity between a true and an estimated matrix,
# treating any non-zero entry as an edge. Returns list(sensitivity, specificity),
# NA when there are no true edges / non-edges.
senspec <- function(orig_mat, est_mat) {
  true_pos <- (orig_mat != 0) & (est_mat != 0)
  true_neg <- (orig_mat == 0) & (est_mat == 0)
  
  n_pos <- sum(orig_mat != 0)
  n_neg <- sum(orig_mat == 0)
  sensitivity <- if (n_pos > 0) sum(true_pos) / n_pos else NA_real_
  specificity <- if (n_neg > 0) sum(true_neg) / n_neg else NA_real_

  return(list("sensitivity" = sensitivity,
              "specificity" = specificity
              )
         ) 
}

