# Check stability of Matrix
check_stability <- function(Matrix) {
  eigen <- eigen(Matrix)
  stability <- ifelse(all(abs(eigen$values) < 1), "Yes", "No")
  
  return(stability)
}