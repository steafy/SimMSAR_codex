# TRUE if a VAR coefficient matrix is stable (all eigenvalue moduli < 1).
check_stability <- function(Matrix) {
  eigen <- eigen(Matrix)
  stability <- all(abs(eigen$values) < 1)
  
  return(stability)
}
