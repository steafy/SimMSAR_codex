# Average controllability of temporal network A: per-node sum, over a T_ac-step
# horizon, of the squared reachability contributions colSums((A^t)^2). Returns a
# length-N vector.
average_controllability <- function(A, T_ac = 25) {
  n <- nrow(A)
  ac <- numeric(n)
  P <- diag(n)

  for (t in 0:(T_ac - 1)) {
    if (t > 0) {
      P <- A %*% P
    }
    ac <- ac + colSums(P^2)
  }

  return(ac)
}
