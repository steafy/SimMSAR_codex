# Build an M x M Markov transition matrix: each row's self-transition probability
# is drawn from U(remain_lower, remain_upper), the remainder split randomly over
# the off-diagonal entries so each row sums to 1.
generate_transmat <- function(M, remain_lower, remain_upper) {
  # Initialize transmat 
  transmat <- matrix(0, nrow = M, ncol = M)
  
  for (i in 1:M) {
    # Set probabilities to remain in current regime
    remain_val <- runif(1, remain_lower, remain_upper)
    transmat[i, i] <- remain_val
    
    # Set probabilities for regime switching in next timestep
    switch_val <- 1 - remain_val
    if (M > 1) {
      rest <- runif(M - 1)
      rest <- rest / sum(rest) * switch_val
      transmat[i, -i] <- rest
    }
  }
  
  return(transmat)
}

