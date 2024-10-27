# Function to generate transmat

generate_transmat <- function(M) {
  # Initialize transmat 
  transmat <- matrix(0, nrow = M, ncol = M)
  
  for (i in 1:M) {
    # Set probabilities to remain in current regime
    remain_val <- runif(1, 0.92, 0.95)
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

