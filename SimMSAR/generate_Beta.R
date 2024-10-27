## Function to generate temporal network dynamics (Beta)

generate_Beta <- function(N, Density, min_edg_val, max_edg_val) {
  
  # Generate Beta with non-zero elements between min_edg_val & max_edg_val
  Beta <- matrix(0, nrow = N, ncol = N)
  num_nonzero <- round(Density * N * N)
  indices <- sample(N*N, num_nonzero, replace = FALSE)
  signs <- sample(c(-1, 1), num_nonzero, replace = TRUE)
  values <- runif(num_nonzero, min_edg_val, max_edg_val) * signs
  Beta[indices] <- values
  
  # Scale Beta to ensure stability
  eigvals <- eigen(Beta)$values
  spectral_radius <- max(abs(eigvals))
  scaling_factor <- ifelse(spectral_radius >= 1, 0.99 / spectral_radius, 1)
  Beta <- Beta * scaling_factor
  
  # Adjust non-zero elements to between min_edg_val & max_edg_val
  adjust_non_zero_elements <- function(Beta, min_edg_val, max_edg_val) {
    non_zero_elements <- which(Beta != 0, arr.ind = TRUE)
    for (idx in 1:nrow(non_zero_elements)) {
      i <- non_zero_elements[idx, 1]
      j <- non_zero_elements[idx, 2]
      Beta[i, j] <- sign(Beta[i, j]) * max(min(abs(Beta[i, j]), max_edg_val), min_edg_val)
    }
    return(Beta)
  }
  
  Beta <- adjust_non_zero_elements(Beta, min_edg_val, max_edg_val)
  
  # Check stability of Beta
  Beta_stability <- check_stability(Beta)
  
  # Calculate sd of Beta
  sd_Beta <- sd(as.vector(Beta))
  
  # Calculate node strength of Beta
  Beta_str <- rowSums(Beta)
  Beta_mean_str <- mean(Beta_str)
 
  # Calculate empirical density of Beta
  Beta_dens <- mean(Beta != 0)
  Beta_weigh_dens <- mean(abs(Beta))
  
 
  return(list(
    Beta = Beta,
    Beta_stability = Beta_stability,
    sd_Beta = sd_Beta,
    Beta_str = Beta_str,
    Beta_mean_str = Beta_mean_str,
    Beta_dens = Beta_dens,
    Beta_weigh_dens = Beta_weigh_dens
  ))
  
}

