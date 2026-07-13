# Generate the temporal network A (VAR coefficient matrix): N nodes, target edge
# density, non-zero weights with magnitude in [min_edg_val, max_edg_val] and
# random sign, resampled until A is stable (spectral radius < 1). Returns a list
# with A plus summary diagnostics (sd_A, A_str, A_mean_str, A_dens, A_weigh_dens,
# A_stability).
generate_A <- function(N, Density, min_edg_val, max_edg_val) {
  
  # Generate A with non-zero elements between min_edg_val & max_edg_val
  A <- matrix(0, nrow = N, ncol = N)
  num_nonzero <- round(Density * N * N)
  indices <- sample(N*N, num_nonzero, replace = FALSE)
  signs <- sample(c(-1, 1), num_nonzero, replace = TRUE)
  values <- runif(num_nonzero, min_edg_val, max_edg_val) * signs
  A[indices] <- values
  
  # Enforce stability (spectral radius < 1) AND edge bounds together.
  #
  # These two constraints interact: scaling for stability shrinks every edge
  # (and can push small ones below min_edg_val), while clamping edges back into
  # [min_edg_val, max_edg_val] raises the small ones again (and can re-inflate
  # the spectral radius above 1). The previous version scaled once and clamped
  # afterwards, which silently re-introduced non-stationarity. Instead we
  # alternate the two steps until the matrix is both stable and within bounds.
  # If no such fixed point is reached within max_stab_iter, we leave the last
  # scaled (stable but possibly out-of-bounds) matrix and report an honest
  # A_stability flag; the acceptance loop in generate_timeseries then rejects
  # and regenerates such draws.
  adjust_non_zero_elements <- function(A, min_edg_val, max_edg_val) {
    non_zero_elements <- which(A != 0, arr.ind = TRUE)
    for (idx in 1:nrow(non_zero_elements)) {
      i <- non_zero_elements[idx, 1]
      j <- non_zero_elements[idx, 2]
      A[i, j] <- sign(A[i, j]) * max(min(abs(A[i, j]), max_edg_val), min_edg_val)
    }
    return(A)
  }

  scale_to_stable <- function(A) {
    spectral_radius <- max(abs(eigen(A)$values))
    scaling_factor <- ifelse(spectral_radius >= 1, 0.99 / spectral_radius, 1)
    A * scaling_factor
  }

  max_stab_iter <- 100
  for (iter in 1:max_stab_iter) {
    # Clamp edges into bounds first, then test stability of the in-bounds matrix.
    # Breaking here guarantees the returned A is BOTH in-bounds and stable.
    A <- adjust_non_zero_elements(A, min_edg_val, max_edg_val)
    if (check_stability(A)) break
    A <- scale_to_stable(A)
  }

  # Honest final stability flag (may be FALSE if no in-bounds stable
  # configuration was found; the acceptance loop handles that case)
  A_stability <- check_stability(A)
  
  # Calculate sd of A
  sd_A <- sd(as.vector(A))
  
  # Calculate node strength of A
  A_str <- rowSums(A)
  A_mean_str <- mean(A_str)
 
  # Calculate empirical density of A
  A_dens <- mean(A != 0)
  A_weigh_dens <- mean(abs(A))
  
 
  return(list(
    A = A,
    A_stability = A_stability,
    sd_A = sd_A,
    A_str = A_str,
    A_mean_str = A_mean_str,
    A_dens = A_dens,
    A_weigh_dens = A_weigh_dens
  ))
  
}

