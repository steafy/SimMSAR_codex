# Dependencies are loaded centrally via R/dependencies.R
# Required packages: Matrix

# Load functions
source("R/generation/generate_random.R")
source("R/generation/generate_A.R")
source("R/generation/generate_K.R")
source("R/generation/check_stability.R")
source("R/utils/average_controllability.R")


# Assemble one regime's network dynamics: a stable temporal network A
# (generate_A), a contemporaneous precision matrix K with Sigma = K^{-1}
# (generate_K), the residual covariance Sigma, and A's average controllability.
# Returns a named list (mu, A, AC, sigma, K, plus per-network diagnostics).
generate_netdyn <- function (N, Density, min_edg_val, max_edg_val){

  # Generate means for each regime
  mu <- runif(N, min = 0, max = 5)

  ## Generate lag-1 regression matrix (A) -- this is the temporal network
  A_list <- generate_A(N, Density, min_edg_val, max_edg_val)
  A <- A_list$A

  # Calculate average controllability for A (temporal network only)
  AC <- average_controllability(A)


  ## Generate contemporaneous network dynamics
  # Make precision matrix K -- this is the contemporaneous network
  K_list <- generate_K(N, Density, min_edg_val, max_edg_val)

  # Calculate sigma from K
  K <- K_list$K
  sigma <- solve(K)

  # Check if all off-diagonal non-zero elements of K are above threshold
  # (K is symmetric, so the lower triangle represents all off-diagonal edges)
  K_offdiag <- K[lower.tri(K, diag = FALSE)]
  K_thresh <- if (any(K_offdiag != 0)) min(abs(K_offdiag[K_offdiag != 0])) else NA_real_

  # Calculate empirical off-diagonal density of K
  K_dens <- mean(K[lower.tri(K, diag = FALSE)] != 0)
  K_weigh_dens <- mean(abs(K[lower.tri(K, diag = FALSE)]))


  return(list(mu = mu,
              A = A_list$A,
              A_sd = A_list$sd_A,
              A_strength = A_list$A_str,
              A_mean_strength = A_list$A_mean_str,
              A_density = A_list$A_dens,
              A_weighted_density = A_list$A_weigh_dens,
              A_stability = A_list$A_stability,
              AC = AC,
              sigma = sigma,
              K = K_list$K,
              K_pos.definit = K_list$K_posdef,
              K_thresh = K_thresh,
              K_density = K_dens,
              K_weighted_density = K_weigh_dens
              ))
  }
