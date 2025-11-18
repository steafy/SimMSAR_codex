## Function to generate dynamics for each regime

# Dependencies are loaded centrally via R/dependencies.R
# Required packages: netcontrol, Matrix

# Load functions
source("R/generation/generate_random.R")
source("R/generation/generate_Beta.R")
source("R/generation/generate_kappa.R")
source("R/generation/check_stability.R")


generate_netdyn <- function (N, Density, min_edg_val, max_edg_val){

  # Generate means for each regime
  mu <- runif(N, min = 0, max = 5)
 
  ## Generate lag-1 regression matrix (Beta)
  Beta_list <- generate_Beta(N, Density, min_edg_val, max_edg_val)
  Beta <- Beta_list$Beta
  
  
  ## Generate contemporaneous network dynamics
  # Make precision matrix kappa
  kappa_list <- generate_kappa(N, Density, min_edg_val, max_edg_val)

  # Calculate sigma from kappa
  kappa <- kappa_list$kappa
  sigma <- solve(kappa)

  # Calculate matrix of partial correlations (Wcont) from kappa
  Wcont <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        Wcont[i, j] <- -kappa[i, j] / sqrt(kappa[i, i] * kappa[j, j])
      }
    }
  }

  # Check if all correlations in Wconr are above 0.1
  Wcont_thresh <- min(abs(Wcont[Wcont != 0]))
  
  # Check stablity of Wcont
  Wcont_stability <- check_stability(Wcont)

  # Calculate node strength
  Wcont_str <- rowSums(Wcont)
  Wcont_mean_str <- mean(Wcont_str)

  # Calculates empirical density of Wcont
  # Wcont_dens <- sum(Wcont[lower.tri(Wcont, diag = FALSE)] != 0) / ((N*(N-1))/2)
  Wcont_dens <- mean(Wcont[lower.tri(Wcont, diag = FALSE)] != 0)
  Wcont_weigh_dens <- mean(abs(Wcont[lower.tri(Wcont, diag = FALSE)]))


  # Calculate average controllability for Wcont
  Wcont_ac <- ave_control_centrality(Wcont)

  # Calculate Wtemp (parial correlation matrix) from Beta, kappa and sigma
  Wtemp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {   # if raus if (i != j) {
      #if (i != j) {
        Wtemp[i, j] <- Beta[i, j] / sqrt(sigma[i, i] %*% kappa[j ,j] + Beta[i, j]^2)
     # }
    }
  }

  # Calculate empirical density of Wtemp
  Wtemp_dens <- mean(Wtemp != 0)
  Wtemp_weigh_dens <- mean(abs(Wtemp))

  # Check if all correlations in Wtemp are above 0.1
  Wtemp_thresh <- min(abs(Wtemp[Wtemp != 0]))
  
  # Calculate average controllability for Wtemp
  Wtemp_ac <- ave_control_centrality(Wtemp)


  return(list(mu = mu,
              Beta = Beta_list$Beta,
              Beta_sd = Beta_list$sd_Beta,
              Beta_strength = Beta_list$Beta_str,
              Beta_mean_strength = Beta_list$Beta_mean_str,
              Beta_density = Beta_list$Beta_dens,
              Beta_weighted_density = Beta_list$Beta_weigh_dens,
              Beta_stability = Beta_list$Beta_stability,
              Wtemp = Wtemp,
              Wtemp_thresh = Wtemp_thresh,
              Wtemp_density = Wtemp_dens,
              Wtemp_weighted_density = Wtemp_weigh_dens,
              Wtemp_ac = Wtemp_ac,
              sigma = sigma,
              kappa = kappa_list$kappa,
              Wcont = Wcont,
              Wcont_thresh = Wcont_thresh,
              Wcont_strength = Wcont_str,
              Wcont_mean_strength = Wcont_mean_str,
              Wcont_density = Wcont_dens,
              Wcont_weighted_density = Wcont_weigh_dens,
              kappa_pos.definit = kappa_list$kappa_posdef,
              Wcont_stability = Wcont_stability,
              Wcont_ac = Wcont_ac
              ))
  }
