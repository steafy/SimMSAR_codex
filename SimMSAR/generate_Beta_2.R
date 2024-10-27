# install.packages("memoise")
library(memoise)

# source("check_stability.R")

## Function to generate temporal network dynamics
#
# Generates a random (NxN) Matrix "Beta" which is:
# Stable: max(abs(eigen(Beta)$values)) < 1,
# Has a density near the given "Density": mean(Beta != 0) ~= Density,
# And coefficients between "min_edg_val" and "max_edg_val".
#
# Test: generate_Beta_2(10,0.9,0.1,1,TRUE)
generate_Beta_2 <-
  function(N,
           Density,
           min_edg_val,
           max_edg_val,
           verbose = FALSE) {
    repeat {
      Beta <- generate_Beta_maybe_unstable (N,
                                            Density,
                                            min_edg_val,
                                            max_edg_val,
                                            verbose = FALSE)
      
      # Check stability of Beta
      Beta_stability <- check_stability_2(Beta)
      
      # Calculate empirical density of Beta
      Beta_dens <- mean(Beta != 0)
      
      if (verbose) {
        message(paste(
          "Beta_stability:",
          Beta_stability,
          ", Beta_dens:",
          Beta_dens
        ))
      }
      
      # browser()
      
      if (Beta_stability == TRUE) {
        # Beta is stable (AND the empirical dens is near "Density" (TODO))
        # Beta is good -> return
        
        # Calculate sd of Beta
        sd_Beta <- sd(as.vector(Beta))
        
        # Calculate node strength of Beta
        Beta_str <- rowSums(Beta)
        Beta_mean_str <- mean(Beta_str)
        
        # Calculate empirical density of Beta (rest)
        Beta_weigh_dens <- mean(abs(Beta))
        
        if (verbose) {
          message(paste("->return"))
        }
        
        return(
          list(
            Beta = Beta,
            Beta_stability = Beta_stability,
            # TODO "YES" / "NO" ?
            sd_Beta = sd_Beta,
            Beta_str = Beta_str,
            Beta_mean_str = Beta_mean_str,
            Beta_dens = Beta_dens,
            Beta_weigh_dens = Beta_weigh_dens
          )
        )
        
      }
      
      # here: Beta is not good -> repeat
      if (verbose) {
        message(paste("->repeat"))
      }
      
    } # end repeat
  }

# Generates a maybe unstable Beta with a maybe wrong "Density".
generate_Beta_maybe_unstable <-
  function(N,
           Density,
           min_edg_val,
           max_edg_val,
           verbose = FALSE) {
    # we generate Beta with generate_Beta_intern(...) and a "num_nonzero"
    # of (Density * N * N) + avg_cuts; we expect that generate_Beta_intern(...)
    # cuts "avg_cuts" of the "num_nonzero"s to zero, and so expect that the
    # returned Beta at the end has a "Density" of the given "Density".
    
    avg_cuts <-
      cached_avg_generate_Beta_cuts(N, Density, min_edg_val, max_edg_val)
    
    num_nonzero <- round(Density * N * N) + avg_cuts
    
    results <-
      generate_Beta_intern(N, num_nonzero, min_edg_val, max_edg_val, verbose)
    
    # Beta <- results$Beta
    # cut_count <- results$cut_count
    
    return (results$Beta)
  }

# Generates a Beta (maybe unstable) and returns these "Beta" and a "cut_count";
# NOTE: the second parameter is "num_nonzero" and NOT "Density" !!!
# The returned "Beta" may have less than "num_nonzero"s non-zeros.
# The returned "cut_count" reflects these lack of non-zeros.
generate_Beta_intern <-
  function(N,
           num_nonzero,
           min_edg_val,
           max_edg_val,
           verbose = FALSE) {
    # Generate Beta with non-zero elements between min_edg_val & max_edg_val
    Beta <- matrix(0, nrow = N, ncol = N)
    # num_nonzero == round(Density * N * N) + ...
    indices <- sample(N * N, num_nonzero, replace = FALSE)
    signs <- sample(c(-1, 1), num_nonzero, replace = TRUE)
    values <- runif(num_nonzero, min_edg_val, max_edg_val) * signs
    Beta[indices] <- values
    
    # Scale Beta to ensure stability
    eigvals <- eigen(Beta)$values
    spectral_radius <- max(abs(eigvals))
    scaling_factor <-
      ifelse(spectral_radius >= 1, 0.99 / spectral_radius, 1)
    Beta <- Beta * scaling_factor
    
    # here we have a "good" Beta but the koeffs may be to small (<min_edg_val).
    # adjust_non_zero_elements(...) cuts these to zeros which may result
    # in an unstable Beta and of cause results in a less densed Beta ...
    
    results <-
      adjust_non_zero_elements(Beta, min_edg_val, max_edg_val, verbose)
    
    # Beta <- results$Beta
    # cut_count <- results$cut_count
    
    return(results)
  }

# Adjust non-zero elements to between min_edg_val & max_edg_val
adjust_non_zero_elements <-
  function(Beta, min_edg_val, max_edg_val, verbose = FALSE) {
    cut_count <- 0
    non_zero_elements <- which(Beta != 0, arr.ind = TRUE)
    
    for (idx in 1:nrow(non_zero_elements)) {
      i <- non_zero_elements[idx, 1]
      j <- non_zero_elements[idx, 2]
      
      beta_i_j <- Beta[i, j]
      abs_beta_i_j <- abs(Beta[i, j])
      
      if (abs_beta_i_j > max_edg_val) {
        # case never actually occurs if max_edg_val == 1
        Beta[i, j] <- sign(beta_i_j) * max_edg_val
      }
      else if (abs_beta_i_j < min_edg_val) {
        if (verbose) {
          message(paste("Info:", "Cut Beta[", i, ",", j, "] to Zero, was", beta_i_j))
        }
        Beta[i, j] <- 0
        cut_count <- cut_count + 1
      }
    }
    
    return (list(Beta = Beta, cut_count = cut_count))
  }

# Check stability of Matrix
check_stability_2 <- function(Matrix) {
  eigen <- eigen(Matrix)
  stability <- ifelse(all(abs(eigen$values) < 1), TRUE, FALSE)
  return(stability)
}

# Determine the average “avg_cut_count” for the given parameters
# uncached, use cached_avg_generate_Beta_cuts !!
avg_generate_Beta_cuts <-
  function(N,
           Density,
           min_edg_val,
           max_edg_val,
           trials = 1000) {
    num_nonzero <- round(Density * N * N)
    
    cut_count <- 0
    for (i in 1:trials) {
      results <-
        generate_Beta_intern(N, num_nonzero, min_edg_val, max_edg_val, verbose = FALSE)
      cut_count <- cut_count + results$cut_count
    }
    avg_cut_count <- cut_count / trials
    return(avg_cut_count)
  }

# Determine the average “avg_cut_count” for the given parameters, cached!.
cached_avg_generate_Beta_cuts <- memoise(avg_generate_Beta_cuts)

