## Function to generate random sets of dynamic matrices (sigma, kappa, Wtemp, Wcont, means)
## for various numbers of regimes. 
## Generates time-series with random starting state vector from dynamic matrices.
##
## Parameters are:
## Density <- c(0.2, 0.5, 0.8)      : Edge density 
## min_edg_val <- 0.08              : Minimal absolute edge value
## min_edg_val <- 1                 : Maximal absolute edge value
## N <- c(5, 8, 10)                 : No. of nodes
## M <- c(1, 2, 5)                  : No. of regimes
## T <- 1000                        : No. of measurements
## n_ts <- 50                       : No. of time-series to be generated
## warmup = 100                     : Warmup vectors (omitted from time-esries)
## totTime <- T + warmup            : Set total measurements for time-series incl. warmup
## mean_rep <- 10                   : Set mean length of regime repetitions in Rseq
## sd_rep <- 2                      : Set sd for mean length of regime repetitions in Rseq


# Load libraries
library(netcontrol)
library(graphicalVAR)
library(abind)
library(mvtnorm)
library(progress)
source("generate_netdyn.R")
source("generate_transmat.R")

generate_timeseries <- function(Density,
                                N,
                                M,
                                T
                                n_ts,
                                warmup,
                                totTime,
                                mean_rep,
                                sd_rep,
                                min_edg_val,
                                max_edg_val,
                                min_kappa_edg_val,
                                max_kappa_edg_val) {
  
  # Setup progressbar
  print("Generating networks", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )
  
  
  # Generate set of dynamics for each regime
  #set.seed(123)
  DynamicsMatrices_list <- list()
  for (i in seq_along(Density)) {
    dyn_level2 <- list()
    for (j in seq_along(N)) {
      dyn_level3 <- list()
      for (k in 1:n_ts) {
        dyn_level4 <- list()
        for (l in seq_along(M)) {
          dyn_level5 <- list()
          for (m in 1:M[l]) {
            repeat {
              W <- generate_netdyn(N[j], Density[i], min_edg_val, max_edg_val)
              if (all(abs(W[["Wtemp"]][W[["Wtemp"]] != 0]) >= min_edg_val) &&
                  all(abs(W[["Wcont"]][W[["Wcont"]] != 0]) >= min_edg_val) &&
                  any(W[["Wcont"]] != 0)) {
                break
                }
             }
            dyn_level5[[paste0("Regime", m)]] <- W
          }
          pb$tick()
          
          dyn_level4[[paste0(M[l], "_Regimes")]] <- dyn_level5[1:M[l]]
        }
        dyn_level3[[paste0("Set_no.", k)]] <- dyn_level4
      }
      dyn_level2[[paste0(N[j], "_Nodes")]] <- dyn_level3
    }
    DynamicsMatrices_list[[paste0("Density_", Density[i])]] <- dyn_level2
  }
 
  
  # Setup progressbar
  print("Generating timeseries", quote = FALSE)
  pb <- progress_bar$new(
    format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
    total = length(T) * length(Density) * length(N) * length(M) * n_ts,
    clear = FALSE,
    width = 80
  )
  
  # Generate timeseries data
  #set.seed(456)
  Timeseries_data <- list()
  for (t in seq_along(totTime)) {
    ts_level1 <- list()
    for (i in seq_along(Density)) {
      ts_level2 <- list()
      for (j in seq_along(N)) {
        ts_level3 <- list()
        for (k in seq_along(M)) {
          ts_level4 <- list()
          Rseq <- numeric()
          all_Wtemp_sd <- numeric()
          for (l in 1:n_ts) {
            
            # Generate initial vector for timeseries
            init <- runif(N[j], min = 0, max = 5)
            
            # Get dynamics for current timeseries
            dynamics <- DynamicsMatrices_list[[i]][[j]][[l]][[k]]
            
            
            # ## Generate TS from Rseq
            # # Generate sequence of regime indices
            # Rseq <- generate_Rseq(M[k], totTime[t], mean_rep, sd_rep)
            # 
            # # Calculate transition matrix from Rseq
            # transmat <- calc_transmat(Rseq, M[k])
            # 
            # # Initialize state matrix
            # X <- matrix(init, nrow = totTime[t], ncol = N[j])
            # 
            # # Get dynamics for current regime
            # for (m in 2:totTime[t]) {
            # curreg_mu <- dynamics[[Rseq[m]]][["mu"]]
            # # curreg_W_temp <- dynamics[[Rseq[m]]][["W_temp"]]
            # curreg_Beta <- dynamics[[Rseq[m]]][["Beta"]]
            # curreg_sigma <- dynamics[[Rseq[m]]][["sigma"]]

            
            
            ## Generate TS from Transmat
            # Generate transmat
            transmat <- generate_transmat(M[k])

            # Initialize state matrix
            X <- matrix(init, nrow = totTime[t], ncol = N[j])
            for (m in 2:totTime[t]) {

              # Set starting regime index
              if (m == 2) {
                reg_index <- sample.int(M[k], 1, replace = TRUE)
                Rseq <- reg_index
              }

              # Determine regime in next timestep from current timestep and transmat
              reg_index <- sample.int(M[k], 1, prob = transmat[reg_index, ])
              Rseq <- c(Rseq, reg_index)

              # Get dynamics for current regime
              curreg_mu <- dynamics[[reg_index]][["mu"]]
              curreg_W_temp <- dynamics[[reg_index]][["W_temp"]]
              curreg_Beta <- dynamics[[reg_index]][["Beta"]]
              curreg_sigma <- dynamics[[reg_index]][["sigma"]]
        
        
              # Generate next state vector from previous x, Beta and residuals from sigma
              X[m, ] <- curreg_mu + 
                t(curreg_Beta %*% (X[m - 1, ] - curreg_mu)) + 
                mvtnorm::rmvnorm(1, rep(0, N[j]), curreg_sigma)
              
            }
            pb$tick()
            
            ts_level4[[paste0("Timeseries_", l)]] <- list("Timeseries" = X[-seq_len(warmup), , drop = FALSE],
                                                          "Regime_sequence" = Rseq,
                                                          "Regime_dynamics" = dynamics,
                                                          "Transmat" = transmat)
          }
          ts_level3[[paste0(M[k], "_Regimes")]] <- ts_level4
        }
        ts_level2[[paste0(N[j], "_Nodes")]] <- ts_level3
      }
      ts_level1[[paste0("Density_", Density[i]*100, "%")]] <- ts_level2
    }
    Timeseries_data[[paste0(T[t], "_Timesteps")]] <- ts_level1
  }
  return(Timeseries_data = Timeseries_data)
}

