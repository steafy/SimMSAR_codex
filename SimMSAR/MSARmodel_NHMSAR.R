# Generate dummy data
library(netcontrol)
library(huge)
source("generate_netdyn.R")
source("generate_timeseries.R")
source("generate_Rseq.R")
source("calc_transmat.R")

# Load functions from NHMSAR
library(NHMSAR)
source("fit.MSAR_revised.R")
source("init.theta.MSAR_revised.R")
source("as.thetaMSAR_revised.R")
source("Mstep.hh.lasso.MSAR_patched.R")
source("Mstep.hh.reduct.MSAR_patched.R")

# Load function to asign and compare regimes
source("asign_regimes.R")

Density <- c(0.2, 0.5, 0.8)
min_edg_val <- 0
max_edg_val <- 1
M <- c(2, 3, 4)
N <- c(3, 4, 5)
mean_rep <- 10
sd_rep <- 2
n_ts <- 50
T <- 5000
warmup <- 500
totTime <- T + warmup


Timeseries_data <- generate_timeseries(Density = Density,
                                  min_edg_val = min_edg_val,
                                  max_edg_val = max_edg_val,
                                  M = M, 
                                  N = N, 
                                  sd_rep = sd_rep, 
                                  mean_rep = mean_rep, 
                                  warmup = warmup,
                                  T = T,
                                  totTime = totTime,
                                  n_ts = n_ts
                                  )

# Extract ts, Rseq, dynamics. Normalise ts
dummy_ts <- dummy_data[["Density_30%"]][[paste0(N, "_Nodes")]][[paste0(M, "_Regimes")]][["Timeseries_1"]][["Timeseries"]]
dummy_R_seq <- dummy_data[["Density_30%"]][[paste0(N, "_Nodes")]][[paste0(M, "_Regimes")]][["Timeseries_1"]][["Regime_sequence"]]
dummy_dynamics <- dummy_data[["Density_30%"]][[paste0(N, "_Nodes")]][[paste0(M, "_Regimes")]][["Timeseries_1"]][["Regime_dynamics"]]
dummy_ts_norm <- huge.npn(dummy_ts)


# Make array from data
T <- nrow(dummy_ts_norm)  # No. of timesteps
d <- ncol(dummy_ts_norm)  # No. of variables
N.samples <- 1            # No. of measurement series
dummy_ts_norm_array <- array(data = dummy_ts_norm, dim = c(T, N.samples, d))


# Set parameters for NHMSAR
oder <- 1

print( "START" )
# Initialize model
model_init <- init.theta.MSAR_revised(dummy_ts_norm_array, M = M, order = 1, label = "HH")

#debug(fit.MSAR_revised, browser)
# Fit model
model_fit <- fit.MSAR_revised(data = dummy_ts_norm_array, theta = model_init, MaxIter = 1000, penalty = "LASSO")
print( "FERTIG" )


Sys.sleep(3)

# Delete coefficients in estimated Wtemp below threshold (0.05)
for(i in 1:M) {
  for(j in 1:N) {
    for(k in 1:N) {
      if(abs(model_fit[["theta"]][["A"]][[i]][["A1"]][j,k]) < 0.05){
        model_fit[["theta"]][["A"]][[i]][["A1"]][j,k] <- 0
        }
      if(abs(model_fit[["theta"]][["sigma"]][[i]][j,k]) < 0.05){
        model_fit[["theta"]][["sigma"]][[i]][j,k] <- 0
      }
    }
  }
}


# Extract original Wtemp
get_org_Wtemp <- function(dummy_dynamics) {
  vectors <- list()
  for (regime in names(dummy_dynamics)) {
    vec_name <- paste(regime, "W_temp", sep = "_")
    vectors[[vec_name]] <- as.vector(dummy_dynamics[[regime]]$W_temp)
  }
  return(vectors)
}

# Extract estimated Wtemp
get_est_Wtemp  <- function(model_fit) {
  vectors <- list()
  for (regime in names(model_fit$theta$A)) {
    vec_name <- paste(regime, "W_temp", sep = "_")
    vectors[[vec_name]] <- as.vector(model_fit$theta$A[[regime]]$A1)
  }
  return(vectors)
}

# Make matrix of vectors from original and estimated Wtemp
vecs_org_Wtemp <- get_org_Wtemp(dummy_dynamics)
vecs_est_Wtemp <- get_est_Wtemp(model_fit)

# Calculate correlations for each pair of original and estimated Wtemp
cor_results <- outer(names(vecs_org_Wtemp), names(vecs_est_Wtemp), Vectorize(function(n1, n2) {
  cor(vecs_org_Wtemp[[n1]], vecs_est_Wtemp[[n2]])
}))


asigned_regimes <- asign_regimes(cor_results)





