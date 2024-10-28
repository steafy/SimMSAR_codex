# Load functions for data generation
source("generate_timeseries.R")

# Load functions for model estimation
source("estimate_MSAR.R")

# Unused functions
#source("as.thetaMSAR_revised.R")
#source("Mstep.hh.lasso.MSAR_patched.R")
#source("Mstep.hh.reduct.MSAR_patched.R")
# source("generate_Rseq.R")
# source("calc_transmat.R")


# Set parameter values
Density <- c(0.25, 0.5, 0.75)
min_edg_val <- 0.1
max_edg_val <- 1
min_kappa_edg_val <- 1
max_kappa_edg_val <- 1
M <- c(1,2,3,4)
N <- c(4,6,10)
# mean_rep <- 10 # for Rseq
# sd_rep <- 3    # for Rseq
n_ts <- 30
T <- c(100,200,500,1500,2000,2500,3000,3500)
warmup <- 50
totTime <- T + warmup
order <- 1
MaxIter <- 200 ### !!!
verbose <- FALSE
# seed <- 83742



app_parms <- paste(
  "Density=", Density, #
  ", M=", M, #
  ", N=", N, #
  ", n_ts=", n_ts, #
  ", T=", T, #
  ", order=", order, #
  #", seed=", seed, #
  sep = ""
  )

print(paste("fit parms:", app_parms), quote = FALSE)


## Generate timeseries data (timeseries, regimes, regime sequence, transition matrix)
Timeseries_data <- generate_timeseries(Density = Density,
                                       min_edg_val = min_edg_val,
                                       max_edg_val = max_edg_val,
                                       min_kappa_edg_val = min_kappa_edg_val,
                                       max_kappa_edg_val = max_kappa_edg_val,
                                       M = M, 
                                       N = N, 
                                       sd_rep = sd_rep, 
                                       mean_rep = mean_rep, 
                                       warmup = warmup,
                                       T = T,
                                       totTime = totTime,
                                       n_ts = n_ts
                                       )


## Estimate network dynamics from timeseries data using NHMSAR
MSAR_dynamics_list <- estimate_MSAR(Density = Density,
                                    M = M, 
                                    N = N, 
                                    T = T,
                                    n_ts = n_ts,
                                    order = order,
                                    MaxIter = MaxIter,
                                    verbose = verbose,
                                    min_edg_val = min_edg_val,
                                    Timeseries_data
                                    )
