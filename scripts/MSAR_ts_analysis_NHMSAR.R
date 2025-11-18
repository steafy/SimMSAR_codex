# Load functions for data generation
source("R/generation/generate_timeseries.R")

# Load functions for model estimation
source("R/estimation/estimate_MSAR.R")

# Unused functions
#source("as.thetaMSAR_revised.R")
#source("Mstep.hh.lasso.MSAR_patched.R")
#source("Mstep.hh.reduct.MSAR_patched.R")
# source("generate_Rseq.R")
# source("calc_transmat.R")


# Set parameter values
Density <- c(0.25)  # Set network density
min_edg_val <- 0.05           # Set min. edge value
max_edg_val <- 1               # Set max edge value
M <- c(2)             # Set number of regimes
N <- c(4)                # Set number of nodes
# mean_rep <- 10 # for Rseq
# sd_rep <- 3    # for Rseq
n_ts <- 30                     # Set number of timeseries for each combination of factors
T <- c(3500) # Set number of time steps
warmup <- 50                   # Set number of time steps to warm up (bein omitted)
remain_lower <- 0.33            # Set lower bound for probability to stay in regime (0.92)
remain_upper <- 0.66            # Set upper bound for probability to stay in regime (0.95)
totTime <- T + warmup
order <- 1                     # Set lag of autoregression
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
                                       M = M, 
                                       N = N, 
                                       sd_rep = sd_rep, 
                                       mean_rep = mean_rep, 
                                       warmup = warmup,
                                       T = T,
                                       totTime = totTime,
                                       n_ts = n_ts,
                                       remain_lower = remain_lower,
                                       remain_upper = remain_upper
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


#saveRDS(Timeseries_data, "Timeseries_data.rds")

#saveRDS(MSAR_dynamics_list, "MSAR_models.rds")
