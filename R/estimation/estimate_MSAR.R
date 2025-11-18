# Dependencies are loaded centrally via R/dependencies.R
# Required packages: huge, NHMSAR, dplyr, progress

# Load functions from NHMSAR
source("R/estimation/fit.MSAR_revised_2.R")
source("R/estimation/init.theta.MSAR_revised_2.R")
source("R/estimation/init_and_fit.MSAR_Lasso_2.R")

# Load functions to assign and compare regime dynamics
source("R/utils/asign_regimes.R")
source("R/utils/senspec.R")
source("R/utils/summarize_cor.R")
source("R/utils/calculate_MAE.R")



## Estimate network dynamics from timeseries data using NHMSAR
estimate_MSAR <- function(Density, N, M, T, n_ts, order, MaxIter, verbose, min_edg_val, Timeseries_data) {

# Setup progressbar
print("Estimating MSAR models from timeseries", quote = FALSE)
pb <- progress_bar$new(
  format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining",
  total = length(T) * length(Density) * length(N) * length(M) * n_ts,
  clear = FALSE,
  width = 80
)

# markov regime-switiching autoregression models with expectation maximization method
MSAR_dynamics_list <- list()
for (t in seq_along(T)) {
  MSAR_level1 <- list()
  for (i in seq_along(Density)) {
    MSAR_level2 <- list()
    for (j in seq_along(N)) {
      MSAR_level3 <- list()
      for (k in seq_along(M)) {
        MSAR_level4 <- list()
        vec_all_Wtemp_cor <- numeric()
        vec_all_Wtemp_ac_cor <-numeric()
        vec_all_Wcont_cor <- numeric() 
        vec_all_Wcont_ac_cor <- numeric()
        vec_all_Wtemp_sen <- numeric()
        vec_all_Wtemp_spec <- numeric()
        vec_all_Wcont_sen <- numeric()
        vec_all_Wcont_spec <- numeric()
        vec_all_Wtemp_MAE <- numeric()
        vec_all_Wcont_MAE <- numeric()
        for (l in 1:n_ts) {
          
          ### Jan
          # reset !!!
          #set.seed(seed)
          
          loc <-
            paste("i=", i, ", j=", j, ", k=", k, ", l=", l, sep = "")
          
          if (verbose)
          {
            print("", quote = FALSE)
            print("", quote = FALSE)
            print(paste("Timeseries:", loc), quote = FALSE)
            print(paste("fit parms:", app_parms), quote = FALSE)
            print("", quote = FALSE)
          }
          ### Jan end
          
          # Get complete data for current timeseries
          current_data <- Timeseries_data[[t]][[i]][[j]][[k]][[l]]
          
          # Extract timeseries data
          current_ts <- current_data[["Timeseries"]]
          
          # Extract regimes data
          current_regimes <- current_data[["Regime_dynamics"]]
          
          # Normalize timeseries data
          current_ts_norm <- huge.npn(current_ts, verbose = verbose)
          #current_ts_norm <- current_ts
          
          # Make array from timeseries data
          timesteps <- T[t]                   # Set no. of timesteps
          d <- ncol(current_ts_norm)          # Set no. of nodes
          N.samples <- 1                      # Set no. of time series
          current_ts_norm_array <- array(data = current_ts_norm,
                                         dim = c(timesteps, N.samples, d))
          
          # Estimate MSAR models
          result <- init_and_fit.MSAR_Lasso_2(
            data = current_ts_norm_array,
            M = M[k],
            order = order,
            MaxIter = MaxIter,
            retry = 5,
            verbose = verbose
          )
          
          model_fit <- result[["fit"]]
          error <- result[["error"]]
          
          # browser()
          
          if (is.null(model_fit)) {
            message("ignore fit! ", "\n\tTimeseries_data: ", loc, "\n\terror: ", error)
            next
          }
          
          Sys.sleep(1)
          
          # Extract estimated Betas for all regimes from MSAR model
          est_Betas <- model_fit[["theta"]][["A"]]
          
          # Extract estimated sigma for all regimes
          est_sigmas <- model_fit[["theta"]][["sigma"]]
          
          # Make kappa from estimated sigma for all regimes
          est_kappas <- list()
          for(m in 1:M[k]) {
            est_kappas[[paste0("Regime", m)]] <- solve(est_sigmas[[m]])
          }
          
          
          # Make estimated Wcont from estimated kappa
          est_Wconts <- list()
          for(m in 1:M[k]) {
            est_Wcont <- matrix(0, nrow = N[j], ncol = N[j])
            est_kappa <- est_kappas[[m]]
            for(n in 1:N[j]) {
              for(o in 1:N[j]) {
                if(n != o) {
                  est_Wcont[n, o] <- -est_kappa[n, o] / sqrt(est_kappa[n, n] * est_kappa[o, o])
                  est_Wconts[[paste0("Regime", m)]] <- est_Wcont
                }
              }
            }
          }
          
          # Make estimated Wtemps from estimates of Beta, kappa and sigma
          est_Wtemps <- list()
          for(m in 1:M[k]) {
            est_Wtemp <- matrix(0, nrow = N[j], ncol = N[j])
            est_Beta <- est_Betas[[m]][["A1"]]
            est_kappa <- est_kappas[[m]]
            est_sigma <- est_sigmas[[m]]
            for(n in 1:N[j]) {
              for(o in 1:N[j]) { # raus: if(n != o) {
                # if(n != o) {
                #est_Wtemp[n, o] <- est_Beta[n, o] / sqrt(est_sigma[i, i] %*% est_kappa[j ,j] + est_Beta[i, j]^2)
                est_Wtemp[n, o] <- est_Beta[n, o] / sqrt(est_sigma[n, n] %*% est_kappa[o, o] + est_Beta[n, o]^2)
                est_Wtemps[[paste0("Regime", m)]] <- est_Wtemp
                # }
              }
            }
          }
          
          # Set coefficients in estimated Wtemp and Wcont below threshold to 0
          for(m in 1:M[k]) {
            for(n in 1:N[j]) {
              for(o in 1:N[j]) {
                if(abs(est_Wtemps[[m]][n, o]) < min_edg_val){
                  est_Wtemps[[m]][n, o] <- 0
                }
                if(abs(est_Wconts[[m]][n, o]) < min_edg_val){
                  est_Wconts[[m]][n, o] <- 0
                }
              }
            }
          }
          
          # Extract original Wtemp
          get_org_Wtemp <- function(current_regimes) {
            vectors <- list()
            for (regime in names(current_regimes)) {
              vec_name <- paste(regime, "Wtemp", sep = "_")
              vectors[[vec_name]] <- as.vector(current_regimes[[regime]]$Wtemp)
            }
            return(vectors)
          }
          
          # Extract estimated Wtemp
          get_est_Wtemp  <- function(model_fit) {  # vorher function(model_fit)
            vectors <- list()
            for (regime in names(est_Wtemps)) {
              vec_name <- paste(regime, "Wtemp", sep = "_")
              vectors[[vec_name]] <- as.vector(est_Wtemps[[regime]])
            }
            return(vectors)
          }
          
          # Make matrix of vectors from original and estimated Wtemp
          vecs_org_Wtemp <- get_org_Wtemp(current_regimes)
          vecs_est_Wtemp <- get_est_Wtemp(model_fit)
          
          # Check if estimated Wtemp has sd == 0
          if (any(sapply(vecs_est_Wtemp,function(vec) sd(vec) == 0))) {
            next
          }
          
          # Calculate correlations for each pair of original and estimated Wtemp
          cor_results <- outer(names(vecs_org_Wtemp),
                               names(vecs_est_Wtemp),
                               Vectorize(function(n1, n2) {
                                 cor(vecs_org_Wtemp[[n1]],
                                     vecs_est_Wtemp[[n2]])
                               }
                               )
          )
          
          
          # Asign original and estimated regimes based on highest correlations of Wtemp
          if (M[k] > 1) {
            asigned_regimes <- asign_regimes(cor_results)
            
            # Add highest correlations to overall vector
            vec_all_Wtemp_cor <- c(vec_all_Wtemp_cor,
                                   asigned_regimes[ ,3]
            )
            
            # Make regimewise pairs of original and estimated Wtemp based on highest correlations
            pair_Wtemp <- list()
            for (m in 1:nrow(asigned_regimes)) {
              orig_Wtemp = current_regimes[[asigned_regimes[[m, 1]]]][["Wtemp"]]
              est_Wtemp <- est_Wtemps[[asigned_regimes[[m, 2]]]]
              Wtemp_senspec <- senspec(orig_Wtemp,
                                       est_Wtemp
              )
              vec_all_Wtemp_sen <- c(vec_all_Wtemp_sen,
                                     Wtemp_senspec[["sensitivity"]]
              )
              vec_all_Wtemp_spec <- c(vec_all_Wtemp_spec,
                                      Wtemp_senspec[["specificity"]]
              )
              MAE_Wtemp <- calculate_MAE(orig_Wtemp, est_Wtemp)
              vec_all_Wtemp_MAE <- c(vec_all_Wtemp_MAE, MAE_Wtemp)
              est_Wtemp_ac <- ave_control_centrality(est_Wtemp)
              orig_Wtemp_ac <- current_regimes[[asigned_regimes[[m, 1]]]][["Wtemp_ac"]]
              cor_Wtemp_ac <- cor(orig_Wtemp_ac,
                                  est_Wtemp_ac,
                                  method = "pearson"
              )
              orig_Wcont <- current_regimes[[asigned_regimes[[m, 1]]]][["Wcont"]]
              est_Wcont <- est_Wconts[[asigned_regimes[[m, 2]]]]
              Wcont_senspec <- senspec(orig_Wcont,
                                       est_Wcont
              )
              vec_all_Wcont_sen <- c(vec_all_Wcont_sen,
                                     Wcont_senspec[["sensitivity"]]
              )
              vec_all_Wcont_spec <- c(vec_all_Wcont_spec,
                                      Wcont_senspec[["specificity"]]
              )
              vec_all_Wtemp_ac_cor <- c(vec_all_Wtemp_ac_cor, cor_Wtemp_ac
              )
              cor_Wcont <- cor(as.vector(orig_Wcont),
                               as.vector(est_Wcont),
                               method = "pearson"
              )
              MAE_Wcont <- calculate_MAE(orig_Wcont, est_Wcont)
              vec_all_Wcont_cor <- c(vec_all_Wcont_cor,
                                     cor_Wcont
              )
              vec_all_Wcont_MAE <- c(vec_all_Wcont_MAE, MAE_Wcont)
              est_Wcont_ac <- ave_control_centrality(est_Wcont)
              orig_Wcont_ac <- current_regimes[[asigned_regimes[[m, 1]]]][["Wcont_ac"]]
              cor_Wcont_ac <- cor(orig_Wcont_ac,
                                  est_Wcont_ac,
                                  method = "pearson"
              )
              vec_all_Wcont_ac_cor <- c(vec_all_Wcont_ac_cor,
                                        cor_Wcont_ac
              )
              
              pair_Wtemp[[paste0("Regime", m)]] <- list("orig. Wtemp" = orig_Wtemp,
                                                        "est. Wtemp" = est_Wtemp,
                                                        "sen. Wtemp" = Wtemp_senspec[["sensitivity"]],
                                                        "spec. Wtemp" = Wtemp_senspec[["specificity"]],
                                                        "corr. Wtemp" = asigned_regimes[[m, 3]],
                                                        "MAE Wtemp" = MAE_Wtemp,
                                                        "orig. Wtemp ac" = orig_Wtemp_ac,
                                                        "est. Wtemp ac" = est_Wtemp_ac,
                                                        "corr. Wtemp ac" = cor_Wtemp_ac,
                                                        "orig. Wcont" = orig_Wcont,
                                                        "est. Wcont" = est_Wcont,
                                                        "sen. Wcont" = Wcont_senspec[["sensitivity"]],
                                                        "spec. Wcont" = Wcont_senspec[["specificity"]],
                                                        "corr. Wcont" = cor_Wcont,
                                                        "MAE Wcont" = MAE_Wcont,
                                                        "orig. Wcont ac" = orig_Wcont_ac,
                                                        "est. Wcont ac" = est_Wcont_ac,
                                                        "corr. Wcont ac" = cor_Wcont_ac
              )
            }
          } else {
            pair_Wtemp <- list()
            orig_Wtemp <- current_regimes[["Regime1"]][["Wtemp"]]
            est_Wtemp <- est_Wtemps[["Regime1"]]
            Wtemp_senspec <- senspec(orig_Wtemp,
                                     est_Wtemp
            )
            vec_all_Wtemp_sen <- c(vec_all_Wtemp_sen,
                                   Wtemp_senspec[["sensitivity"]]
            )
            vec_all_Wtemp_spec <- c(vec_all_Wtemp_spec,
                                    Wtemp_senspec[["specificity"]]
            )
            cor_Wtemp <- cor(as.vector(orig_Wtemp),
                             as.vector(est_Wtemp),
                             method = "pearson"
            )
            MAE_Wtemp <- calculate_MAE(orig_Wtemp, est_Wtemp)
            vec_all_Wtemp_MAE <- c(vec_all_Wtemp_MAE, MAE_Wtemp)
            est_Wtemp_ac <- ave_control_centrality(est_Wtemp)
            orig_Wtemp_ac <- current_regimes[["Regime1"]][["Wtemp_ac"]]
            cor_Wtemp_ac <- cor(orig_Wtemp_ac,
                                est_Wtemp_ac,
                                method = "pearson"
            )
            orig_Wcont <- current_regimes[["Regime1"]][["Wcont"]]
            est_Wcont <- est_Wconts[["Regime1"]]
            cor_Wcont <- cor(as.vector(orig_Wcont),
                             as.vector(est_Wcont),
                             method = "pearson"
            )
            vec_all_Wtemp_cor <- c(vec_all_Wtemp_cor,
                                   cor_Wtemp
            )
            vec_all_Wtemp_ac_cor <- c(vec_all_Wtemp_ac_cor,
                                      cor_Wtemp_ac
            )
            vec_all_Wcont_cor <- c(vec_all_Wcont_cor,
                                   cor_Wcont
            )
            est_Wcont_ac <- ave_control_centrality(est_Wcont)
            orig_Wcont_ac <- current_regimes[["Regime1"]][["Wcont_ac"]]
            Wcont_senspec <- senspec(orig_Wcont,
                                     est_Wcont
            )
            vec_all_Wcont_sen <- c(vec_all_Wcont_sen,
                                   Wcont_senspec[["sensitivity"]]
            )
            vec_all_Wcont_spec <- c(vec_all_Wcont_spec,
                                    Wcont_senspec[["specificity"]]
            )
            cor_Wcont_ac <- cor(orig_Wcont_ac,
                                est_Wcont_ac,
                                method = "pearson"
            )
            vec_all_Wcont_ac_cor <- c(vec_all_Wcont_ac_cor,
                                      cor_Wcont_ac
            )
            MAE_Wcont <- calculate_MAE(orig_Wcont, est_Wcont)
            vec_all_Wcont_MAE <- c(vec_all_Wcont_MAE, MAE_Wcont)
            
            pair_Wtemp[[paste0("Regime", m)]] <- list("orig. Wtemp" = orig_Wtemp,
                                                      "est. Wtemp" = est_Wtemp,
                                                      "sen. Wtemp" = Wtemp_senspec[["sensitivity"]],
                                                      "spec. Wtemp" = Wtemp_senspec[["specificity"]],
                                                      "corr. Wtemp" = cor_Wtemp,
                                                      "MAE Wtemp" = MAE_Wtemp,
                                                      "orig. Wtemp ac" = orig_Wtemp_ac,
                                                      "est. Wtemp ac" = est_Wtemp_ac,
                                                      "corr. Wtemp ac" = cor_Wtemp_ac,
                                                      "orig. Wcont" = orig_Wcont,
                                                      "est. Wcont" = est_Wcont,
                                                      "sen. Wcont" = Wcont_senspec[["sensitivity"]],
                                                      "spec. Wcont" = Wcont_senspec[["specificity"]],
                                                      "corr. Wcont" = cor_Wcont,
                                                      "MAE Wcont" = MAE_Wcont,
                                                      "orig. Wcont ac" = orig_Wcont_ac,
                                                      "est. Wcont ac" = est_Wcont_ac,
                                                      "corr. Wcont ac" = cor_Wcont_ac
            )
          }
          
          # Update progressbar
          pb$tick()
          
          MSAR_level4[[paste0("Timeseries_", l)]] <- pair_Wtemp
        }
        
        
        vec_all_Wcont_cor <- na.omit(vec_all_Wcont_cor)
        vec_all_Wcont_ac_cor <- na.omit(vec_all_Wcont_ac_cor)
        
        stats <- list(Wtemp_corr = summarize_cor(vec_all_Wtemp_cor),
                      Wtemp_MAE = summarize_cor(vec_all_Wtemp_MAE),
                      Wtemp_sensitivity = summarize_cor(vec_all_Wtemp_sen),
                      Wtemp_specificity = summarize_cor(vec_all_Wtemp_spec),
                      Wtemp_ac_corr = summarize_cor(vec_all_Wtemp_ac_cor),
                      Wcont_corr = summarize_cor(vec_all_Wcont_cor),
                      Wcont_MAE = summarize_cor(vec_all_Wcont_MAE),
                      Wcont_sensitivity = summarize_cor(vec_all_Wcont_sen),
                      Wcont_specificity = summarize_cor(vec_all_Wcont_spec),
                      Wcont_ac_corr = summarize_cor(vec_all_Wcont_ac_cor)
        )
        
        MSAR_level3[[paste0(M[k], "_Regimes")]] <- list(MSAR_models = MSAR_level4, Stats = stats)
      }
      MSAR_level2[[paste0(N[j], "_Nodes")]] <- MSAR_level3  
    }
    MSAR_level1[[paste0("Density_", Density[i]*100, "%")]] <- MSAR_level2
  }
  MSAR_dynamics_list[[paste0(T[t], "_Timesteps")]] <- MSAR_level1
}

return(MSAR_dynamics_list)

}
