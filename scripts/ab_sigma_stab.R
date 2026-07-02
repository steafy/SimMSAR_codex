# =============================================================================
# A/B (uncommitted): Sigma stabilization vs baseline. ONE (cell, config) per
# invocation, so a segfault inside a degenerate baseline fit's compiled solve()
# (which tryCatch cannot catch) only loses that one combination, not the run.
#   Rscript scripts/ab_sigma_stab.R <cell_index> <config_name>
# Data is seeded PER CELL (independent of other cells), init PER (cell,rep), so
# every config sees identical data + identical init for a given cell. Writes
# output/diagnostics/ab_part_<cell>_<config>.csv. Combine with ab_sigma_combine.R.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
ci   <- as.integer(args[[1]])
cfnm <- args[[2]]

suppressMessages({
  source("R/dependencies.R")
  source("R/generation/generate_timeseries.R")
  source("R/estimation/fit_msar.R")
  source("R/estimation/init_theta_msar.R")
  source("R/estimation/match_regimes.R")
  source("R/utils/assign_regimes.R")
  source("R/utils/senspec.R")
  source("R/utils/undirected_metrics.R")
})

base_opts <- list(
  simmsar_lasso_engine="cvglmnet", simmsar_lasso_reselect=TRUE,
  simmsar_lasso_lambda="1se", simmsar_lasso_refit=TRUE, simmsar_lasso_adaptive=FALSE,
  simmsar_lasso_nfolds=5, simmsar_lasso_nlambda=50, simmsar_lasso_fixedfolds=TRUE,
  simmsar_lasso_reselect_iters=Inf, simmsar_lasso_reselect_every=1)

configs <- list(
  baseline  = list(simmsar_sigma_stab="none"),
  ridge_1e3 = list(simmsar_sigma_stab="ridge", simmsar_sigma_stab_lambda=1e-3),
  ridge_1e2 = list(simmsar_sigma_stab="ridge", simmsar_sigma_stab_lambda=1e-2),
  floor_1e3 = list(simmsar_sigma_stab="floor", simmsar_sigma_stab_floor=1e-3),
  floor_1e2 = list(simmsar_sigma_stab="floor", simmsar_sigma_stab_floor=1e-2))

cells <- data.frame(
  N   = c(8,   8,   8,   6,   4,   6),
  D   = c(0.75,0.5, 0.75,0.75,0.25,0.5),
  M   = c(4,   4,   3,   3,   2,   2),
  T   = c(200, 400, 400, 200, 800, 800),
  kind= c("hard","hard","hard","hard","healthy","healthy"),
  stringsAsFactors=FALSE)
n_ts <- 4; MaxIter <- 200; eps <- 1e-5; order <- 1
KAPPA_COND_MAX <- 1e6; MAG <- 10; MEV <- 0.05
cc <- cells[ci,]

score_fit <- function(fit, orig_regimes, M) {
  est_A <- lapply(seq_len(M), function(m) fit$theta$A[[m]][[1]])
  est_S <- lapply(seq_len(M), function(m) as.matrix(fit$theta$sigma[[m]]))
  names(est_A) <- paste0("Regime", seq_len(M))
  orig_B <- lapply(seq_len(M), function(m) orig_regimes[[paste0("Regime",m)]][["Beta"]])
  names(orig_B) <- paste0("Regime", seq_len(M))
  am <- match_regimes(orig_B, est_A); if (is.null(am)) return(NULL)
  rows <- vector("list", M)
  for (m in seq_len(M)) {
    oi <- am[[m,1]]; ei <- am[[m,2]]; oreg <- orig_regimes[[paste0("Regime",oi)]]
    oB <- oreg[["Beta"]]; oK <- oreg[["kappa"]]
    eB <- est_A[[ei]]; eBt <- eB; eBt[abs(eBt)<MEV] <- 0
    S  <- est_S[[ei]]; n_exact_zero_S <- sum(S == 0)
    rc <- tryCatch(rcond(S), error=function(e) 0); cond <- if (is.finite(rc)&&rc>0) 1/rc else Inf
    inv <- tryCatch(solve(S), error=function(e) NULL)
    kcor<-ksen<-kspec<-knrmse<-koffmax<-NA_real_; kempty<-NA; ill<-!is.finite(cond)||cond>KAPPA_COND_MAX
    if (!is.null(inv)) {
      K <- symmetrize_matrix(inv); off <- !diag(TRUE,nrow(K)); koffmax <- max(abs(K[off]))
      Kt <- K; Kt[off & abs(K)<MEV] <- 0; kempty <- (sum(Kt[off]!=0)==0)
      kss <- tryCatch(senspec_upper_tri(oK,Kt,diag=FALSE), error=function(e) list(sensitivity=NA,specificity=NA))
      ksen<-kss[["sensitivity"]]; kspec<-kss[["specificity"]]
      kcor <- tryCatch(cor(vectorize_upper_tri(oK,FALSE), vectorize_upper_tri(Kt,FALSE)), error=function(e) NA_real_)
      pos <- vectorize_upper_tri(oK,FALSE)!=0
      knrmse <- if (any(pos)) sqrt(mean((vectorize_upper_tri(oK,FALSE)[pos]-vectorize_upper_tri(Kt,FALSE)[pos])^2))/2 else NA_real_
    }
    bss <- senspec(oB, eBt)
    rows[[m]] <- data.frame(regime=oi,
      Beta_corr=cor(as.vector(oB),as.vector(eBt)), Beta_sen=bss[["sensitivity"]], Beta_spec=bss[["specificity"]],
      Kappa_corr=kcor, Kappa_sen=ksen, Kappa_spec=kspec, NRMSE_Kappa=knrmse,
      cond=cond, kappa_offdiag_max=koffmax, ill_conditioned=ill,
      mag_flagged=is.finite(koffmax)&&koffmax>MAG, kappa_empty=isTRUE(kempty), sigma_exact_zeros=n_exact_zero_S)
  }
  do.call(rbind, rows)
}

set.seed(700 + ci)   # per-cell data seed (independent of other cells / configs)
ts <- generate_timeseries(Density=cc$D, min_edg_val=0.05, max_edg_val=1, M=cc$M, N=cc$N,
                          warmup=50, T=cc$T, totTime=cc$T+50, n_ts=n_ts,
                          remain_lower=0.85, remain_upper=0.85, workers=1)
out <- list(); k <- 0
for (l in 1:n_ts) {
  row <- ts[ts$ts_id==l,]; dmat <- row$timeseries_data[[1]]; d <- ncol(dmat)
  arr <- array(dmat, dim=c(cc$T,1,d)); oreg <- row$regime_dynamics[[1]]
  set.seed(1000*ci + l)
  th0 <- tryCatch(init_theta_msar(arr, M=cc$M, order=order, label="HH"), error=function(e) NULL)
  if (is.null(th0)) next
  do.call(options, c(base_opts, configs[[cfnm]]))
  t0 <- Sys.time()
  fit <- tryCatch(fit_msar(arr, th0, penalty="LASSO", MaxIter=MaxIter, eps=eps, verbose=FALSE),
                  error=function(e) list(.err=conditionMessage(e)),
                  warning=function(w) list(.err=conditionMessage(w)))
  el <- as.numeric(difftime(Sys.time(), t0, units="secs"))
  meta <- data.frame(cell=ci, kind=cc$kind, N=cc$N, D=cc$D, M=cc$M, T=cc$T, ts_id=l, config=cfnm)
  if (!is.null(fit$.err)) {
    k<-k+1; out[[k]] <- cbind(meta, regime=NA, failed=TRUE, iters=NA, secs=el,
      Beta_corr=NA,Beta_sen=NA,Beta_spec=NA,Kappa_corr=NA,Kappa_sen=NA,Kappa_spec=NA,NRMSE_Kappa=NA,
      cond=NA,kappa_offdiag_max=NA,ill_conditioned=NA,mag_flagged=NA,kappa_empty=NA,sigma_exact_zeros=NA,
      err=substr(fit$.err,1,60)); next
  }
  sc <- score_fit(fit, oreg, cc$M); if (is.null(sc)) next
  k<-k+1; out[[k]] <- cbind(meta[rep(1,nrow(sc)),], sc, failed=FALSE, iters=fit$Iter, secs=el, err=NA)
}
if (length(out) > 0) {
  res <- dplyr::bind_rows(out)
  dir.create("output/diagnostics", showWarnings=FALSE, recursive=TRUE)
  write.csv(res, sprintf("output/diagnostics/ab_part_%d_%s.csv", ci, cfnm), row.names=FALSE)
  cat(sprintf("cell %d (%s N%d D%.2f M%d T%d) config %s: %d rows, mean secs %.1f\n",
      ci, cc$kind, cc$N, cc$D, cc$M, cc$T, cfnm, nrow(res), mean(unique(res[,c("ts_id","secs")])$secs)))
} else cat(sprintf("cell %d config %s: NO ROWS\n", ci, cfnm))
