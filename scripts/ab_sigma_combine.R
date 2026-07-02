# Combine per-(cell,config) A/B parts and print the comparison summary.
suppressMessages({ library(dplyr) })
parts <- list.files("output/diagnostics", pattern="^ab_part_.*\\.csv$", full.names=TRUE)
res <- bind_rows(lapply(parts, read.csv, stringsAsFactors=FALSE))
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(res, sprintf("output/diagnostics/ab_sigma_stab_%s.csv", stamp), row.names=FALSE)
configs <- c("baseline","ridge_1e3","ridge_1e2","floor_1e3","floor_1e2")

cat("################ A/B: Sigma stabilization ################\n")
cat(sprintf("parts: %d | rows: %d\n", length(parts), nrow(res)))

fitlevel <- unique(res[,c("cell","kind","ts_id","config","failed","iters","secs")])
cat("\n-- fit failures / iters / runtime by config (fit-level) --\n")
a <- aggregate(cbind(fail_rate=failed, mean_iters=iters, mean_secs=secs) ~ config, data=fitlevel,
               FUN=function(v) mean(v, na.rm=TRUE), na.action=na.pass)
print(a[match(configs,a$config),], row.names=FALSE)
# non-convergence (hit MaxIter=200)
cat("\n-- fits hitting MaxIter=200 (by config) --\n")
print(tapply(fitlevel$iters, fitlevel$config, function(v) sum(v>=200, na.rm=TRUE))[configs])

ok <- res[!res$failed & !is.na(res$config),]
cat("\n-- degeneracy flags (per-regime) & density check by config --\n")
for (cf in configs) { s <- ok[ok$config==cf,]; if(nrow(s)==0) next
  cat(sprintf("%-10s n=%3d | ill %5.2f%% | mag %5.2f%% | empty %5.2f%% | cond med %.2e | max|Koff| med %.2e p90 %.2e | SigmaExactZeros=%d\n",
      cf, nrow(s), 100*mean(s$ill_conditioned,na.rm=TRUE), 100*mean(s$mag_flagged,na.rm=TRUE),
      100*mean(s$kappa_empty,na.rm=TRUE), median(s$cond[is.finite(s$cond)]),
      median(s$kappa_offdiag_max,na.rm=TRUE), quantile(s$kappa_offdiag_max,.9,na.rm=TRUE),
      sum(s$sigma_exact_zeros,na.rm=TRUE))) }

mets <- c("Beta_corr","Beta_sen","Beta_spec","Kappa_corr","Kappa_sen","Kappa_spec","NRMSE_Kappa")
for (kd in c("hard","healthy")) {
  cat("\n-- recovery by config [",kd,"] (mean over regimes) --\n")
  tab <- do.call(rbind, lapply(configs, function(cf){ s <- ok[ok$config==cf & ok$kind==kd,]
    setNames(round(sapply(mets, function(m) mean(s[[m]], na.rm=TRUE)),3), mets) }))
  rownames(tab) <- configs; print(tab)
}

cat("\n-- HARD worst-case severity (cond & |Koff| max) by config --\n")
h <- ok[ok$kind=="hard",]
for (cf in configs) { s <- h[h$config==cf,]; if(nrow(s)==0) next
  cat(sprintf("%-10s cond max %.2e | max|Koff| max %.2e\n", cf,
      max(s$cond[is.finite(s$cond)], na.rm=TRUE), max(s$kappa_offdiag_max, na.rm=TRUE))) }
cat("\nSaved output/diagnostics/ab_sigma_stab_", stamp, ".csv\nDONE\n", sep="")
