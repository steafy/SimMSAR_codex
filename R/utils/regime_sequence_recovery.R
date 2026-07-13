# Hard (argmax) regime sequence from a matrix of smoothed regime probabilities:
# one regime label per time point.
get_hard_regime_sequence <- function(smoothedprob) {
  apply(smoothedprob, 1, which.max)
}


# Cohen's kappa (unweighted) between two equal-length categorical sequences,
# implemented manually to avoid a package dependency; NA when chance agreement is
# 1. This is Cohen's kappa for regime-sequence recovery (RQ4) -- unrelated to the
# K precision network.
cohens_kappa_manual <- function(true_seq, est_seq) {
  levels_all <- sort(unique(c(true_seq, est_seq)))
  tab <- table(factor(true_seq, levels = levels_all),
               factor(est_seq, levels = levels_all))
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum((rowSums(tab) / n) * (colSums(tab) / n))
  if (isTRUE(all.equal(pe, 1))) return(NA_real_)
  (po - pe) / (1 - pe)
}
