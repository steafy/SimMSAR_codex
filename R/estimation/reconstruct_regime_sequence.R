# Turn smoothed regime probabilities into a hard sequence and relabel it to the
# true regime labels using the assign_regimes() map. Returns list(raw, mapped).
reconstruct_regime_sequence <- function(smoothedprob, assigned_regimes) {
  est_seq_raw <- get_hard_regime_sequence(smoothedprob)
  # col 1 = true label, col 2 = estimated label that maps to it
  label_map <- setNames(assigned_regimes[, 1], assigned_regimes[, 2])
  est_seq_mapped <- as.integer(label_map[as.character(est_seq_raw)])
  list(raw = est_seq_raw, mapped = est_seq_mapped)
}
