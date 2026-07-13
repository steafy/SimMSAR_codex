# Match estimated temporal networks (est_As) to the true ones (orig_As) by
# vectorised correlation + Hungarian assignment (assign_regimes). Degenerate
# (zero-variance) estimated regimes are handled explicitly: all degenerate ->
# NULL (unmatchable); none -> full M x M match; some -> partial match of the
# healthy regimes, leaving the affected true regime(s) unassigned.
match_regimes <- function(orig_As, est_As) {

  vecs_org <- lapply(orig_As, as.vector)
  vecs_est <- lapply(est_As,  as.vector)

  M <- length(vecs_est)

  # A zero-variance estimated A makes cor() undefined (NA). Identify the set D
  # of degenerate estimated regimes (indices into est_As).
  degenerate <- which(vapply(vecs_est, function(v) sd(v) == 0, logical(1)))

  # |D| == M: every regime degenerate -> fully unmatchable, unchanged behaviour.
  if (length(degenerate) == M) {
    return(NULL)
  }

  # |D| == 0: existing full M x M matching (true regimes as rows).
  if (length(degenerate) == 0) {
    cor_results <- outer(
      seq_along(vecs_org),
      seq_along(vecs_est),
      Vectorize(function(a, b) cor(vecs_org[[a]], vecs_est[[b]]))
    )
    return(assign_regimes(cor_results))
  }

  # 0 < |D| < M: PARTIAL match. Match only the healthy estimated regimes, leaving
  # the true regime(s) whose best estimated counterpart was degenerate unassigned.
  #
  # Orientation matters: clue::solve_LSAP() requires nrow <= ncol, so the
  # M_healthy healthy estimated regimes are the ROWS and all M true regimes the
  # COLUMNS. solve_LSAP then assigns each healthy estimated regime to a DISTINCT
  # true regime; the unassigned true regimes are exactly the ones we want left
  # out (their estimated counterpart is degenerate and carries no usable signal).
  healthy <- setdiff(seq_len(M), degenerate)

  cor_results <- outer(
    seq_along(healthy),          # rows: healthy estimated regimes
    seq_along(vecs_org),         # cols: all true regimes
    Vectorize(function(a, b) cor(vecs_est[[healthy[a]]], vecs_org[[b]]))
  )

  # solve_LSAP minimises cost; convert correlations to costs as assign_regimes does.
  cost_matrix <- max(cor_results) - cor_results
  assignment  <- clue::solve_LSAP(cost_matrix, maximum = FALSE)

  # Build the standard 3-column assign matrix in TRUE-regime order 1..M, with NA
  # in the est_Reg_No / Value columns for any true regime left unmatched.
  assign <- matrix(
    NA_real_, nrow = M, ncol = 3,
    dimnames = list(NULL, c("orig_Reg_No", "est_Reg_No", "Value"))
  )
  assign[, 1] <- seq_len(M)
  for (a in seq_along(healthy)) {
    true_b <- as.integer(assignment[a])       # true regime matched to healthy est a
    assign[true_b, 2] <- healthy[a]           # estimated regime index (into est_As)
    assign[true_b, 3] <- cor_results[a, true_b]
  }
  unmatched_true <- which(is.na(assign[, 2]))

  structure(
    list(
      assign         = assign,
      degenerate_est = as.integer(degenerate),
      unmatched_true = as.integer(unmatched_true),
      n_regimes      = M,
      n_healthy      = length(healthy)
    ),
    class = "partial_regime_match"
  )
}
