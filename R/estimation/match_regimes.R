#' Match Estimated Regimes to True Regimes (Hungarian Algorithm)
#'
#' Solves the regime label-switching problem by matching each estimated regime
#' to a true regime via the globally optimal one-to-one assignment that
#' maximises the correlation of the vectorised lag-1 Beta matrices. This is the
#' single Beta-based mapping that is subsequently applied consistently to the
#' per-regime raw-estimate table and to the regime-sequence / TPM relabelling.
#'
#' The matching is deterministic (no RNG) and intentionally lives on the
#' estimation side: the matched permutation is needed to label both the
#' per-regime table and the decoded regime sequence while the EM fit object is
#' still in memory. It does NOT compute any scored outcome metric -- it only
#' produces the permutation. (All scoring -- Beta_corr, Kappa_corr, sens/spec,
#' NRMSE, AC correlations -- happens later in the analysis pipeline from the
#' stored raw matrices.)
#'
#' \strong{Degenerate (zero-variance) estimated Betas.} A zero-variance
#' estimated Beta makes its matching correlations undefined (\code{cor()} = NA),
#' so it cannot be placed by the assignment. Rather than dropping the WHOLE fit
#' whenever ANY single regime is degenerate (which loses the healthy regimes'
#' network-recovery rows too), the behaviour depends on how many regimes are
#' degenerate:
#' \itemize{
#'   \item \strong{All M degenerate} -> \code{NULL} (fully failed fit, unchanged;
#'     the caller logs \code{zero_var_beta} and drops the whole fit).
#'   \item \strong{Some but not all degenerate} (0 < |D| < M) -> a
#'     \code{partial_regime_match} object: the \code{M - |D|} HEALTHY estimated
#'     regimes are matched optimally to their best true regimes, and the true
#'     regime(s) whose estimated counterpart was degenerate are left UNMATCHED
#'     (\code{est_Reg_No = NA}). The caller salvages the matched regimes for the
#'     per-regime recovery outcomes (RQ1-3) and treats the fit as missing for the
#'     regime-sequence outcome (RQ4), which cannot be defined from a partial map.
#'   \item \strong{None degenerate} (|D| = 0) -> the usual full M x M matching.
#' }
#'
#' @param orig_Betas Named list of true lag-1 Beta matrices, one per regime,
#'   in regime order (e.g. \code{Regime1}, ..., \code{RegimeM}).
#' @param est_Betas Named list of estimated lag-1 Beta matrices, same length
#'   and order as \code{orig_Betas}.
#'
#' @return One of three things:
#'   \describe{
#'     \item{Full match}{A matrix as returned by \code{\link{assign_regimes}}
#'       with columns \code{orig_Reg_No}, \code{est_Reg_No}, \code{Value} (one
#'       row per true regime, ordered 1..M) -- when no estimated Beta is
#'       degenerate.}
#'     \item{Partial match}{An object of class \code{"partial_regime_match"}: a
#'       list with \code{$assign} (the same 3-column matrix, one row per true
#'       regime 1..M, but with \code{est_Reg_No} and \code{Value} = \code{NA}
#'       for any true regime left unmatched), \code{$degenerate_est} (integer
#'       indices of the zero-variance estimated regimes), \code{$unmatched_true}
#'       (integer indices of the true regimes left unassigned), \code{$n_regimes}
#'       (M) and \code{$n_healthy} (M - |D|) -- when 0 < |D| < M.}
#'     \item{\code{NULL}}{when EVERY estimated Beta is degenerate (|D| = M): the
#'       fit is fully unmatchable and dropped by the caller as before.}
#'   }
#'
#' @seealso \code{\link{assign_regimes}} for the underlying Hungarian solver.
#' @export
match_regimes <- function(orig_Betas, est_Betas) {

  vecs_org <- lapply(orig_Betas, as.vector)
  vecs_est <- lapply(est_Betas,  as.vector)

  M <- length(vecs_est)

  # A zero-variance estimated Beta makes cor() undefined (NA). Identify the set D
  # of degenerate estimated regimes (indices into est_Betas).
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
    assign[true_b, 2] <- healthy[a]           # estimated regime index (into est_Betas)
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
