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
#' @param orig_Betas Named list of true lag-1 Beta matrices, one per regime,
#'   in regime order (e.g. \code{Regime1}, ..., \code{RegimeM}).
#' @param est_Betas Named list of estimated lag-1 Beta matrices, same length
#'   and order as \code{orig_Betas}.
#'
#' @return A matrix as returned by \code{\link{assign_regimes}} with columns
#'   \code{orig_Reg_No}, \code{est_Reg_No}, \code{Value} (one row per regime,
#'   ordered by true regime 1..M), or \code{NULL} if any estimated Beta is
#'   degenerate (zero variance), which makes the matching correlations
#'   undefined.
#'
#' @seealso \code{\link{assign_regimes}} for the underlying Hungarian solver.
#' @export
match_regimes <- function(orig_Betas, est_Betas) {

  vecs_org <- lapply(orig_Betas, as.vector)
  vecs_est <- lapply(est_Betas,  as.vector)

  # A zero-variance estimated Beta makes cor() undefined (NA); the fit cannot be
  # matched and is reported as degenerate to the caller.
  if (any(vapply(vecs_est, function(v) sd(v) == 0, logical(1)))) {
    return(NULL)
  }

  cor_results <- outer(
    seq_along(vecs_org),
    seq_along(vecs_est),
    Vectorize(function(a, b) cor(vecs_org[[a]], vecs_est[[b]]))
  )

  assign_regimes(cor_results)
}
