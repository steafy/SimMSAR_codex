#' Reconstruct the Estimated Regime Sequence from Smoothed Probabilities
#'
#' Hard-decodes the smoothed regime probabilities (xi_{t|T}) of a fitted MSAR
#' model into an integer regime sequence and relabels it with the Beta-based
#' matching permutation, so that sequence index \code{m} refers to the same
#' regime as \code{regime_id == m} in the per-regime raw-estimate table.
#'
#' This lives on the estimation side because the smoothed probabilities are a
#' transient \code{(T - 1) x M} array inside the fit object that we do not want
#' to serialise for thousands of fits just to defer the decoding. Decoding is
#' deterministic; only the resulting integer sequence (raw and relabelled) is
#' stored.
#'
#' @param smoothedprob The \code{(T - 1) x M} matrix of smoothed probabilities
#'   from the fit (\code{model_fit[["smoothedprob"]]}).
#' @param assigned_regimes The matching matrix from \code{\link{match_regimes}}
#'   (columns \code{orig_Reg_No} = true label, \code{est_Reg_No} = estimated
#'   label that maps to it).
#'
#' @return A list with:
#'   \describe{
#'     \item{raw}{Integer vector: the hard (argmax) estimated regime labels in
#'       the model's own (estimated) labelling.}
#'     \item{mapped}{Integer vector: the same sequence relabelled into the true
#'       regime labelling via the matching permutation.}
#'   }
#'
#' @seealso \code{\link{get_hard_regime_sequence}}, \code{\link{match_regimes}}
#' @export
reconstruct_regime_sequence <- function(smoothedprob, assigned_regimes) {
  est_seq_raw <- get_hard_regime_sequence(smoothedprob)
  # col 1 = true label, col 2 = estimated label that maps to it
  label_map <- setNames(assigned_regimes[, 1], assigned_regimes[, 2])
  est_seq_mapped <- as.integer(label_map[as.character(est_seq_raw)])
  list(raw = est_seq_raw, mapped = est_seq_mapped)
}
