#' Extract Hard Regime Sequence from Smoothed Probabilities
#'
#' Converts smoothed regime probabilities (from the EM forward-backward step)
#' into a hard (argmax) regime sequence, one regime label per time point.
#'
#' @param smoothedprob Matrix of smoothed probabilities as returned by
#'   \code{fit_msar} (\code{FB$probS}) for the M > 1 case used throughout this
#'   pipeline (one time series per fit, i.e. N.samples = 1). For N.samples = 1,
#'   \code{fit_msar} drops the (size-1) sample dimension when reordering
#'   regimes by ascending residual variance, so \code{FB$probS} is a plain
#'   (T - 1) x M matrix (one row fewer than the input series, since the first
#'   time point of an order-1 fit only serves as the initial AR lag and is not
#'   assigned a regime probability).
#'
#' @return Integer vector of length T - 1 with the most probable regime
#'   (1 to M) at each time point.
#'
#' @export
get_hard_regime_sequence <- function(smoothedprob) {
  apply(smoothedprob, 1, which.max)
}


#' Cohen's Kappa (Manual, Unweighted)
#'
#' Computes the chance-corrected agreement (Cohen's kappa) between two
#' categorical sequences of equal length. Implemented manually to avoid
#' adding a new package dependency.
#'
#' @param true_seq Integer/factor vector. Reference (true) sequence.
#' @param est_seq Integer/factor vector. Comparison (estimated) sequence,
#'   same length as \code{true_seq}.
#'
#' @return Numeric. Cohen's kappa, or \code{NA} if the expected agreement
#'   under chance is 1 (degenerate case).
#'
#' @details
#' \deqn{\kappa = \frac{p_o - p_e}{1 - p_e}}
#' where \eqn{p_o} is the observed proportion of agreement and \eqn{p_e} is
#' the proportion of agreement expected by chance, based on the marginal
#' distributions of both sequences.
#'
#' @export
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
