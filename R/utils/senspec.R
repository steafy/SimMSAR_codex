#' Calculate Sensitivity and Specificity for Network Recovery
#'
#' Computes edge detection accuracy metrics by comparing original (true) and
#' estimated network matrices. Treats edge detection as a binary classification
#' problem where edges are "positive" and non-edges are "negative".
#'
#' @param orig_mat Matrix. True network (ground truth).
#' @param est_mat Matrix. Estimated network (same dimensions as orig_mat).
#'
#' @return List with two components:
#'   \describe{
#'     \item{sensitivity}{True Positive Rate = TP / (TP + FN).
#'       Proportion of true edges correctly identified.}
#'     \item{specificity}{True Negative Rate = TN / (TN + FP).
#'       Proportion of true non-edges correctly identified.}
#'   }
#'
#' @details
#' The function classifies each matrix element as:
#' \itemize{
#'   \item **True Positive (TP)**: Edge in both original and estimated (orig != 0 & est != 0)
#'   \item **True Negative (TN)**: Non-edge in both (orig == 0 & est == 0)
#'   \item **False Positive (FP)**: Non-edge in original, edge in estimated (orig == 0 & est != 0)
#'   \item **False Negative (FN)**: Edge in original, non-edge in estimated (orig != 0 & est == 0)
#' }
#'
#' Then computes:
#' \itemize{
#'   \item Sensitivity = TP / (TP + FN) = TP / (all true edges)
#'   \item Specificity = TN / (TN + FP) = TN / (all true non-edges)
#' }
#'
#' @note
#' \itemize{
#'   \item High sensitivity means few missed edges (low false negative rate)
#'   \item High specificity means few spurious edges (low false positive rate)
#'   \item Both metrics range from 0 to 1, with 1 being perfect
#'   \item For sparse networks, specificity is typically very high due to many true non-edges
#' }
#'
#' @examples
#' \dontrun{
#' # True 3x3 network
#' true_net <- matrix(c(0, 0.5, 0, 0.3, 0, 0.7, 0, 0, 0), 3, 3)
#'
#' # Estimated network (one false positive, one false negative)
#' est_net <- matrix(c(0, 0.4, 0.2, 0, 0, 0.6, 0, 0, 0), 3, 3)
#'
#' metrics <- senspec(true_net, est_net)
#' print(metrics$sensitivity)  # 2/3 = 0.67 (missed one edge)
#' print(metrics$specificity)  # 5/6 = 0.83 (one false positive)
#' }
#'
#' @export
senspec <- function(orig_mat, est_mat) {
  true_pos <- (orig_mat != 0) & (est_mat != 0)
  true_neg <- (orig_mat == 0) & (est_mat == 0)
  
  n_pos <- sum(orig_mat != 0)
  n_neg <- sum(orig_mat == 0)
  sensitivity <- if (n_pos > 0) sum(true_pos) / n_pos else NA_real_
  specificity <- if (n_neg > 0) sum(true_neg) / n_neg else NA_real_

  return(list("sensitivity" = sensitivity,
              "specificity" = specificity
              )
         ) 
}

