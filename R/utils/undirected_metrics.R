#' Utilities for Undirected (Symmetric) Network Evaluation
#'
#' Helpers to evaluate undirected networks without double-counting edges.
#'
#' @keywords internal

vectorize_upper_tri <- function(mat, diag = FALSE) {
  mat[upper.tri(mat, diag = diag)]
}

symmetrize_matrix <- function(mat) {
  (mat + t(mat)) / 2
}

senspec_upper_tri <- function(orig_mat, est_mat, diag = FALSE) {
  orig_vec <- vectorize_upper_tri(orig_mat, diag = diag)
  est_vec <- vectorize_upper_tri(est_mat, diag = diag)

  true_pos <- (orig_vec != 0) & (est_vec != 0)
  true_neg <- (orig_vec == 0) & (est_vec == 0)

  sensitivity <- sum(true_pos) / sum(orig_vec != 0)
  specificity <- sum(true_neg) / sum(orig_vec == 0)

  list(sensitivity = sensitivity, specificity = specificity)
}

mae_upper_tri <- function(orig_mat, est_mat, diag = FALSE) {
  orig_vec <- vectorize_upper_tri(orig_mat, diag = diag)
  est_vec <- vectorize_upper_tri(est_mat, diag = diag)
  mean(abs(orig_vec - est_vec))
}
