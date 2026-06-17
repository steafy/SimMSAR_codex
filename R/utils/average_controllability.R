#' Average Controllability (Finite Horizon)
#'
#' Computes node-wise average controllability for a directed system matrix,
#' following the finite-horizon definition used in the theory section of this
#' thesis: the squared Euclidean norm of \code{A^t \%*\% e_i}, summed over a
#' finite set of discrete time steps.
#'
#' @param A Numeric square matrix. System matrix under the convention
#'   \code{x_t = A x_{t-1}}, i.e. \code{A[i, j]} is the effect of node j at
#'   time t-1 on node i at time t. For this project, A is the raw (unstandardized)
#'   temporal Beta matrix; it is used as-is, without transposition.
#' @param T_ac Integer. Number of time steps over which controllability is
#'   accumulated (t = 0, ..., T_ac - 1). Default: 25.
#'
#' @return Numeric vector of length \code{nrow(A)}, one average controllability
#'   value per node.
#'
#' @details
#' Node-wise average controllability of node i is defined as
#' \deqn{AC_i = \sum_{t=0}^{T_{ac}-1} \| A^t e_i \|^2}
#' which equals the squared Euclidean norm of the i-th column of \code{A^t}.
#' \code{A^t} is computed iteratively (\code{P_0 = I}, \code{P_t = A \%*\% P_{t-1}})
#' rather than by repeated exponentiation, for efficiency. No normalization by
#' \code{(1 + |lambda_max|)} is applied; this is the finite, unnormalized variant.
#'
#' @note
#' Validation (not run automatically): for symmetric matrices, this
#' implementation is expected to agree with \code{netcontrol::ave_control_centrality}
#' in node ranking, up to differences in horizon and normalization convention.
#' netcontrol remains a project dependency (see R/dependencies.R) for this kind
#' of manual consistency check, but is no longer called in the pipeline itself.
#'
#' @examples
#' \dontrun{
#' A <- matrix(c(0.3, 0.1, 0.2, 0.4), 2, 2)
#' average_controllability(A, T_ac = 25)
#' }
#'
#' @export
average_controllability <- function(A, T_ac = 25) {
  n <- nrow(A)
  ac <- numeric(n)
  P <- diag(n)

  for (t in 0:(T_ac - 1)) {
    if (t > 0) {
      P <- A %*% P
    }
    ac <- ac + colSums(P^2)
  }

  return(ac)
}
