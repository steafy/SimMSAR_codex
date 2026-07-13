#' Cheap in-EM stabilization of the residual covariance Sigma
#'
#' Optional, default-OFF regularization of the regime-weighted residual
#' covariance produced by the M-step, to prevent it from becoming pathologically
#' ill-conditioned when a regime is assigned very few effective observations
#' (small \code{postmix[j]}) relative to its \eqn{d(d+1)/2} covariance parameters.
#' See \code{docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md}: the degeneracy is an
#' under-determination problem (near-zero smallest eigenvalue), which is exactly
#' what a ridge lift or an eigenvalue floor targets.
#'
#' It deliberately does NOT sparsify Sigma. Both methods keep the matrix fully
#' DENSE (no exact zeros introduced), so the frozen-support logic in
#' \code{mstep_hh_reduct_msar} (\code{w = which(abs(theta$sigma[[j]]) > 0)})
#' remains the no-op it currently is. This is asserted, not assumed, by
#' \code{tests}/the A/B harness.
#'
#' Controlled entirely by options so call sites don't change and it can be
#' A/B-tested per parallel worker (options are re-applied inside each worker via
#' the same \code{lasso_control} mechanism as the LASSO knobs):
#' \itemize{
#'   \item \code{simmsar_sigma_stab} = \code{"none"} (default) | \code{"ridge"} |
#'     \code{"floor"}. \code{"none"} returns Sigma untouched (production default).
#'   \item \code{simmsar_sigma_stab_lambda} (ridge only, default \code{1e-3}):
#'     \eqn{\Sigma \gets \Sigma + \lambda\,\bar\sigma\, I} with
#'     \eqn{\bar\sigma = \mathrm{tr}(\Sigma)/d} (scale-free: the lift is a fraction
#'     of the matrix's own average variance). Lifts EVERY eigenvalue by the same
#'     amount, so it cannot create a zero and cannot lower the smallest eigenvalue.
#'   \item \code{simmsar_sigma_stab_floor} (floor only, default \code{1e-3}):
#'     eigen-decompose, clip any eigenvalue below \code{floor * max(eig)} up to that
#'     floor, reconstruct. Because the floor is a fraction of the LARGEST
#'     eigenvalue, this bounds the condition number at exactly
#'     \code{1 / simmsar_sigma_stab_floor} (e.g. 1e-3 -> cond <= 1000) while leaving
#'     any Sigma already better-conditioned than that completely untouched.
#' }
#'
#' Both act on EVERY M-step (the ill-conditioning is present from iteration 1 and
#' persists; it also feeds A's re-estimation via \code{S.th} in the reduct
#' step), never only at convergence.
#'
#' @param S A \eqn{d \times d} covariance matrix (one regime).
#' @return The stabilized covariance (same dimensions), symmetric and dense.
#' @keywords internal
#' @export
stabilize_sigma <- function(S) {
  method <- getOption("simmsar_sigma_stab", "none")
  if (is.null(method) || identical(method, "none")) return(S)

  S <- as.matrix(S)
  d <- nrow(S)
  if (is.null(d) || is.na(d) || d < 1) return(S)

  if (identical(method, "ridge")) {
    lam <- getOption("simmsar_sigma_stab_lambda", 1e-3)
    avg_var <- sum(diag(S)) / d
    if (!is.finite(avg_var) || avg_var <= 0) avg_var <- 1
    return(S + diag(lam * avg_var, d))
  }

  if (identical(method, "floor")) {
    fr <- getOption("simmsar_sigma_stab_floor", 1e-3)
    ev <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
    if (is.null(ev)) return(S)          # non-finite Sigma: leave to existing guards
    max_ev <- max(ev$values)
    if (!is.finite(max_ev) || max_ev <= 0) return(S)
    floor_val <- fr * max_ev
    vals <- pmax(ev$values, floor_val)
    Sout <- ev$vectors %*% (vals * t(ev$vectors))   # V diag(vals) V^T
    return((Sout + t(Sout)) / 2)                    # enforce exact symmetry
  }

  S
}
