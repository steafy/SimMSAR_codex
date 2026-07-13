# Optional eigenvalue stabilisation of an estimated residual covariance Sigma,
# selected by option simmsar_sigma_stab: "none" (default, no-op), "ridge" (add
# lambda * mean-variance to the diagonal) or "floor" (raise eigenvalues to a
# fraction of the largest). Opt-in; see docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md.
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
