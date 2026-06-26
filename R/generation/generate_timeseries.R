#' Generate MSAR Time Series Data
#'
#' Generates multivariate time series data from Markov-Switching Autoregressive (MSAR)
#' models with network dynamics. Creates random network dynamics for each regime and
#' simulates time series using regime-switching behavior controlled by a transition matrix.
#'
#' @param Density Numeric vector of network edge densities (0 to 1).
#'   Defines the proportion of possible edges that exist in the network.
#' @param N Integer vector of node counts. Number of variables/nodes in the network.
#' @param M Integer vector of regime counts. Number of distinct network states.
#' @param T Integer vector of time steps. Length of each generated time series (excluding warmup).
#' @param n_ts Integer. Number of time series to generate per condition.
#' @param warmup Integer. Number of initial time steps to discard (allows dynamics to stabilize).
#' @param totTime Integer vector. Total time steps including warmup (typically T + warmup).
#' @param mean_rep Numeric. Mean length of regime repetitions in sequence (currently unused).
#' @param sd_rep Numeric. Standard deviation for regime repetition lengths (currently unused).
#' @param min_edg_val Numeric. Minimum absolute edge weight. Edges below this are set to zero.
#' @param max_edg_val Numeric. Maximum absolute edge weight.
#' @param remain_lower Numeric. Lower bound for probability of remaining in same regime (0 to 1).
#' @param remain_upper Numeric. Upper bound for probability of remaining in same regime (0 to 1).
#'
#' @return A tibble (timeseries_data object) with the following columns:
#'   \describe{
#'     \item{timesteps}{Number of time steps in the series (excluding warmup)}
#'     \item{density}{Edge density of the network}
#'     \item{nodes}{Number of nodes in the network}
#'     \item{regimes}{Number of regimes in the model}
#'     \item{ts_id}{Time series ID (1 to n_ts)}
#'     \item{timeseries_data}{List-column: Matrix of time series data (T × N)}
#'     \item{regime_sequence}{List-column: Vector of regime indices for each time step}
#'     \item{regime_dynamics}{List-column: List of network dynamics for each regime}
#'     \item{transmat}{List-column: Transition probability matrix between regimes}
#'   }
#'
#'   Each regime_dynamics list contains:
#'   \itemize{
#'     \item mu: Mean vector
#'     \item Beta: Autoregressive coefficient matrix (temporal network)
#'     \item sigma: Residual covariance matrix
#'     \item kappa: Precision matrix (contemporaneous network)
#'     \item Beta_ac: Average controllability for Beta
#'   }
#'
#' @param workers Integer. Number of parallel worker processes (via
#'   \code{future::multisession}, portable across Windows/Mac/Linux). Defaults to
#'   \code{max(1, parallel::detectCores(logical = FALSE) - 1)}. Pass
#'   \code{workers = 1} for fully sequential execution (debugging / sanity check).
#'
#' @details
#' Each \emph{generation cell} -- one \code{(density, nodes, regimes, set_id)}
#' combination -- is fully independent and is the unit dispatched to a parallel
#' worker. Within a cell the function:
#'
#' 1. **Generates network dynamics** ONCE (temporal + contemporaneous networks,
#'    covariance matrices) via \code{generate_netdyn}, repeating until the minimum
#'    edge-value constraints are satisfied, plus a transition matrix.
#'
#' 2. **Simulates one VAR(1) series per requested \code{T}**, reusing those
#'    dynamics (dynamics do not depend on \code{T}):
#'    \itemize{
#'      \item Initial state sampled from uniform(0, 5)
#'      \item Regime transitions follow a Markov process via transition matrix
#'      \item State evolution: X[t] = mu + Beta(X[t-1] - mu) + epsilon
#'      \item epsilon ~ MVN(0, sigma) of the regime active at step t
#'    }
#'
#' PARALLELISATION / PERFORMANCE (2026-06):
#' \itemize{
#'   \item \strong{Across cells}: \code{future.apply::future_lapply()} over
#'     \code{future::multisession} (portable, unlike \code{multicore} which is
#'     sequential on Windows). \code{future.seed = TRUE} gives one independent
#'     L'Ecuyer-CMRG stream per cell, so results are reproducible run-to-run for a
#'     fixed grid AND invariant to \code{workers} (worker count does not change
#'     the numbers). Cell-isolated streams also make the rejection-sampling retry
#'     count in one cell unable to perturb any other cell.
#'   \item \strong{Cholesky residuals}: the per-regime residual covariance is
#'     Cholesky-factored ONCE per cell (\eqn{R'R = \Sigma}); each timestep then
#'     draws \eqn{\epsilon = R' z,\ z \sim N(0, I)}, distributionally identical to
#'     the old per-step \code{mvtnorm::rmvnorm()} but without re-factorising
#'     \eqn{\Sigma} on every one of the (up to thousands of) steps.
#' }
#' Because both changes alter the RNG-draw mechanism, exact numeric values differ
#' from any pre-2026-06 SERIAL run with the same \code{set.seed()}; the
#' distribution is unchanged and this is a one-time re-baseline.
#'
#' A live cross-process progress bar is shown via \pkg{progressr} when installed.
#'
#' @note
#' Dependencies are loaded centrally via R/dependencies.R. The driver uses
#' \pkg{tibble} (output assembly), \pkg{future} / \pkg{future.apply} (parallelism)
#' and \pkg{progressr} (optional progress bar). The per-worker generation path
#' itself is pure base/stats R (no attached packages required in workers).
#'
#' @seealso
#' \code{\link{generate_netdyn}} for network dynamics generation
#' \code{\link{generate_transmat}} for transition matrix generation
#'
#' @examples
#' \dontrun{
#' # Generate simple 2-regime time series
#' ts_data <- generate_timeseries(
#'   Density = 0.3,
#'   N = 4,
#'   M = 2,
#'   T = 1000,
#'   n_ts = 10,
#'   warmup = 50,
#'   totTime = 1050,
#'   mean_rep = 10,
#'   sd_rep = 2,
#'   min_edg_val = 0.05,
#'   max_edg_val = 1,
#'   remain_lower = 0.33,
#'   remain_upper = 0.66
#' )
#'
#' # Access specific time series
#' ts_data %>% filter(timesteps == 1000, density == 0.3, ts_id == 1)
#' }
#'
#' @export
# Dependencies are loaded centrally via R/dependencies.R
# Driver: tibble, future, future.apply, progressr (optional), parallel.
# Per-worker generation path: base/stats only.

source("R/generation/generate_netdyn.R")
source("R/generation/generate_transmat.R")

# -----------------------------------------------------------------------------
# NOTE on worker packages: unlike estimate_MSAR (EM + LASSO needs the heavy
# stat/LASSO stack), the entire generation path is pure base/stats R --
# generate_netdyn/Beta/kappa/random, check_stability, average_controllability,
# generate_transmat, and the simulation use only matrix(), runif(), rnorm(),
# eigen(), solve(), chol(), sample.int(). The worker functions are auto-exported
# by future's globals detection, so generation workers deliberately load NO
# packages: calling load_packages() here would attach ~30 unused packages
# (lme4/glmnet/plotly/...) per worker and dwarf the (now Cholesky-cheap) compute.

# -----------------------------------------------------------------------------
# Build the list of independent generation cells. One cell per
# (Density, N, M, set_id) combination: it owns its network dynamics and produces
# one simulated series per requested T (dynamics are reused across T). This is
# the unit dispatched to a parallel worker.
build_generation_cells <- function(Density, N, M, n_ts) {
  grid <- expand.grid(
    set_id = seq_len(n_ts), k = seq_along(M), j = seq_along(N), i = seq_along(Density),
    KEEP.OUT.ATTRS = FALSE
  )
  lapply(seq_len(nrow(grid)), function(r) {
    list(
      i = grid$i[r], j = grid$j[r], k = grid$k[r], set_id = grid$set_id[r],
      density_val = Density[grid$i[r]],
      nodes_val   = N[grid$j[r]],
      regimes_val = M[grid$k[r]]
    )
  })
}

# -----------------------------------------------------------------------------
# Simulate ONE VAR(1) regime-switching series of length totTime_t. Sequential in
# time (X[t] depends on X[t-1]); the residual at each step is drawn from the
# covariance of the regime active at that step, using a pre-computed per-regime
# Cholesky factor:  eps = t(R) %*% z  with  R'R = Sigma  =>  Cov(eps) = Sigma.
# Distributionally identical to mvtnorm::rmvnorm(1, 0, Sigma), without
# re-factorising Sigma on every step.
simulate_one_series <- function(dynamics, transmat, chol_sigma, N_j, M_k, totTime_t, warmup) {
  init <- runif(N_j, min = 0, max = 5)
  Rseq <- integer(totTime_t - 1)
  X    <- matrix(init, nrow = totTime_t, ncol = N_j)

  reg_index <- sample.int(M_k, 1, replace = TRUE)
  Rseq[1]   <- reg_index
  for (m in 2:totTime_t) {
    # Regime in next timestep from current regime and transition matrix
    reg_index    <- sample.int(M_k, 1, prob = transmat[reg_index, ])
    Rseq[m - 1]  <- reg_index

    curreg_mu   <- dynamics[[reg_index]][["mu"]]
    curreg_Beta <- dynamics[[reg_index]][["Beta"]]
    eps         <- crossprod(chol_sigma[[reg_index]], rnorm(N_j))  # t(R) %*% z ~ MVN(0, Sigma)

    X[m, ] <- curreg_mu + curreg_Beta %*% (X[m - 1, ] - curreg_mu) + eps
  }

  list(timeseries_data = X[-seq_len(warmup), , drop = FALSE], regime_sequence = Rseq)
}

# -----------------------------------------------------------------------------
# Run ONE generation cell: build dynamics once (rejection sampling + transmat +
# per-regime Cholesky), then simulate one series per T. Returns one assembled
# row-list per T, each tagged with the historical nested-loop order key so the
# final tibble can be restored to the original row order.
generate_one_cell <- function(cell, T, totTime, n_T, n_D, n_N, n_M, n_ts,
                              warmup, min_edg_val, max_edg_val,
                              remain_lower, remain_upper, max_attempts,
                              progress_fun = NULL) {
  N_j <- cell$nodes_val
  M_k <- cell$regimes_val
  D_i <- cell$density_val

  # --- Network dynamics for all regimes (rejection sampling), generated ONCE
  regime_dynamics <- vector("list", M_k)
  for (mm in seq_len(M_k)) {
    attempts <- 0
    repeat {
      attempts <- attempts + 1
      W <- generate_netdyn(N_j, D_i, min_edg_val, max_edg_val)
      kappa_offdiag <- W[["kappa"]][lower.tri(W[["kappa"]], diag = FALSE)]
      if (isTRUE(W[["Beta_stability"]]) &&
          all(abs(W[["Beta"]][W[["Beta"]] != 0]) >= min_edg_val) &&
          all(abs(kappa_offdiag[kappa_offdiag != 0]) >= min_edg_val) &&
          any(kappa_offdiag != 0)) {
        break
      }
      if (attempts >= max_attempts) {
        stop("Failed to generate valid network dynamics after ", max_attempts,
             " attempts. Check Density/min_edg_val/max_edg_val settings.")
      }
    }
    regime_dynamics[[mm]] <- W
  }
  names(regime_dynamics) <- paste0("Regime", seq_len(M_k))

  transmat <- generate_transmat(M_k, remain_lower, remain_upper)

  # Pre-factor each regime's residual covariance ONCE (chol: R'R = Sigma),
  # reused across every timestep and every T-series in this cell.
  chol_sigma <- lapply(regime_dynamics, function(W) chol(W[["sigma"]]))

  # --- One series per requested T, reusing the same dynamics ---------------
  rows <- vector("list", n_T)
  for (t in seq_len(n_T)) {
    sim <- simulate_one_series(regime_dynamics, transmat, chol_sigma,
                               N_j, M_k, totTime[t], warmup)

    # Historical row order: nested loops were t > i(Density) > j(N) > k(M) > set_id
    order_key <- ((((t - 1) * n_D + (cell$i - 1)) * n_N + (cell$j - 1)) * n_M +
                    (cell$k - 1)) * n_ts + cell$set_id

    rows[[t]] <- list(
      order_key       = order_key,
      timesteps       = T[t],
      density         = D_i,
      nodes           = N_j,
      regimes         = M_k,
      ts_id           = cell$set_id,
      timeseries_data = sim$timeseries_data,
      regime_sequence = sim$regime_sequence,
      regime_dynamics = regime_dynamics,
      transmat        = transmat
    )
    if (!is.null(progress_fun)) progress_fun()
  }

  rows
}

generate_timeseries <- function(Density,
                                N,
                                M,
                                T,
                                n_ts,
                                warmup,
                                totTime,
                                mean_rep,
                                sd_rep,
                                min_edg_val,
                                max_edg_val,
                                remain_lower,
                                remain_upper,
                                max_attempts = 1000,
                                workers = NULL
                                ) {

  if (is.null(workers)) {
    workers <- max(1, parallel::detectCores(logical = FALSE) - 1)
  }

  print(sprintf("Generating MSAR timeseries (%d parallel worker%s)",
                workers, if (workers == 1) "" else "s"), quote = FALSE)

  # future::multisession works identically on Windows/Mac/Linux (unlike
  # future::multicore, which silently falls back to sequential on Windows).
  # workers = 1 collapses to fully sequential execution.
  old_plan <- future::plan()
  if (workers <= 1) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  on.exit(future::plan(old_plan), add = TRUE)

  cells <- build_generation_cells(Density, N, M, n_ts)
  n_D <- length(Density); n_N <- length(N); n_M <- length(M); n_T <- length(totTime)
  n_total <- length(cells) * n_T

  run_cells <- function(progress_fun = NULL) {
    future.apply::future_lapply(
      cells,
      generate_one_cell,
      T = T, totTime = totTime, n_T = n_T, n_D = n_D, n_N = n_N, n_M = n_M,
      n_ts = n_ts, warmup = warmup,
      min_edg_val = min_edg_val, max_edg_val = max_edg_val,
      remain_lower = remain_lower, remain_upper = remain_upper,
      max_attempts = max_attempts, progress_fun = progress_fun,
      future.seed = TRUE,
      future.scheduling = Inf   # one cell per dispatch: cost is heterogeneous
                                # across cells (high N/M dynamics + long T), so
                                # fine-grained dispatch balances workers better
                                # than pre-chunking.
    )
  }

  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers(progressr::handler_progress(
      format = "[:bar] :percent :elapsedfull Elapsed, :eta Remaining"
    ))
    results <- progressr::with_progress({
      p <- progressr::progressor(steps = n_total)
      run_cells(progress_fun = p)
    })
  } else {
    message("Package 'progressr' not installed -- running without a live progress bar.\n",
            "Install it (install.packages(\"progressr\")) for progress reporting across workers.")
    results <- run_cells(progress_fun = NULL)
  }

  # ---------------------------------------------------------------------------
  # Flatten the per-cell row-lists into one tibble, restored to historical order.
  # ---------------------------------------------------------------------------
  all_rows <- unlist(results, recursive = FALSE)
  ord      <- order(vapply(all_rows, function(r) r$order_key, numeric(1)))
  all_rows <- all_rows[ord]

  # Column types match the original schema exactly: timesteps/nodes/regimes are
  # numeric (the original pre-allocated them as integer but assigned the double
  # design-grid values T[t]/N[j]/M[k], promoting each column to double), while
  # ts_id stays integer (assigned the integer loop index). Kept identical so this
  # is a drop-in replacement for any downstream code that inspects the schema.
  timeseries_tibble <- tibble::tibble(
    timesteps       = vapply(all_rows, function(r) as.numeric(r$timesteps), numeric(1)),
    density         = vapply(all_rows, function(r) as.numeric(r$density), numeric(1)),
    nodes           = vapply(all_rows, function(r) as.numeric(r$nodes), numeric(1)),
    regimes         = vapply(all_rows, function(r) as.numeric(r$regimes), numeric(1)),
    ts_id           = vapply(all_rows, function(r) as.integer(r$ts_id), integer(1)),
    timeseries_data = lapply(all_rows, `[[`, "timeseries_data"),
    regime_sequence = lapply(all_rows, `[[`, "regime_sequence"),
    regime_dynamics = lapply(all_rows, `[[`, "regime_dynamics"),
    transmat        = lapply(all_rows, `[[`, "transmat")
  )

  # Add S3 class
  class(timeseries_tibble) <- c("timeseries_data", class(timeseries_tibble))

  return(timeseries_tibble)
}


#' Print Method for timeseries_data
#'
#' @param x A timeseries_data object
#' @param ... Additional arguments (unused)
#'
#' @export
print.timeseries_data <- function(x, ...) {
  cat("MSAR Timeseries Data\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Total time series: %d\n", nrow(x)))
  cat(sprintf("Conditions:\n"))
  cat(sprintf("  Timesteps: %s\n", paste(unique(x$timesteps), collapse = ", ")))
  cat(sprintf("  Density: %s\n", paste(unique(x$density), collapse = ", ")))
  cat(sprintf("  Nodes: %s\n", paste(unique(x$nodes), collapse = ", ")))
  cat(sprintf("  Regimes: %s\n", paste(unique(x$regimes), collapse = ", ")))
  cat(sprintf("  Time series per condition: %d\n", max(x$ts_id)))
  cat("════════════════════════════════════════════════════════════════\n")
  cat("\nUse dplyr::filter() to subset by conditions.\n")
  cat("Access data: $timeseries_data, $regime_dynamics, etc.\n")
  NextMethod()
}
