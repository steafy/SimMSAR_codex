#' Bootstrap Stability Selection for MSAR Network Edges
#'
#' Uses bootstrap resampling to identify stable edges in estimated networks.
#' Edges that appear consistently across bootstrap samples are more likely
#' to be true edges rather than noise.
#'
#' @details
#' This implements stability selection (Meinshausen & Bühlmann, 2010) for
#' MSAR network estimation. Instead of using a fixed threshold (min_edg_val),
#' edges are kept if they are selected in a sufficient proportion of bootstrap
#' samples.
#'
#' Benefits:
#' - Data-driven edge selection
#' - Controls false positive rate
#' - Provides selection probabilities for each edge
#'
#' @references
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
#' Journal of the Royal Statistical Society: Series B, 72(4), 417-473.


#' Block Bootstrap for Time Series
#'
#' Creates a bootstrap sample preserving temporal dependence using
#' non-overlapping blocks.
#'
#' @param data Matrix of time series data (T x N)
#' @param block_size Integer. Size of blocks. Default: 20.
#'
#' @return Bootstrap sample matrix (same dimensions as input)
#'
#' @keywords internal
block_bootstrap <- function(data, block_size = 20) {
  T_len <- nrow(data)
  n_blocks <- ceiling(T_len / block_size)

  # Number of complete blocks we can sample from

  n_available <- floor(T_len / block_size)

  if (n_available < 2) {
    # Not enough data for block bootstrap, use regular bootstrap
    boot_indices <- sample(1:T_len, T_len, replace = TRUE)
  } else {
    # Sample block starting positions
    block_starts <- (sample(1:n_available, n_blocks, replace = TRUE) - 1) * block_size + 1

    # Concatenate blocks
    boot_indices <- unlist(lapply(block_starts, function(s) {
      end_idx <- min(s + block_size - 1, T_len)
      s:end_idx
    }))

    # Trim or pad to original length
    if (length(boot_indices) > T_len) {
      boot_indices <- boot_indices[1:T_len]
    } else if (length(boot_indices) < T_len) {
      # Pad with random samples
      extra <- sample(1:T_len, T_len - length(boot_indices), replace = TRUE)
      boot_indices <- c(boot_indices, extra)
    }
  }

  return(data[boot_indices, , drop = FALSE])
}


#' Bootstrap Stability Selection for Single Time Series
#'
#' Fits MSAR model on multiple bootstrap samples and computes edge
#' selection frequencies.
#'
#' @param data_norm Normalized time series data (T x N matrix)
#' @param M Integer. Number of regimes.
#' @param order Integer. AR order (typically 1).
#' @param MaxIter Integer. Maximum EM iterations per bootstrap.
#' @param n_bootstrap Integer. Number of bootstrap samples. Default: 50.
#' @param block_size Integer. Block size for block bootstrap. Default: 20.
#' @param threshold Numeric. Selection threshold (0-1). Default: 0.6.
#'   Edges selected in > threshold proportion of bootstraps are kept.
#' @param verbose Logical. Print progress. Default: FALSE.
#'
#' @return List with:
#'   \describe{
#'     \item{Beta_selected}{List of logical matrices indicating selected edges per regime (full matrix)}
#'     \item{Kappa_selected}{List of logical matrices indicating selected off-diagonal edges per regime}
#'     \item{Beta_probs}{List of selection probability matrices per regime}
#'     \item{Kappa_probs}{List of selection probability matrices per regime}
#'   }
#'
#' @export
bootstrap_stability_selection <- function(data_norm, M, order, MaxIter,
                                           n_bootstrap = 50, block_size = 20,
                                           threshold = 0.6, verbose = FALSE) {

  N <- ncol(data_norm)
  T_len <- nrow(data_norm)

  # Initialize selection count matrices
  Beta_counts <- list()
  Kappa_counts <- list()
  for (m in 1:M) {
    Beta_counts[[m]] <- matrix(0, nrow = N, ncol = N)
    Kappa_counts[[m]] <- matrix(0, nrow = N, ncol = N)
  }

  successful_boots <- 0

  if (verbose) {
    message("Running ", n_bootstrap, " bootstrap samples for stability selection...")
  }

  for (b in 1:n_bootstrap) {
    if (verbose && b %% 10 == 0) {
      message("  Bootstrap ", b, "/", n_bootstrap)
    }

    # Create bootstrap sample
    boot_data <- block_bootstrap(data_norm, block_size)

    # Make array for MSAR
    boot_array <- array(boot_data, dim = c(T_len, 1, N))

    # Fit MSAR on bootstrap sample
    boot_result <- tryCatch({
      init_and_fit_msar_lasso(
        data = boot_array,
        M = M,
        order = order,
        MaxIter = MaxIter,
        retry = 2,  # Fewer retries for speed
        verbose = FALSE
      )
    }, error = function(e) {
      return(list(fit = NULL))
    })

    if (is.null(boot_result$fit)) {
      next  # Skip failed fits
    }

    successful_boots <- successful_boots + 1

    # Extract estimates and count selections (directly on raw Beta/Kappa)
    boot_Betas <- boot_result$fit$theta$A
    boot_sigmas <- boot_result$fit$theta$sigma

    for (m in 1:M) {
      boot_Beta <- boot_Betas[[m]]$A1
      boot_sigma <- boot_sigmas[[m]]
      boot_kappa <- tryCatch(solve(boot_sigma), error = function(e) NULL)

      if (is.null(boot_kappa)) next

      # Count selections (non-zero edges)
      # Using a small threshold to account for numerical precision
      Beta_counts[[m]] <- Beta_counts[[m]] + (abs(boot_Beta) > 1e-6)
      Kappa_counts[[m]] <- Kappa_counts[[m]] + (abs(boot_kappa) > 1e-6)
    }
  }

  if (successful_boots == 0) {
    warning("No successful bootstrap fits")
    return(NULL)
  }

  if (verbose) {
    message("  Successful bootstraps: ", successful_boots, "/", n_bootstrap)
  }

  # Compute selection probabilities
  Beta_probs <- list()
  Kappa_probs <- list()
  Beta_selected <- list()
  Kappa_selected <- list()

  for (m in 1:M) {
    Beta_probs[[m]] <- Beta_counts[[m]] / successful_boots
    Kappa_probs[[m]] <- Kappa_counts[[m]] / successful_boots

    Beta_selected[[m]] <- Beta_probs[[m]] >= threshold
    Kappa_selected[[m]] <- Kappa_probs[[m]] >= threshold
  }

  return(list(
    Beta_selected = Beta_selected,
    Kappa_selected = Kappa_selected,
    Beta_probs = Beta_probs,
    Kappa_probs = Kappa_probs,
    n_successful = successful_boots
  ))
}


#' Apply Bootstrap Selection to Estimated Networks
#'
#' Zeros out edges that were not selected by bootstrap stability selection.
#' The diagonal of Kappa is never zeroed (it is structurally dominant and
#' required for positive definiteness / invertibility).
#'
#' @param est_Betas List of estimated Beta matrices
#' @param est_Kappas List of estimated Kappa matrices
#' @param stability_result Output from bootstrap_stability_selection
#'
#' @return List with filtered est_Betas and est_Kappas
#'
#' @keywords internal
apply_bootstrap_selection <- function(est_Betas, est_Kappas, stability_result) {

  M <- length(est_Betas)

  for (m in 1:M) {
    # Zero out non-selected Beta edges (full matrix)
    est_Betas[[m]][!stability_result$Beta_selected[[m]]] <- 0

    # Zero out non-selected Kappa edges, but never the diagonal
    kappa_not_selected <- !stability_result$Kappa_selected[[m]]
    diag(kappa_not_selected) <- FALSE
    est_Kappas[[m]][kappa_not_selected] <- 0
  }

  return(list(
    est_Betas = est_Betas,
    est_Kappas = est_Kappas
  ))
}
