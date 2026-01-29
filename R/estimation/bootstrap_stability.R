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
#'     \item{Wtemp_selected}{List of logical matrices indicating selected edges per regime}
#'     \item{Wcont_selected}{List of logical matrices indicating selected edges per regime}
#'     \item{Wtemp_probs}{List of selection probability matrices per regime}
#'     \item{Wcont_probs}{List of selection probability matrices per regime}
#'   }
#'
#' @export
bootstrap_stability_selection <- function(data_norm, M, order, MaxIter,
                                           n_bootstrap = 50, block_size = 20,
                                           threshold = 0.6, verbose = FALSE) {

  N <- ncol(data_norm)
  T_len <- nrow(data_norm)

  # Initialize selection count matrices
  Wtemp_counts <- list()
  Wcont_counts <- list()
  for (m in 1:M) {
    Wtemp_counts[[m]] <- matrix(0, nrow = N, ncol = N)
    Wcont_counts[[m]] <- matrix(0, nrow = N, ncol = N)
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

    # Extract estimates and count selections
    boot_Betas <- boot_result$fit$theta$A
    boot_sigmas <- boot_result$fit$theta$sigma

    for (m in 1:M) {
      # Get Beta and compute Wtemp
      boot_Beta <- boot_Betas[[m]]$A1
      boot_sigma <- boot_sigmas[[m]]
      boot_kappa <- tryCatch(solve(boot_sigma), error = function(e) NULL)

      if (is.null(boot_kappa)) next

      # Compute Wtemp
      boot_Wtemp <- matrix(0, nrow = N, ncol = N)
      for (n in 1:N) {
        for (o in 1:N) {
          denom <- sqrt(boot_sigma[n, n] * boot_kappa[o, o] + boot_Beta[n, o]^2)
          if (denom > 1e-10) {
            boot_Wtemp[n, o] <- boot_Beta[n, o] / denom
          }
        }
      }

      # Compute Wcont
      boot_Wcont <- matrix(0, nrow = N, ncol = N)
      for (n in 1:N) {
        for (o in 1:N) {
          if (n != o) {
            denom <- sqrt(boot_kappa[n, n] * boot_kappa[o, o])
            if (denom > 1e-10) {
              boot_Wcont[n, o] <- -boot_kappa[n, o] / denom
            }
          }
        }
      }

      # Count selections (non-zero edges)
      # Using a small threshold to account for numerical precision
      Wtemp_counts[[m]] <- Wtemp_counts[[m]] + (abs(boot_Wtemp) > 1e-6)
      Wcont_counts[[m]] <- Wcont_counts[[m]] + (abs(boot_Wcont) > 1e-6)
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
  Wtemp_probs <- list()
  Wcont_probs <- list()
  Wtemp_selected <- list()
  Wcont_selected <- list()

  for (m in 1:M) {
    Wtemp_probs[[m]] <- Wtemp_counts[[m]] / successful_boots
    Wcont_probs[[m]] <- Wcont_counts[[m]] / successful_boots

    Wtemp_selected[[m]] <- Wtemp_probs[[m]] >= threshold
    Wcont_selected[[m]] <- Wcont_probs[[m]] >= threshold
  }

  return(list(
    Wtemp_selected = Wtemp_selected,
    Wcont_selected = Wcont_selected,
    Wtemp_probs = Wtemp_probs,
    Wcont_probs = Wcont_probs,
    n_successful = successful_boots
  ))
}


#' Apply Bootstrap Selection to Estimated Networks
#'
#' Zeros out edges that were not selected by bootstrap stability selection.
#'
#' @param est_Wtemps List of estimated Wtemp matrices
#' @param est_Wconts List of estimated Wcont matrices
#' @param stability_result Output from bootstrap_stability_selection
#'
#' @return List with filtered est_Wtemps and est_Wconts
#'
#' @keywords internal
apply_bootstrap_selection <- function(est_Wtemps, est_Wconts, stability_result) {

  M <- length(est_Wtemps)

  for (m in 1:M) {
    # Zero out non-selected Wtemp edges
    est_Wtemps[[m]][!stability_result$Wtemp_selected[[m]]] <- 0

    # Zero out non-selected Wcont edges
    est_Wconts[[m]][!stability_result$Wcont_selected[[m]]] <- 0
  }

  return(list(
    est_Wtemps = est_Wtemps,
    est_Wconts = est_Wconts
  ))
}
