#' Initialize and Fit MSAR Model with LASSO and Retry Logic
#'
#' Wrapper function that combines parameter initialization and model fitting for
#' MSAR models with automatic retry on failure. Implements robust estimation by
#' re-initializing and retrying if fitting fails due to numerical issues.
#'
#' @param data 3D array of time series data with dimensions (time, samples, variables).
#'   Typically a single time series, so samples = 1.
#' @param M Integer. Number of regimes (latent states) in the MSAR model.
#' @param order Integer. Autoregressive order (lag). Typically 1 for VAR(1).
#' @param MaxIter Integer. Maximum iterations for EM algorithm. Default: 200.
#' @param retry Integer. Number of retry attempts if fitting fails. Default: 4.
#' @param eps Numeric. Convergence criterion (epsilon). Default: 1e-5.
#' @param verbose Logical. Print detailed progress messages. Default: FALSE.
#'
#' @return List with two components:
#'   \describe{
#'     \item{fit}{Fitted MSAR model object (from \code{fit_msar}), or NULL if all attempts failed}
#'     \item{error}{Error message if fitting failed, or NULL if successful}
#'   }
#'
#' @details
#' The function implements a robust fitting procedure:
#'
#' 1. **Initialization**: Calls \code{init_theta_msar} to generate
#'    initial parameter estimates using hierarchical clustering (HH method)
#'
#' 2. **Fitting**: Calls \code{fit_msar} with LASSO penalty to fit
#'    the MSAR model via EM algorithm
#'
#' 3. **Error Handling**: Wraps fitting in \code{tryCatch} to catch:
#'    \itemize{
#'      \item Numerical errors (singular matrices, convergence failures)
#'      \item Warnings (ill-conditioned problems)
#'    }
#'
#' 4. **Retry Logic**: If fitting fails:
#'    \itemize{
#'      \item Re-initializes parameters (different random initialization)
#'      \item Retries up to \code{retry + 1} total attempts
#'      \item Returns NULL fit if all attempts fail
#'    }
#'
#' This retry mechanism is critical because:
#' \itemize{
#'   \item EM algorithms can fail with poor initialization
#'   \item Numerical issues can occur with ill-conditioned data
#'   \item Different initializations may converge to different local optima
#' }
#'
#' @note
#' \itemize{
#'   \item Each retry gets a fresh random initialization
#'   \item Retry messages are printed to console via \code{message()}
#'   \item Error/warning messages are captured and returned in \code{error} field
#'   \item Successful fits break the retry loop immediately
#' }
#'
#' @seealso
#' \code{\link{init_theta_msar}} for parameter initialization
#' \code{\link{fit_msar}} for EM algorithm fitting
#' \code{\link{estimate_MSAR}} which calls this function
#'
#' @examples
#' \dontrun{
#' # Prepare data array
#' ts_data <- matrix(rnorm(1000 * 4), 1000, 4)
#' data_array <- array(ts_data, dim = c(1000, 1, 4))
#'
#' # Fit 2-regime MSAR model with retry
#' result <- init_and_fit_msar_lasso(
#'   data = data_array,
#'   M = 2,
#'   order = 1,
#'   MaxIter = 200,
#'   retry = 5,
#'   verbose = TRUE
#' )
#'
#' if (!is.null(result$fit)) {
#'   print("Fitting succeeded!")
#' } else {
#'   print(paste("Fitting failed:", result$error))
#' }
#' }
#'
#' @export
init_and_fit_msar_lasso <-
  function(data,
           M,
           order,
           MaxIter = 200,
           retry = 4,
           eps = 1e-5,
           verbose = FALSE) {
    result <- list(fit = NULL, error = NULL)
    
    # browser()
    
    for (try in 1:(retry + 1)) {
      if (try > 1) {
        message("retry: ", try, " of ", retry)
      }
      
      
      # Initialize model for NHMSAR
      model_init <-
        init_theta_msar(
          data = data,
          M = M,
          order = order,
          label = "HH",
          verbose = verbose
        )
      
      # try
      result <- tryCatch({
        # browser()
        
        # Fit MSAR model with NHMSAR
        model_fit <-
          fit_msar(
            data = data,
            theta = model_init,
            penalty = "LASSO",
            MaxIter = MaxIter,
            eps = eps,
            verbose = verbose
          )
        
        list(fit = model_fit, error = NULL)
        
      }, error = function(e) {
        message("fit_msar: can't fit: \n\tError: ", e$message)
        list(fit = NULL, error = e$message)
        # brower()
      }, warning = function(w) {
        message("fit_msar: can't fit: \n\tWarning: ", w$message)
        list(fit = NULL, error = w$message)
        # brower()
      })
      
      fit <- result[["fit"]]
      if (!is.null(fit)) {
        # we have a fit! -> break
        break
      }
      
    } # retry
    
    return(result)
  }