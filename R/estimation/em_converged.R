#' Check EM Algorithm Convergence
#'
#' Determines if the EM algorithm has converged by comparing current and previous
#' log-likelihood values. Also detects decreasing log-likelihood (theoretical violation).
#'
#' @param loglik Numeric. Current iteration log-likelihood value.
#' @param previous_loglik Numeric. Previous iteration log-likelihood value.
#' @param threshold Numeric. Relative convergence threshold. Default: 1e-4.
#' @param verbose Logical. Print warning if log-likelihood decreased. Default: FALSE.
#'
#' @return List with two components:
#'   \describe{
#'     \item{converged}{Integer. 1 if converged, 0 otherwise}
#'     \item{decrease}{Integer. 1 if log-likelihood decreased, 0 otherwise}
#'   }
#'
#' @details
#' **Convergence Criterion**:
#'
#' The algorithm is considered converged when:
#'
#' \deqn{\frac{|loglik - previous\_loglik|}{avg\_loglik} < threshold}
#'
#' where \eqn{avg\_loglik = \frac{|loglik| + |previous\_loglik| + threshold}{2}}
#'
#' This relative criterion scales with the magnitude of the log-likelihood,
#' making it robust across different data scales.
#'
#' **Decrease Detection**:
#'
#' EM algorithm guarantees non-decreasing log-likelihood. If loglik decreases
#' by more than 0.01 (allowing small numerical imprecision):
#' \itemize{
#'   \item Sets decrease = 1
#'   \item Prints warning if verbose = TRUE
#'   \item May indicate numerical issues or bugs in M-step
#' }
#'
#' **Special Cases**:
#' \itemize{
#'   \item If previous_loglik = -Inf (first iteration): No convergence check performed
#'   \item Always returns converged=0 and decrease=0 for first iteration
#' }
#'
#' @note
#' \itemize{
#'   \item Typical threshold: 1e-4 to 1e-6
#'   \item Smaller threshold = stricter convergence = more iterations
#'   \item Log-likelihood decrease suggests numerical instability
#'   \item Used internally by \code{fit_msar}
#' }
#'
#' @seealso
#' \code{\link{fit_msar}} which uses this function to check convergence
#'
#' @examples
#' \dontrun{
#' # Simulate EM iterations
#' ll_history <- c(-1500, -1350, -1320, -1315, -1314.5, -1314.3)
#'
#' for (i in 2:length(ll_history)) {
#'   result <- EM_converged_patched_2(
#'     loglik = ll_history[i],
#'     previous_loglik = ll_history[i-1],
#'     threshold = 1e-4,
#'     verbose = TRUE
#'   )
#'   if (result$converged == 1) {
#'     print(paste("Converged at iteration", i))
#'     break
#'   }
#' }
#' }
#'
#' @export
em_converged <-
function(loglik, previous_loglik, threshold = 1e-4, verbose = FALSE) {

  converged = 0;
  decrease = 0;
  if(!(previous_loglik==-Inf)){
  if (loglik - previous_loglik < -1e-2) # allow for a little imprecision 
    {
      if( verbose ) { 
        print(paste("******likelihood decreased from ",previous_loglik," to ", loglik,sep=""),quote = FALSE)
      }
      decrease = 1;
    }

  delta_loglik = abs(loglik - previous_loglik);
  avg_loglik = (abs(loglik) + abs(previous_loglik) + threshold)/2;
  bb = ((delta_loglik/avg_loglik) < threshold)
  if (bb) {converged = 1}
  }

  res <- NULL
  res$converged <- converged
  res$decrease <- decrease
  return(res)

}
