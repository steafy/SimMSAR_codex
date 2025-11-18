#' Generate Regime Transition Matrix
#'
#' Creates a stochastic transition matrix for Markov-switching between regimes.
#' Each row represents transition probabilities from one regime to all others.
#'
#' @param M Integer. Number of regimes.
#' @param remain_lower Numeric (0 to 1). Lower bound for probability of staying in current regime.
#' @param remain_upper Numeric (0 to 1). Upper bound for probability of staying in current regime.
#'
#' @return A transition matrix (M × M) where:
#'   \itemize{
#'     \item Each row sums to 1
#'     \item Diagonal elements (stay probabilities) are between remain_lower and remain_upper
#'     \item Off-diagonal elements (switch probabilities) are randomly distributed
#'   }
#'
#' @details
#' For each regime i, the transition probabilities are generated as:
#'
#' 1. **Stay probability**: Sample p[i,i] ~ Uniform(remain_lower, remain_upper)
#'
#' 2. **Switch probabilities**: Remaining probability (1 - p[i,i]) is divided randomly
#'    among other regimes:
#'    \itemize{
#'      \item Sample M-1 random values
#'      \item Normalize to sum to (1 - p[i,i])
#'      \item Assign to off-diagonal elements
#'    }
#'
#' This ensures each row is a valid probability distribution over next states.
#'
#' @note
#' Higher values of remain_lower and remain_upper create more persistent regimes
#' (longer average regime duration). Lower values create more frequent switching.
#'
#' @seealso
#' \code{\link{generate_timeseries}} which uses this function
#'
#' @examples
#' \dontrun{
#' # Generate transition matrix for 3 regimes with moderate persistence
#' transmat <- generate_transmat(M = 3, remain_lower = 0.33, remain_upper = 0.66)
#' print(rowSums(transmat))  # Should all be 1.0
#' }
#'
#' @export
generate_transmat <- function(M, remain_lower, remain_upper) {
  # Initialize transmat 
  transmat <- matrix(0, nrow = M, ncol = M)
  
  for (i in 1:M) {
    # Set probabilities to remain in current regime
    remain_val <- runif(1, remain_lower, remain_upper)
    transmat[i, i] <- remain_val
    
    # Set probabilities for regime switching in next timestep
    switch_val <- 1 - remain_val
    if (M > 1) {
      rest <- runif(M - 1)
      rest <- rest / sum(rest) * switch_val
      transmat[i, -i] <- rest
    }
  }
  
  return(transmat)
}

