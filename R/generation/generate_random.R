#' Generate Random Edge Weight with Random Sign
#'
#' Samples a random edge weight from a specified range with equal probability
#' of being positive or negative. Used for generating network edge weights.
#'
#' @param min_edg_val Numeric. Minimum absolute value for the edge weight.
#' @param max_edg_val Numeric. Maximum absolute value for the edge weight.
#'
#' @return Numeric value in either [-max_edg_val, -min_edg_val] or
#'   [min_edg_val, max_edg_val] with equal probability (0.5 each).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Randomly chooses sign (positive or negative) with probability 0.5
#'   \item Samples uniformly from appropriate range based on sign
#'   \item Ensures absolute value is between min_edg_val and max_edg_val
#' }
#'
#' This is useful for generating network edges that can represent both
#' excitatory (positive) and inhibitory (negative) connections.
#'
#' @examples
#' \dontrun{
#' # Generate 10 random edge weights
#' weights <- replicate(10, generate_random(0.1, 1.0))
#' }
#'
#' @export
generate_random <- function(min_edg_val, max_edg_val) {
  if (runif(1) < 0.5) {
    # Sample radom value between -max & -min
    return(runif(1, min = -max_edg_val, max = -min_edg_val))
  } else {
    # Sample radom value between min & max
    return(runif(1, min = min_edg_val, max_edg_val))
  }
}