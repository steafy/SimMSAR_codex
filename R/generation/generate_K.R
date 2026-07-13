#' Generate Precision Matrix (K)
#'
#' Generates a positive definite precision matrix for the contemporaneous network
#' with specified edge density. Ensures positive definiteness through diagonal dominance.
#'
#' @param N Integer. Number of nodes in the network.
#' @param Density Numeric (0 to 1). Target edge density for off-diagonal elements.
#' @param min_edg_val Numeric. Minimum absolute value for non-zero edges.
#' @param max_edg_val Numeric. Maximum absolute value for non-zero edges.
#'
#' @return A list containing:
#'   \describe{
#'     \item{K}{Precision matrix (N × N), symmetric and positive definite}
#'     \item{K_posdef}{"Yes" or "No" indicating positive definiteness}
#'   }
#'
#' @details
#' The function generates K as follows:
#'
#' 1. **Off-diagonal elements**: Randomly selects edges from lower triangle based on
#'    Density parameter, assigns random values from [min_edg_val, max_edg_val] with
#'    random signs, and mirrors to upper triangle for symmetry
#'
#' 2. **Diagonal dominance**: Sets diagonal elements to sum of absolute row values + 0.1:
#'    \deqn{K[i,i] = \sum_{j \neq i} |K[i,j]| + 0.1}
#'
#' 3. **Verification**: Checks positive definiteness via eigenvalue analysis
#'
#' The precision matrix K is inverted to obtain the covariance matrix sigma.
#' K itself (with its dominant diagonal) is used directly as the
#' contemporaneous network.
#'
#' @note
#' The diagonal dominance approach guarantees positive definiteness, which ensures
#' K can be inverted to obtain a valid covariance matrix.
#'
#' @seealso
#' \code{\link{generate_random}} for random value generation
#' \code{\link{generate_netdyn}} which uses this function
#'
#' @examples
#' \dontrun{
#' # Generate precision matrix for 4 nodes
#' K_result <- generate_K(N = 4, Density = 0.3,
#'                                min_edg_val = 0.05, max_edg_val = 1)
#' print(K_result$K_posdef)  # Should be "Yes"
#' }
#'
#' @export
generate_K <- function(N, Density, min_edg_val, max_edg_val) {
K <- matrix(0, nrow = N, ncol = N)

# Set number of non null elements of K
num_elements <- round(Density * ((N * (N - 1)) / 2), 0)

# Generate K with random values
indices <- combn(N, 2)
chosen_indices <- indices[, sample(ncol(indices), num_elements, replace = FALSE)]

for (index in 1:num_elements) {
  row <- chosen_indices[1, index]
  col <- chosen_indices[2, index]
  value <- generate_random(min_edg_val, max_edg_val)
  K[row, col] <- value
  K[col, row] <- value
}

# Make K positive definite by ensuring diagonal dominance
for(i in 1:N) {
  K[i, i] <- sum(abs(K[i, -i])) + 0.1
}

# Check if K is positive definite
K_eigval <- eigen(K)$values
K_posdef <- ifelse(all(K_eigval >= 0), "Yes", "No")

return(list(K = K,
            K_posdef = K_posdef
            ))

}
