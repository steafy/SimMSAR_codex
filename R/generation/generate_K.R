# Generate the contemporaneous precision matrix K: N nodes, target off-diagonal
# density, made symmetric positive-definite by diagonal dominance
# (K[i,i] = sum_{j!=i} |K[i,j]| + 0.1). The residual covariance is Sigma = K^{-1}.
# Returns a list with K plus validity/diagnostic fields (K_posdef, ...).
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
