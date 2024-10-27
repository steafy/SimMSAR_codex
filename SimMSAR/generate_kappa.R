# Function to generate precision matrix (kappa)

generate_kappa <- function(N, Density, min_edg_val, max_edg_val) {
kappa <- matrix(0, nrow = N, ncol = N)

# Set number of non null elements of kappa
num_elements <- round(Density * ((N * (N - 1)) / 2), 0)

# Generate kappa with random values
indices <- combn(N, 2)
chosen_indices <- indices[, sample(ncol(indices), num_elements, replace = FALSE)]

for (index in 1:num_elements) {
  row <- chosen_indices[1, index]
  col <- chosen_indices[2, index]
  value <- generate_random(min_edg_val, max_edg_val)
  kappa[row, col] <- value
  kappa[col, row] <- value
}

# Make kappa positive definite by ensuring diagonal dominance
for(i in 1:N) {
  kappa[i, i] <- sum(abs(kappa[i, -i])) + 0.1
}

# Check if kappa is positive definite
kappa_eigval <- eigen(kappa)$values
kappa_posdef <- ifelse(all(kappa_eigval >= 0), "Yes", "No")

return(list(kappa = kappa,
            kappa_posdef = kappa_posdef
            ))

}


# 
# generate_kappa <- function(N, Density, min_edg_val, max_edg_val) {
#   kappa <- matrix(0, nrow = N, ncol = N)
#   
#   # Set number of non-null elements of kappa
#   num_elements <- round(Density * ((N * (N - 1)) / 2), 0)
#   
#   # Generate kappa with random values
#   indices <- combn(N, 2)  
#   chosen_indices <- indices[, sample(ncol(indices), num_elements, replace = FALSE)]
#   
#   for (index in 1:num_elements) {
#     row <- chosen_indices[1, index]
#     col <- chosen_indices[2, index]
#     
#     # Generate random value for kappa with min constraint for Wcont
#     diag_sqrt <- sqrt((sum(abs(kappa[row, -row])) + 0.1) * (sum(abs(kappa[col, -col])) + 0.1))
#     min_kappa_value <- 0.1 * diag_sqrt
#     value <- generate_random(max(min_edg_val, min_kappa_value), max_edg_val)
#     
#     kappa[row, col] <- value
#     kappa[col, row] <- value
#   }
#   
#   # Ensure diagonal dominance for positive definiteness
#   for (i in 1:N) {
#     kappa[i, i] <- sum(abs(kappa[i, -i])) + 0.1
#   }
#   
#   # Check if kappa is positive definite
#   kappa_eigval <- eigen(kappa)$values
#   kappa_posdef <- ifelse(all(kappa_eigval >= 0), "Yes", "No")
#   
#   return(list(kappa = kappa,
#               kappa_posdef = kappa_posdef))
# }
