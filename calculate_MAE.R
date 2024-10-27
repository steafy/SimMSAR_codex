# Function to calculate mean absolute error (MAE)
calculate_MAE <- function(A, B) {
 
  # Number of elements from matrix A
  n_elements <- prod(dim(A))
  
  # Absolute error for all nodes between matrices
  abs_errors <- abs(A - B)
  
  # Calculate MAE
  MAE <- sum(abs_errors) / n_elements
  
  return(MAE)
}