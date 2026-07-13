# Draw one random edge weight with magnitude in [min_edg_val, max_edg_val] and a
# random sign (equal probability positive/negative).
generate_random <- function(min_edg_val, max_edg_val) {
  if (runif(1) < 0.5) {
    # Sample radom value between -max & -min
    return(runif(1, min = -max_edg_val, max = -min_edg_val))
  } else {
    # Sample radom value between min & max
    return(runif(1, min = min_edg_val, max_edg_val))
  }
}