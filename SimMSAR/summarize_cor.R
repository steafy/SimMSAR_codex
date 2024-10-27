
# Summarizes vector of correlations
summarize_cor <- function(cor_vec) {
  # Make dataframe from vector
  df_cor_vec <- data.frame(cor_vec = cor_vec)
  
  # Calculate statistics
  stats_cor <- df_cor_vec %>% 
    summarise(N = length(cor_vec),
              Mean = mean(cor_vec),
              Sd = sd(cor_vec),
              Median = median(cor_vec),
              Min = min(cor_vec),
              Max = max(cor_vec)
    )
  
  # Return statistics
  return(stats_cor)
}
