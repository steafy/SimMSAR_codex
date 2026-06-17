#' Summarize Vector of Correlation Values
#'
#' Computes descriptive statistics for a vector of correlation coefficients or
#' other numeric metrics. Provides comprehensive summary including central tendency,
#' dispersion, and range.
#'
#' @param cor_vec Numeric vector of correlation values (or other metrics like MAE,
#'   sensitivity, etc.).
#'
#' @return Data frame (1 row) with summary statistics:
#'   \describe{
#'     \item{N}{Sample size (length of vector)}
#'     \item{Mean}{Arithmetic mean}
#'     \item{Sd}{Standard deviation}
#'     \item{Median}{50th percentile}
#'     \item{Min}{Minimum value}
#'     \item{Max}{Maximum value}
#'   }
#'
#' @details
#' This function is used throughout the estimation pipeline to aggregate
#' performance metrics across multiple time series or regimes. For example:
#' \itemize{
#'   \item Summarizing Beta correlations across all time series in a condition
#'   \item Aggregating MAE values for all estimated regimes
#'   \item Computing average sensitivity/specificity
#' }
#'
#' The output format (data frame) allows easy combination with other summaries
#' and integration into larger results structures.
#'
#' @note
#' NA values are not explicitly handled. If cor_vec contains NAs, statistics
#' will also be NA. Use \code{na.omit(cor_vec)} before calling if needed.
#'
#' @seealso
#' \code{\link{estimate_MSAR}} which uses this function to aggregate results
#'
#' @examples
#' \dontrun{
#' # Summarize network recovery correlations
#' correlations <- c(0.85, 0.72, 0.91, 0.68, 0.88)
#' summary <- summarize_cor(correlations)
#' print(summary$Mean)  # 0.808
#' print(summary$Sd)    # ~0.095
#' }
#'
#' @export
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
