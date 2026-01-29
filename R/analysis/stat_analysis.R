# =============================================================================
# Statistical Analysis Script for MSAR Models
# =============================================================================
# This is a standalone analysis script that requires all packages to be loaded.
# It should be run AFTER the main MSAR estimation pipeline.

# Load dependencies
if (!exists("ALL_PACKAGES") || !all(c("dplyr", "ggplot2") %in% loadedNamespaces())) {
  source("R/dependencies.R")
}

### Data preparation
## The data is now already in tibble format from estimate_MSAR()
## No need for extraction loops - just add grouping variables

# Verify we have an msar_results object
if (!inherits(MSAR_dynamics_list, "msar_results")) {
  stop("MSAR_dynamics_list is not an msar_results object. Please re-run estimate_MSAR().")
}

# Create corr_results with additional grouping variables
corr_results <- MSAR_dynamics_list %>%
  rename(
    Timesteps = timesteps,
    Density = density,
    Nodes = nodes,
    Regimes = regimes
  ) %>%
  # Add SimID and RegimeIndex
  group_by(Timesteps, Density, Nodes, Regimes) %>%
  mutate(
    SimID = ts_id,
    RegimeIndex = regime_id,
    N = n_distinct(SimID)
  ) %>%
  ungroup() %>%
  mutate(Condition = interaction(Timesteps, Density, Nodes, Regimes, drop = TRUE)) %>%
  select(SimID, Timesteps, Density, Nodes, Regimes, RegimeIndex, N, everything())

# Transform to factors
corr_results <- corr_results %>%
  mutate(
    Timesteps = factor(Timesteps),
    Density = factor(Density),
    Nodes = factor(Nodes),
    Regimes = factor(Regimes)
  )




## Extract stats using the new get_stats() function
descript_stats_summary <- get_stats(MSAR_dynamics_list)

# Convert to the old format for compatibility with existing code
# (Creates separate data frames for each metric)
descript_stats <- list()

metric_mapping <- list(
  "Wtemp_corr" = "Wtemp_corr",
  "Wtemp_MAE" = "Wtemp_MAE",
  "Wtemp_sensitivity" = "Wtemp_sen",
  "Wtemp_specificity" = "Wtemp_spec",
  "Wtemp_ac_corr" = "Wtemp_ac_corr",
  "Wcont_corr" = "Wcont_corr",
  "Wcont_MAE" = "Wcont_MAE",
  "Wcont_sensitivity" = "Wcont_sen",
  "Wcont_specificity" = "Wcont_spec",
  "Wcont_ac_corr" = "Wcont_ac_corr"
)

for (old_name in names(metric_mapping)) {
  new_name <- metric_mapping[[old_name]]

  descript_stats[[old_name]] <- descript_stats_summary %>%
    select(
      Timesteps = timesteps,
      Density = density,
      Nodes = nodes,
      Regimes = regimes,
      Mean = !!sym(paste0(new_name, "_mean")),
      Sd = !!sym(paste0(new_name, "_sd")),
      Median = !!sym(paste0(new_name, "_median")),
      Min = !!sym(paste0(new_name, "_min")),
      Max = !!sym(paste0(new_name, "_max"))
    )
}


####################################################################
### Calculate statistical information for estimation process

## Determine no. of omissions (ts without converging model) (30 - N)
aggr_factors <- corr_results %>%
  dplyr::select(Timesteps, Density, Nodes, Regimes, N) %>%
  dplyr::distinct()

# Calculate expected total (assumes n_ts from the data)
n_ts_per_condition <- max(corr_results$SimID)
n_conditions <- nrow(aggr_factors)
expected_total <- n_ts_per_condition * n_conditions * unique(as.numeric(as.character(corr_results$Regimes)))[1]

omissions <- expected_total - sum(aggr_factors$N)
omissions_per <- (omissions / expected_total) * 100


## Calculate mean no. of estimated models (N) for factorlevels
## Calculate dunn-test to compare N across factorlevels
n_means <- list()
dunn_results <- list()

# First check: Is there any variance in N at all?
if (length(unique(aggr_factors$N)) == 1) {
  message(sprintf("Note: All conditions have the same number of successful fits (N = %d).",
                  unique(aggr_factors$N)))
  message("Kruskal-Wallis and Dunn tests are not applicable - no variance to test.")

  # Still calculate means
  for (i in 1:4) {
    var <- colnames(aggr_factors)[i]
    val <- unique(aggr_factors[[var]])
    for (j in seq_along(val)) {
      x <- val[j]
      subset <- aggr_factors %>% filter(.data[[var]] == x)
      n_means[[paste(x, var)]] <- mean(subset$N)
    }
    dunn_results[[var]] <- list(
      Kruskal = "Not applicable - no variance in N",
      Dunn = "Not applicable - no variance in N"
    )
  }
} else {
  # There is variance - proceed with tests
  for (i in 1:4) {
    var <- colnames(aggr_factors)[i]
    val <- unique(aggr_factors[[var]])

    # Calculate means for each level
    for (j in seq_along(val)) {
      x <- val[j]
      subset <- aggr_factors %>% filter(.data[[var]] == x)
      n_means[[paste(x, var)]] <- mean(subset$N)
    }

    # Only perform statistical tests if there are multiple groups to compare
    if (length(val) > 1) {
      # Check if there's variance within this specific factor
      group_means <- tapply(aggr_factors$N, aggr_factors[[var]], mean)
      if (length(unique(group_means)) == 1) {
        message(sprintf("Skipping tests for '%s' - all groups have same mean N (%.1f)",
                        var, group_means[1]))
        dunn_results[[var]] <- list(
          Kruskal = "Not applicable - no between-group variance",
          Dunn = "Not applicable - no between-group variance"
        )
        next
      }

      # Run tests with error handling
      kw_result <- tryCatch({
        kruskal.test(aggr_factors$N, aggr_factors[[var]])
      }, error = function(e) {
        message(sprintf("Error in Kruskal-Wallis for '%s': %s", var, e$message))
        return(NULL)
      }, warning = function(w) {
        message(sprintf("Warning in Kruskal-Wallis for '%s': %s", var, w$message))
        return(NULL)
      })

      dunn_result <- tryCatch({
        dunn.test(aggr_factors$N, aggr_factors[[var]], method = "bonferroni")
      }, error = function(e) {
        message(sprintf("Error in Dunn test for '%s': %s", var, e$message))
        return(NULL)
      }, warning = function(w) {
        message(sprintf("Warning in Dunn test for '%s': %s", var, w$message))
        return(NULL)
      })

      # Store results (even if NULL, for debugging)
      if (!is.null(kw_result) && !is.null(dunn_result)) {
        dunn_matrix <- data.frame(
          comp = dunn_result$comparisons,
          Z_val = dunn_result$Z,
          p_val = dunn_result$P,
          p.adj = round(dunn_result$P.adjusted, 4)
        )
        dunn_matrix <- dunn_matrix[order(dunn_matrix$p.adj), ]
        dunn_results[[var]] <- list(
          Kruskal = kw_result,
          Dunn = dunn_matrix
        )
      } else {
        dunn_results[[var]] <- list(
          Kruskal = if (is.null(kw_result)) "Test failed" else kw_result,
          Dunn = if (is.null(dunn_result)) "Test failed" else dunn_result
        )
      }
    } else {
      # Only one group - tests not applicable
      message(sprintf("Skipping Kruskal-Wallis and Dunn tests for '%s' - only one group (%s)",
                      var, val))
      dunn_results[[var]] <- list(
        Kruskal = "Not applicable - only one group",
        Dunn = "Not applicable - only one group"
      )
    }
  }
}


###############################################
### Calculate descriptive statistics

desc <- as.data.frame(matrix(nrow = 5, ncol = 10,
                             dimnames = list(c("Mean", "Sd", "Median", "Min", "Max"),
                                             c("Wtemp_corr",
                                               "Wtemp_MAE",
                                               "Wtemp_sensitivity",
                                               "Wtemp_specificity",
                                               "Wtemp_ac_corr",
                                               "Wcont_corr",
                                               "Wcont_MAE",
                                               "Wcont_sensitivity",
                                               "Wcont_specificity",
                                               "Wcont_ac_corr")
                                             )
                             )
                      )

for (i in 1:10) {
  x <- descript_stats[[i]]
  val <- c(mean(x[["Mean"]]),
           sd(x[["Mean"]]),
           mean(x[["Median"]]),
           min(x[["Min"]]),
           max(x[["Max"]])
           )
  desc[ ,i] <- val
}

###################################
# ### Calculate a linear mixed model
# 
# # Check for factors with only one level (will cause issues with scaling and modeling)
# lmm_factor_levels <- sapply(c("Timesteps", "Density", "Nodes", "Regimes"), function(var) {
#   length(unique(corr_results[[var]]))
# })
# 
# if (any(lmm_factor_levels == 1)) {
#   single_level_vars <- names(lmm_factor_levels)[lmm_factor_levels == 1]
#   message(sprintf("Warning: Cannot fit linear mixed models - the following variables have only one level: %s",
#                   paste(single_level_vars, collapse = ", ")))
#   message("Skipping linear mixed model analysis. Run with multiple factor levels to enable this analysis.")
#   linear_mixed_models <- list()
# } else {
#   # Scale factors
#   corr_results <- corr_results %>%
#     mutate(
#       Timesteps_scaled = scale(as.numeric(as.character(Timesteps))),
#       Density_scaled   = scale(as.numeric(as.character(Density))),
#       Nodes_scaled     = scale(as.numeric(as.character(Nodes))),
#       Regimes_scaled   = scale(as.numeric(as.character(Regimes)))
#     )
#   # Set variables to model
#   cols <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")
# 
#   linear_mixed_models <- list()
# 
#   for (col_name in cols) {
# 
#     # Specify model
#     fmla <- as.formula(paste0(col_name, " ~ Timesteps_scaled * Density_scaled * Nodes_scaled * Regimes_scaled +
#                             RegimeIndex +
#                             (1 | Condition/SimID)")
#                      )
# 
#     model <- tryCatch({
#       lmer(fmla, data = corr_results)
#     }, error = function(e) {
#       message(sprintf("Error fitting LMM for %s: %s", col_name, e$message))
#       return(NULL)
#     })
# 
#     if (is.null(model)) {
#       linear_mixed_models[[col_name]] <- "Model failed - see error message above"
#       next
#     }
# 
#     # Summarize model
#     lmm <- as.data.frame(coef(summary(model))) %>%
#       dplyr::mutate(Signif = ifelse(`Pr(>|t|)` < 0.001, "***",
#                                     ifelse(`Pr(>|t|)` < 0.01, "**",
#                                            ifelse(`Pr(>|t|)` < 0.05, "*", ""))))
# 
#     rownames(lmm) <- gsub(":", " × ", gsub("_scaled", "", rownames(lmm)))
# 
#     linear_mixed_models[[col_name]] <- lmm
#   }
# }
# 
# # Export results for each model in the list
# # Only export if we have valid model results (data frames, not error messages)
# if (length(linear_mixed_models) > 0) {
#   for (name in names(linear_mixed_models)) {
#     tab <- linear_mixed_models[[name]]
# 
#     # Skip if this is an error message rather than a model
#     if (is.character(tab)) {
#       message(sprintf("Skipping export for %s: %s", name, tab))
#       next
#     }
# 
#     # Export to HTML
#     tryCatch({
#       kable(tab, caption = paste("Ergebnisse für", name), format = "html") %>%
#         kable_styling(bootstrap_options = c("striped", "hover")) %>%
#         save_kable(file = paste0("Ergebnisse_", name, ".html"))
#     }, error = function(e) {
#       message(sprintf("Error exporting HTML for %s: %s", name, e$message))
#     })
# 
#     # Export to LaTeX
#     tryCatch({
#       print(xtable(tab, caption = paste("Ergebnisse für", name), digits = 3),
#             file = paste0("Ergebnisse_", name, ".tex"))
#     }, error = function(e) {
#       message(sprintf("Error exporting LaTeX for %s: %s", name, e$message))
#     })
#   }
# }


## Calculate a linear mixed model  (reworked: unique SimUID, RegimeIndex as factor, safer optimizer)

# Check for factors with only one level (will cause issues with scaling and modeling)
lmm_factor_levels <- sapply(c("Timesteps", "Density", "Nodes", "Regimes"), function(var) {
  length(unique(corr_results[[var]]))
})

if (any(lmm_factor_levels == 1)) {
  single_level_vars <- names(lmm_factor_levels)[lmm_factor_levels == 1]
  message(sprintf(
    "Warning: Cannot fit linear mixed models - the following variables have only one level: %s",
    paste(single_level_vars, collapse = ", ")
  ))
  message("Skipping linear mixed model analysis. Run with multiple factor levels to enable this analysis.")
  linear_mixed_models <- list()
  
} else {
  
  # Prepare grouping + predictors
  corr_results <- corr_results %>%
    mutate(
      # RegimeIndex should not be treated as numeric trend
      RegimeIndex = factor(RegimeIndex),
      
      # Make a globally unique simulation identifier (ts_id restarts at 1 per condition)
      SimUID = interaction(Condition, SimID, drop = TRUE),
      
      # Scale numeric versions of the factors (treat as quantitative levels)
      Timesteps_scaled = as.numeric(scale(as.numeric(as.character(Timesteps)))),
      Density_scaled   = as.numeric(scale(as.numeric(as.character(Density)))),
      Nodes_scaled     = as.numeric(scale(as.numeric(as.character(Nodes)))),
      Regimes_scaled   = as.numeric(scale(as.numeric(as.character(Regimes))))
    )
  
  # Set variables to model
  cols <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")
  
  linear_mixed_models <- list()
  
  for (col_name in cols) {
    
    # Specify model (drop Condition random intercept to avoid singular fits)
    fmla <- as.formula(paste0(
      col_name,
      " ~ Timesteps_scaled * Density_scaled * Nodes_scaled * Regimes_scaled + ",
      "RegimeIndex + ",
      "(1 | SimUID)"
    ))
    
    model <- tryCatch({
      lmer(
        fmla,
        data = corr_results,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      )
    }, error = function(e) {
      message(sprintf("Error fitting LMM for %s: %s", col_name, e$message))
      return(NULL)
    }, warning = function(w) {
      # Keep the warning visible but continue (you can change this behavior if you want)
      message(sprintf("Warning fitting LMM for %s: %s", col_name, w$message))
      invokeRestart("muffleWarning")
    })
    
    if (is.null(model)) {
      linear_mixed_models[[col_name]] <- "Model failed - see error message above"
      next
    }
    
    # Summarize model
    lmm <- as.data.frame(coef(summary(model))) %>%
      dplyr::mutate(
        Signif = ifelse(`Pr(>|t|)` < 0.001, "***",
                        ifelse(`Pr(>|t|)` < 0.01, "**",
                               ifelse(`Pr(>|t|)` < 0.05, "*", "")))
      )
    
    rownames(lmm) <- gsub(":", " × ", gsub("_scaled", "", rownames(lmm)))
    
    linear_mixed_models[[col_name]] <- lmm
  }
}

# Export results for each model in the list
# Only export if we have valid model results (data frames, not error messages)
if (length(linear_mixed_models) > 0) {
  for (name in names(linear_mixed_models)) {
    tab <- linear_mixed_models[[name]]
    
    # Skip if this is an error message rather than a model
    if (is.character(tab)) {
      message(sprintf("Skipping export for %s: %s", name, tab))
      next
    }
    
    # Export to HTML
    tryCatch({
      kable(tab, caption = paste("Ergebnisse für", name), format = "html") %>%
        kable_styling(bootstrap_options = c("striped", "hover")) %>%
        save_kable(file = paste0("Ergebnisse_", name, ".html"))
    }, error = function(e) {
      message(sprintf("Error exporting HTML for %s: %s", name, e$message))
    })
    
    # Export to LaTeX (digits=3)
    tryCatch({
      print(
        xtable(tab, caption = paste("Ergebnisse für", name), digits = 3),
        file = paste0("Ergebnisse_", name, ".tex")
      )
    }, error = function(e) {
      message(sprintf("Error exporting LaTeX for %s: %s", name, e$message))
    })
  }
}


#######################################
### Make lineplot panels for each variable

# Set variables to plot
cols <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")

for (col_name in cols) {

  title <- switch(
    col_name,
    "Wtemp_corr"      = "Mean correlations for Wtemp",
    "Wtemp_ac_corr"   = "Mean correlations for Wtemp average controllability",
    "Wcont_corr"      = "Mean correlations for Wcont",
    "Wcont_ac_corr"   = "Mean correlations for Wcont average controllability",
  )

  # Summarize data
  summary_data <- corr_results %>%
    group_by(Timesteps, Density, Nodes, Regimes) %>%
    summarise(
      mean_val = mean(.data[[col_name]], na.rm = TRUE),
      sd_val   = sd(.data[[col_name]], na.rm = TRUE),
      .groups  = "drop"
    )

  # Make hover info
  summary_data <- summary_data %>%
    mutate(HoverInfo = paste("Timesteps:", Timesteps,
                             "<br>Nodes:", Nodes,
                             "<br>Mean:", round(mean_val, 2),
                             "<br>Sd:",   round(sd_val, 2)))

  # Position for shifted plots
  dodge <- position_dodge(width = 0.5)

  # Make plot
  p <- ggplot() +
    # Singular values as scatterplot
    geom_jitter(
      data = corr_results,
      aes_string(x = "Timesteps", y = col_name, color = "Nodes"),
      position = dodge, alpha = 0.2, size = 0.1
    ) +
    # Aggregated data as line plots
    geom_line(
      data = summary_data,
      aes(x = as.numeric(Timesteps), y = mean_val, color = Nodes, group = Nodes),
      position = dodge
    ) +
    # Aggregated data as scatterplots
    geom_point(
      data = summary_data,
      aes(x = as.numeric(Timesteps), y = mean_val, color = Nodes, text = HoverInfo),
      position = dodge, size = 1
    ) +
    # Facets
    facet_grid(Density ~ Regimes, labeller = label_value) +
    labs(
      title = title,
      x = "Timesteps",
      y = "Mean correlations",
      color = "Nodes"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
      strip.text.y = element_text(vjust = -0.25),
      panel.spacing = unit(0.2, "in"),
      plot.title = element_text(hjust = 0.5)
    )

  # Additional labels
  label_x_right <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Density", angle = -90, size = 4, hjust = 0)

  label_y_top <- ggplot() +
    theme_void() +
    xlim(0, 1) +
    ylim(0, 1) +
    annotate("text", x = 0.4725, y = 0.5, label = "Regimes", size = 4, hjust = 0.5)

  # Combine plots with `cowplot`
  final_plot <- ggdraw() +
    draw_plot(p, 0, 0, 1, 1) +
    draw_plot(label_x_right, 0.96, 0.08, 0.03, 0.8) +
    draw_plot(label_y_top,   0.1,  0.875, 0.84, 0.05)

  # Store plot
  output_file <- paste0("Plots/", col_name, "_plot_with_labels.pdf")

  ggsave(
    filename = output_file,
    plot     = final_plot,
    width    = 10,
    height   = 8,
    units    = "in",
    dpi      = 600
  )

  message("Gespeichert: ", output_file)
}






#######################################
library(lme4)
# optional: library(performance)  # für r2_nakagawa

outcomes <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")

# Falls noch nicht vorhanden:
corr_results <- corr_results %>%
  mutate(
    SimUID = interaction(Condition, SimID, drop = TRUE),
    RegimeIndex = factor(RegimeIndex)
  )

fits <- list()

for (y in outcomes) {
  
  f_with <- as.formula(paste0(
    y, " ~ Timesteps_scaled * Density_scaled * Nodes_scaled * Regimes_scaled + ",
    "RegimeIndex + (1 | SimUID)"
  ))
  
  f_without <- as.formula(paste0(
    y, " ~ Timesteps_scaled * Density_scaled * Nodes_scaled * Regimes_scaled + ",
    "(1 | SimUID)"
  ))
  
  m_with <- lmer(f_with, data = corr_results,
                 REML = FALSE,
                 control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  
  m_without <- lmer(f_without, data = corr_results,
                    REML = FALSE,
                    control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  
  comp <- anova(m_without, m_with)   # Likelihood-Ratio-Test
  aic  <- AIC(m_without, m_with)
  bic  <- BIC(m_without, m_with)
  
  # optional:
  # r2 <- performance::r2_nakagawa(m_with); r2_wo <- performance::r2_nakagawa(m_without)
  
  fits[[y]] <- list(
    LRT = comp,
    AIC = aic,
    BIC = bic
    #, R2_with = r2, R2_without = r2_wo
  )
}

fits

