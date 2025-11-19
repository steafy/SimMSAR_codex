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
## Extract stats data from MSAR results list

# PERFORMANCE: Pre-allocate list to avoid O(n²) list copying
# Estimate max size conservatively (actual size may be smaller due to failed fits)
max_size <- length(T) * length(Density) * length(N) * length(M) * 30 * max(M)
corr_list <- vector("list", max_size)
list_idx <- 0

for (t in seq_along(T)) {
   for (density in seq_along(Density)) {
    for (nodes in seq_along(N)) {
      for (regimes in seq_along(M)) {
        models <- MSAR_models[[t]][[density]][[nodes]][[regimes]][["MSAR_models"]]
        for (ts in seq_along(models)) {
          for (r in seq_along(models[[ts]])) {
            result <- models[[ts]][[r]]
            list_idx <- list_idx + 1
            temp <- data.frame(
              Timesteps = T[t],
              Density = Density[density],
              Nodes = N[nodes],
              Regimes = M[regimes],
              N = length(models),
              Wtemp_corr = result[["corr. Wtemp"]],
              Wtemp_ac_corr = result[["corr. Wtemp ac"]],
              Wcont_corr = result[["corr. Wcont"]],
              Wcont_ac_corr = result[["corr. Wcont ac"]]
            )
            corr_list[[list_idx]] <- temp
          }
        }
      }
    }
  }
}

# PERFORMANCE: Trim to actual size
corr_list <- corr_list[1:list_idx]


# Combine all dataframes
corr_results <- do.call(rbind, corr_list)

# Remove NA
corr_results <- na.omit(corr_results)

# Group regimes from one timeseries together, add SimID
corr_results <- corr_results %>%
  group_by(Timesteps, Density, Nodes, Regimes) %>%
  mutate(
    SimID = rep(seq_len(ceiling(n() / first(as.numeric(as.character(Regimes))))),
                each = first(as.numeric(as.character(Regimes))), length.out = n()),
    RegimeIndex = rep(seq_len(first(as.numeric(as.character(Regimes)))), length.out = n())
  ) %>%
  ungroup() %>%
  mutate(Condition = interaction(Timesteps, Density, Nodes, Regimes, drop = TRUE)) %>%
  select(SimID, Timesteps, Density, Nodes, Regimes, RegimeIndex, everything())

# Transform to factors
corr_results <- corr_results %>%
  mutate(
    Timesteps = factor(Timesteps),
    Density = factor(Density),
    Nodes = factor(Nodes),
    Regimes = factor(Regimes)
  )




## Extract stats from MSAR_models
descript_stats <- list()

names <- c("Wtemp_corr",
           "Wtemp_MAE",
           "Wtemp_sensitivity",
           "Wtemp_specificity",
           "Wtemp_ac_corr",
           "Wcont_corr",
           "Wcont_MAE",
           "Wcont_sensitivity",
           "Wcont_specificity",
           "Wcont_ac_corr")


# Extract all combinations
all_combos <- expand.grid(
  T = T,
  Density = Density,
  N = N,
  M = M,
  stringsAsFactors = FALSE
)

# Make dataframe for every outcome variable
for (current_name in names) {
  
  # Construct dataframe
  df_temp <- data.frame(
    Timesteps = numeric(0),
    Density   = numeric(0),
    Nodes     = numeric(0),
    Regimes   = numeric(0),
    Value   = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(all_combos))) {
    
    t_idx  <- match(all_combos$T[i], T)
    d_idx  <- match(all_combos$Density[i], Density)
    n_idx  <- match(all_combos$N[i], N)
    m_idx  <- match(all_combos$M[i], M)

    # Extract stats from MSAR_models
    stats <- MSAR_models[[t_idx]][[d_idx]][[n_idx]][[m_idx]][["Stats"]]
    val   <- stats[[current_name]]
    val_df <- as.data.frame(as.list(val), stringsAsFactors = FALSE)
    
    row_df <- data.frame(
      Timesteps = all_combos$T[i],
      Density   = all_combos$Density[i],
      Nodes     = all_combos$N[i],
      Regimes    = all_combos$M[i],
      stringsAsFactors = FALSE
    )
    
    row_df <- cbind(row_df, val_df)
    df_temp <- rbind(df_temp, row_df)

  }

  # Store results
  descript_stats[[current_name]] <- df_temp
}


####################################################################
### Calculate statistical information for estimation process

## Determine no. of omissions (ts without converging model) (30 - N)
aggr_factors <- corr_results %>%
  dplyr::select(Timesteps, Density, Nodes, Regimes, N) %>%
  dplyr::distinct()

omissions <- 9720 - sum(aggr_factors$N)
omissions_per <- (omissions / 9720) * 100


## Calculate mean no. of estimated models (N) for factorlevels
## Calculate dunn-test to compare N across factorlevels
n_means <- list()
dunn_results <- list()
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
    kw <- kruskal.test(aggr_factors$N, aggr_factors[[var]])
    dunn <- dunn.test(aggr_factors$N, aggr_factors[[var]], method = "bonferroni")
    dunn_matrix <- data.frame(
      comp = dunn$comparisons,
      Z_val = dunn$Z,
      p_val = dunn$P,
      p.adj = round(dunn$P.adjusted, 4)
    )
    dunn_matrix <- dunn_matrix[order(dunn_matrix$p.adj), ]
    dunn_results[[var]] <- list(
      Kruskal = kw,
      Dunn = dunn_matrix
    )
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
  means <- c()
  for (j in 6:8) {
    x <- mean(descript_stats[[i]][[j]])
    means <- c(means, x)
    val <- c(means,
             min(descript_stats[[i]][["Min"]]),
             max(descript_stats[[i]][["Max"]])
             )
  }
  desc[ ,i] <- val
}


###############################################
### Calculate inferential statistics to describe models

## Calculate PERMANOVA
# Specify dependent and independent variables
dependent_vars <- colnames(corr_results)[8:11]
independent_vars <- colnames(corr_results)[2:5]

permanova_results <- list()

# Check if there are multiple levels in each factor
factor_levels <- sapply(independent_vars, function(var) {
  length(unique(corr_results[[var]]))
})

if (any(factor_levels == 1)) {
  single_level_vars <- names(factor_levels)[factor_levels == 1]
  message(sprintf("Warning: The following variables have only one level: %s",
                  paste(single_level_vars, collapse = ", ")))
  message("PERMANOVA results may be limited. Consider running with multiple factor levels.")
}

# Loop over dependent variables
for (dep_var in dependent_vars) {
  # Create formula for ANOVA
  formula <- as.formula(paste(dep_var, "~", paste(independent_vars, collapse = " + "), "+",
                              paste(combn(independent_vars, 2, FUN = paste, collapse = ":"), collapse = " + "), "+",
                              paste(combn(independent_vars, 3, FUN = paste, collapse = ":"), collapse = " + "), "+",
                              paste(independent_vars, collapse = ":"), collapse = " "))

  # Calculate permutations ANOVA
  model <- tryCatch({
    aovp(formula, data = corr_results, perm = "Prob", maxIter = 5000)
  }, error = function(e) {
    message(sprintf("Error in PERMANOVA for %s: %s", dep_var, e$message))
    return(NULL)
  })

  # Skip if model failed
  if (is.null(model)) {
    permanova_results[[dep_var]] <- "Model failed - see error message above"
    next
  }

  # Extract summary
  summary_model <- summary(model)
  
  # Extract effects, degrees of freedom, and mean squares.n
  effects <- rownames(summary_model[[1]])
  Df <- summary_model[[1]][, "Df"]
  R_Sum_Sq <- summary_model[[1]][, "R Sum Sq"]
  R_Mean_Sq <- summary_model[[1]][, "R Mean Sq"]
  
  # Extract residual Mean Square and p-values
  effects <- trimws(effects) 
  residual_mean_sq <- R_Mean_Sq[effects == "Residuals"]
  if (length(residual_mean_sq) == 0 || is.na(residual_mean_sq)) {
    stop("Residual Mean Square konnte nicht berechnet werden.")
  }

  # Calculate F-values without residuals
  f_values <- R_Mean_Sq / residual_mean_sq
  f_values[effects == "Residuals"] <- NA  
  
  summary_model[[1]]$F_values <- round(f_values, 3)
  
  # Store results in list
  permanova_results[[dep_var]] <- summary_model
}

# Show results
for (dep_var in names(permanova_results)) {
  cat("\nResults for dependent variable:", dep_var, "\n")
  print(permanova_results[[dep_var]])
}



###################################
### Calculate a linear mixed model

# Check for factors with only one level (will cause issues with scaling and modeling)
lmm_factor_levels <- sapply(c("Timesteps", "Density", "Nodes", "Regimes"), function(var) {
  length(unique(corr_results[[var]]))
})

if (any(lmm_factor_levels == 1)) {
  single_level_vars <- names(lmm_factor_levels)[lmm_factor_levels == 1]
  message(sprintf("Warning: Cannot fit linear mixed models - the following variables have only one level: %s",
                  paste(single_level_vars, collapse = ", ")))
  message("Skipping linear mixed model analysis. Run with multiple factor levels to enable this analysis.")
  linear_mixed_models <- list()
} else {
  # Scale factors
  corr_results <- corr_results %>%
    mutate(
      Timesteps_scaled = scale(as.numeric(as.character(Timesteps))),
      Density_scaled   = scale(as.numeric(as.character(Density))),
      Nodes_scaled     = scale(as.numeric(as.character(Nodes))),
      Regimes_scaled   = scale(as.numeric(as.character(Regimes)))
    )

  cols <- colnames(corr_results)[8:11]

  linear_mixed_models <- list()

  for (col_name in cols) {

    # Specify model
    fmla <- as.formula(paste0(col_name, " ~ Timesteps_scaled * Density_scaled * Nodes_scaled * Regimes_scaled +
                            RegimeIndex +
                            (1 | Condition/SimID)")
                     )

    model <- tryCatch({
      lmer(fmla, data = corr_results)
    }, error = function(e) {
      message(sprintf("Error fitting LMM for %s: %s", col_name, e$message))
      return(NULL)
    })

    if (is.null(model)) {
      linear_mixed_models[[col_name]] <- "Model failed - see error message above"
      next
    }

    # Summarize model
    lmm <- as.data.frame(coef(summary(model))) %>%
      dplyr::mutate(Signif = ifelse(`Pr(>|t|)` < 0.001, "***",
                                    ifelse(`Pr(>|t|)` < 0.01, "**",
                                           ifelse(`Pr(>|t|)` < 0.05, "*", ""))))

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

    # Export to LaTeX
    tryCatch({
      print(xtable(tab, caption = paste("Ergebnisse für", name)),
            file = paste0("Ergebnisse_", name, ".tex"))
    }, error = function(e) {
      message(sprintf("Error exporting LaTeX for %s: %s", name, e$message))
    })
  }
}

#######################################
### Make lineplot panels for each variable

# Set variables to plot
cols <- colnames(corr_results)[8:11]

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


