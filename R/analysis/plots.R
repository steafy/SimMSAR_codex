# =============================================================================
# Graphical visualisations: line plots, sensitivity/specificity, coefficient plot
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PARTs 7, 7.1, 8).

# -----------------------------------------------------------------------------
# PART 7: Line plots with facets (one per correlation outcome)
# -----------------------------------------------------------------------------
make_line_plots <- function(corr_results, cols) {

  for (col_name in cols) {

    title <- switch(
      col_name,
      "Beta_corr"       = "Mean correlations for Beta",
      "Kappa_corr"      = "Mean correlations for Kappa",
      "Beta_ac_corr_pearson"  = "Mean correlations for Beta AC (Pearson)",
      "Beta_ac_corr_spearman" = "Mean correlations for Beta AC (Spearman)"
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
        aes(x = Timesteps, y = .data[[col_name]], color = Nodes),
        position = dodge, alpha = 0.2, size = 0.1
      ) +
      # Aggregated data as line plots
      geom_line(
        data = summary_data,
        aes(x = as.numeric(as.character(Timesteps)), y = mean_val, color = Nodes, group = Nodes),
        position = dodge
      ) +
      # Aggregated data as scatterplots
      geom_point(
        data = summary_data,
        aes(x = as.numeric(as.character(Timesteps)), y = mean_val, color = Nodes, text = HoverInfo),
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

    message("saved: ", output_file)
  }
}

# -----------------------------------------------------------------------------
# PART 7.1: Sensitivity/Specificity (DESCRIPTIVE ONLY -- no LMMs)
# -----------------------------------------------------------------------------
# Sens/Spec are reported descriptively (aggregated at simulation level, mean
# across regimes, consistent with the rest of the aggregation), not modeled
# via LMMs. One facet plot per network (Beta, Kappa), each showing both
# sensitivity and specificity (distinguished by linetype).
make_senspec_plots <- function(corr_results) {

  for (net in c("Beta", "Kappa")) {

    sen_col  <- paste0(net, "_sen")
    spec_col <- paste0(net, "_spec")

    base_cols <- corr_results %>% select(Timesteps, Density, Nodes, Regimes)
    senspec_long <- bind_rows(
      base_cols %>% mutate(Metric = "Sensitivity", Value = corr_results[[sen_col]]),
      base_cols %>% mutate(Metric = "Specificity", Value = corr_results[[spec_col]])
    )

    summary_data <- senspec_long %>%
      group_by(Timesteps, Density, Nodes, Regimes, Metric) %>%
      summarise(
        mean_val = mean(Value, na.rm = TRUE),
        sd_val   = sd(Value, na.rm = TRUE),
        .groups  = "drop"
      )

    dodge <- position_dodge(width = 0.5)

    p <- ggplot() +
      geom_jitter(
        data = senspec_long,
        aes(x = Timesteps, y = Value, color = Nodes),
        position = dodge, alpha = 0.15, size = 0.1
      ) +
      geom_line(
        data = summary_data,
        aes(x = as.numeric(as.character(Timesteps)), y = mean_val,
            color = Nodes, group = interaction(Nodes, Metric), linetype = Metric),
        position = dodge
      ) +
      geom_point(
        data = summary_data,
        aes(x = as.numeric(as.character(Timesteps)), y = mean_val, color = Nodes),
        position = dodge, size = 1
      ) +
      facet_grid(Density ~ Regimes, labeller = label_value) +
      labs(
        title = paste("Mean sensitivity/specificity for", net),
        x = "Timesteps",
        y = "Mean value",
        color = "Nodes",
        linetype = "Metric"
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

    label_x_right <- ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Density", angle = -90, size = 4, hjust = 0)

    label_y_top <- ggplot() +
      theme_void() +
      xlim(0, 1) +
      ylim(0, 1) +
      annotate("text", x = 0.4725, y = 0.5, label = "Regimes", size = 4, hjust = 0.5)

    final_plot <- ggdraw() +
      draw_plot(p, 0, 0, 1, 1) +
      draw_plot(label_x_right, 0.96, 0.08, 0.03, 0.8) +
      draw_plot(label_y_top,   0.1,  0.875, 0.84, 0.05)

    output_file <- paste0("Plots/", net, "_senspec_plot_with_labels.pdf")

    ggsave(
      filename = output_file,
      plot     = final_plot,
      width    = 10,
      height   = 8,
      units    = "in",
      dpi      = 600
    )

    message("saved: ", output_file)
  }
}

# -----------------------------------------------------------------------------
# PART 8: Coefficient plot (full + main-effects-only)
# -----------------------------------------------------------------------------
make_coefficient_plots <- function(all_results, corr_cols) {

  # Readable labels for outcomes
  outcome_labels <- c(
    Beta_corr     = "Beta (Correlation)",
    Kappa_corr    = "Kappa (Correlation)",
    Beta_ac_corr_pearson  = "Beta AC (Pearson)",
    Beta_ac_corr_spearman = "Beta AC (Spearman)"
  )

  # Collect coefficients + CIs from all primary models
  coef_plot_data <- lapply(corr_cols, function(outcome) {
    model        <- all_results[[outcome]]$primary_model
    coef_summary <- summary(model)$coefficients

    data.frame(
      outcome   = outcome_labels[outcome],
      predictor = rownames(coef_summary),
      estimate  = coef_summary[, "Estimate"],
      se        = coef_summary[, "Std. Error"],
      p_adj     = p.adjust(coef_summary[, "Pr(>|t|)"], method = "BH"),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  # Clean predictor names (same as export)
  coef_plot_data$predictor <- gsub(":",    " × ", coef_plot_data$predictor)
  coef_plot_data$predictor <- gsub("_s\\b", "",   coef_plot_data$predictor)
  coef_plot_data$predictor <- gsub("_num\\b", "",  coef_plot_data$predictor)

  # Compute 95% CI and significance flag
  coef_plot_data <- coef_plot_data %>%
    mutate(
      ci_lo  = estimate - 1.96 * se,
      ci_hi  = estimate + 1.96 * se,
      sig    = ifelse(p_adj < 0.05, "p < .05", "n.s."),
      # Order predictors: intercept last, main effects first, then interactions by n of ×
      n_cross = stringr::str_count(predictor, "×")
    )

  # Exclude intercept from plot (usually not of interest)
  coef_plot_data <- coef_plot_data %>%
    filter(predictor != "(Intercept)")

  # Order predictors by complexity then alphabetically
  predictor_order <- coef_plot_data %>%
    arrange(n_cross, predictor) %>%
    pull(predictor) %>%
    unique()

  coef_plot_data$predictor <- factor(coef_plot_data$predictor, levels = rev(predictor_order))

  # --- 8.1: Full coefficient plot (all terms) ---
  p_coef_full <- ggplot(coef_plot_data,
                        aes(x = estimate, y = predictor,
                            color = outcome, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1)) +
    scale_color_viridis_d(option = "D", end = 0.85) +
    facet_wrap(~ n_cross, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(
                 "0" = "Main Effects",
                 "1" = "2-Way Interactions",
                 "2" = "3-Way Interactions"
               ))) +
    labs(
      x      = "Regression Coefficient (95% CI)",
      y      = NULL,
      color  = "Outcome",
      shape  = "Significance"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "bottom",
      legend.box       = "vertical",
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave("Plots/coefficient_plot_full.pdf",
         plot = p_coef_full, width = 10, height = 14, dpi = 600)

  # --- 8.2: Main effects only (compact version for main text) ---
  p_coef_main <- coef_plot_data %>%
    filter(n_cross == 0) %>%
    ggplot(aes(x = estimate, y = predictor,
               color = outcome, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1)) +
    scale_color_viridis_d(option = "D", end = 0.85) +
    labs(
      x      = "Regression Coefficient (95% CI)",
      y      = NULL,
      color  = "Outcome",
      shape  = "Significance",
      title  = "Main Effects"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.box       = "vertical",
      panel.grid.minor = element_blank()
    )

  ggsave("Plots/coefficient_plot_main_effects.pdf",
         plot = p_coef_main, width = 8, height = 5, dpi = 600)

  message("Coefficient plots saved to Plots/")
}
