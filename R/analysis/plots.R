# =============================================================================
# Graphical visualisations: line plots, sensitivity/specificity, coefficient plot
# =============================================================================
# Extracted from the former monolithic stat_analysis.R (PARTs 7, 7.1, 8).

# -----------------------------------------------------------------------------
# HOUSE STYLE: one shared palette + theme for every figure in the thesis
# -----------------------------------------------------------------------------
# Single source of truth: change the palette/theme HERE and it applies to
# every figure the pipeline produces, so the thesis reads as one consistent
# document instead of a patchwork of one-off plot styles.
#
# Palette: viridis. Perceptually uniform AND distinguishable under all common
# forms of colour-vision deficiency (deuteranopia/protanopia/tritanopia) --
# unlike ggplot2's default hue palette, which places a red and a green at
# very similar luminance right next to each other (the classic problem case).
# `begin`/`end` stay inside (0, 1) so the lightest yellow and darkest purple
# (low-contrast against a white page / the panel grid, respectively) are
# excluded.
THESIS_VIRIDIS_RANGE <- c(0.08, 0.85)

scale_color_thesis_d <- function(...) {
  scale_color_viridis_d(begin = THESIS_VIRIDIS_RANGE[1], end = THESIS_VIRIDIS_RANGE[2], ...)
}
scale_fill_thesis_d <- function(...) {
  scale_fill_viridis_d(begin = THESIS_VIRIDIS_RANGE[1], end = THESIS_VIRIDIS_RANGE[2], ...)
}

# Redundant shape coding for Nodes, IN ADDITION to colour: keeps figures
# legible in greyscale printouts/photocopies and gives a second,
# colour-independent cue on top of the colourblind-safe palette -- the
# standard "don't rely on colour alone" accessibility recommendation.
NODE_SHAPES <- c("4" = 16, "6" = 17, "8" = 15)  # filled circle / triangle / square

# One shared theme, used by every figure below. Close to theme_minimal()
# (clean, low ink-to-data ratio, prints well) with slightly bolder
# titles/strips so figures stay legible at thesis page size.
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position   = "top",
      legend.title      = element_text(face = "bold", size = rel(0.95)),
      legend.text       = element_text(size = rel(0.85)),
      strip.text        = element_text(face = "bold", size = rel(0.85)),
      strip.background  = element_rect(fill = "grey93", color = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey88", linewidth = 0.3),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = rel(1.15)),
      plot.subtitle     = element_text(hjust = 0.5, color = "grey35", size = rel(0.9)),
      plot.caption      = element_text(color = "grey45", size = rel(0.75), hjust = 0),
      axis.title        = element_text(face = "bold"),
      axis.text.x       = element_text(angle = 45, hjust = 1)
    )
}

# -----------------------------------------------------------------------------
# PART 7: Line plots with facets (one per correlation outcome)
# -----------------------------------------------------------------------------
make_line_plots <- function(corr_results, cols) {

  for (col_name in cols) {

    title <- switch(
      col_name,
      "Beta_corr"             = "Wiederherstellung des temporalen Netzwerks (Beta)",
      "Kappa_corr"            = "Wiederherstellung des kontempor\u00e4ren Netzwerks (Kappa)",
      "Beta_ac_corr_pearson"  = "Wiederherstellung der durchschnittlichen Kontrollierbarkeit (Pearson)",
      "Beta_ac_corr_spearman" = "Wiederherstellung der durchschnittlichen Kontrollierbarkeit (Spearman)"
    )

    # Boxplots directly on the raw per-regime values: median + IQR (box) +
    # whiskers/outliers, computed empirically with NO transform and NO
    # parametric (Gaussian-shape) assumption. This replaces the earlier
    # mean +/- SD approach entirely -- not just the raw-r version (which was
    # the wrong SHAPE near the boundary) but also the Fisher-z version (which
    # fixed the shape but is highly sensitive to the few near-|1| replicates:
    # atanh(r) -> Inf as r -> 1, so a handful of near-perfect fits can drag
    # the back-transformed mean far above where most of the data actually is
    # -- e.g. one design cell had raw mean 0.50 but z-mean 0.84, because 75%
    # of replicates clustered near 1 and 25% were near zero/negative. Medians
    # and quantiles are immune to this: median(f(x)) = f(median(x)) for any
    # monotonic f, so there's no "which scale do I average on" ambiguity, and
    # a few extreme replicates can only shift which side of the median they
    # fall on, never drag it toward them.
    dodge <- position_dodge(width = 0.7)
    box_width <- 0.7 / 3 * 0.85  # 3 Nodes-Gruppen je Zeitschritt, mit kleinem Zwischenraum

    p <- ggplot(corr_results, aes(x = Timesteps, y = .data[[col_name]],
                                  color = Nodes, fill = Nodes)) +
      geom_boxplot(
        position = dodge, width = box_width, linewidth = 0.4,
        alpha = 0.25, outlier.size = 0.6, outlier.alpha = 0.4
      ) +
      # Thin connecting line through the medians, so the trend across T
      # stays as easy to follow at a glance as in the old line plot.
      stat_summary(
        fun = median, geom = "line", position = dodge,
        aes(group = Nodes), linewidth = 0.6
      ) +
      coord_cartesian(ylim = c(NA, 1)) +
      facet_grid(Density ~ Regimes,
                 labeller = labeller(
                   Density = function(x) paste0("Dichte: ", x),
                   Regimes = function(x) paste0("Regime: ", x)
                 )) +
      scale_color_thesis_d() +
      scale_fill_thesis_d() +
      labs(
        title    = title,
        subtitle = "Boxplot \u00fcber alle Replikationen: Median, IQR, Whisker (1.5\u00d7IQR), Ausrei\u00dfer",
        x        = "Zeitschritte (T)",
        y        = "Korrelation",
        color    = "Nodes",
        fill     = "Nodes"
      ) +
      theme_thesis() +
      theme(
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.spacing = unit(0.15, "in")
      )

    output_file <- paste0("Plots/", col_name, "_plot_with_labels.pdf")

    # cairo_pdf instead of the default pdf() device: the base device's font-
    # metric calculation assumes a single-byte locale and mangles multi-byte
    # UTF-8 characters (umlauts, sharp s, multiplication sign, plus-minus,
    # ...) when R is running under a C/POSIX locale (as in this sandbox) --
    # cairo renders the glyphs directly and isn't affected.
    ggsave(
      filename = output_file,
      plot     = p,
      device   = cairo_pdf,
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
      geom_errorbar(
        data = summary_data,
        aes(x = Timesteps,
            ymin = pmax(mean_val - sd_val, 0), ymax = pmin(mean_val + sd_val, 1),
            color = Nodes, group = interaction(Nodes, Metric)),
        width = 0.3, linewidth = 0.4, position = dodge
      ) +
      geom_line(
        data = summary_data,
        aes(x = Timesteps, y = mean_val,
            color = Nodes, group = interaction(Nodes, Metric), linetype = Metric),
        position = dodge
      ) +
      geom_point(
        data = summary_data,
        aes(x = Timesteps, y = mean_val, color = Nodes, shape = Nodes),
        position = dodge, size = 1.8
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      facet_grid(Density ~ Regimes,
                 labeller = labeller(
                   Density = function(x) paste0("Dichte: ", x),
                   Regimes = function(x) paste0("Regime: ", x)
                 )) +
      scale_color_thesis_d() +
      scale_shape_manual(values = NODE_SHAPES) +
      labs(
        title    = paste0("Sensitivit\u00e4t & Spezifit\u00e4t: ",
                          if (net == "Beta") "temporales Netzwerk (Beta)" else "kontempor\u00e4res Netzwerk (Kappa)"),
        subtitle = "Fehlerbalken: \u00b11 SD (auf [0, 1] gekappt)",
        x        = "Zeitschritte (T)",
        y        = "Mittlerer Wert",
        color    = "Nodes",
        shape    = "Nodes",
        linetype = "Metrik"
      ) +
      theme_thesis() +
      theme(
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.spacing = unit(0.15, "in")
      )

    output_file <- paste0("Plots/", net, "_senspec_plot_with_labels.pdf")

    ggsave(
      filename = output_file,
      plot     = p,
      device   = cairo_pdf,
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

  # Readable labels for outcomes (matches the r_Beta / r_Kappa / r_AC
  # notation used in the thesis text and comparison_all.tex)
  outcome_labels <- c(
    Beta_corr     = "Beta (Korrelation)",
    Kappa_corr    = "Kappa (Korrelation)",
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
  coef_plot_data$predictor <- gsub(":",    " \u00d7 ", coef_plot_data$predictor)
  coef_plot_data$predictor <- gsub("_s\\b", "",   coef_plot_data$predictor)
  coef_plot_data$predictor <- gsub("_num\\b", "",  coef_plot_data$predictor)

  # Compute 95% CI and significance flag
  coef_plot_data <- coef_plot_data %>%
    mutate(
      ci_lo  = estimate - 1.96 * se,
      ci_hi  = estimate + 1.96 * se,
      sig    = ifelse(p_adj < 0.05, "p < .05", "n.s."),
      # Order predictors: intercept last, main effects first, then interactions by number of "x"
      n_cross = stringr::str_count(predictor, "\u00d7")
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
    scale_color_thesis_d() +
    facet_wrap(~ n_cross, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(
                 "0" = "Haupteffekte",
                 "1" = "2-Wege-Interaktionen",
                 "2" = "3-Wege-Interaktionen"
               ))) +
    labs(
      title  = "Regressionskoeffizienten der prim\u00e4ren 3-Wege-Modelle",
      x      = "Regressionskoeffizient (95%-KI)",
      y      = NULL,
      color  = "Outcome",
      shape  = "Signifikanz"
    ) +
    theme_thesis(base_size = 10) +
    theme(
      legend.box       = "vertical",
      panel.grid.minor = element_blank()
    )

  ggsave("Plots/coefficient_plot_full.pdf",
         plot = p_coef_full, device = cairo_pdf, width = 10, height = 14, dpi = 600)

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
    scale_color_thesis_d() +
    labs(
      x      = "Regressionskoeffizient (95%-KI)",
      y      = NULL,
      color  = "Outcome",
      shape  = "Signifikanz",
      title  = "Haupteffekte"
    ) +
    theme_thesis(base_size = 11) +
    theme(
      legend.box       = "vertical",
      panel.grid.minor = element_blank()
    )

  ggsave("Plots/coefficient_plot_main_effects.pdf",
         plot = p_coef_main, device = cairo_pdf, width = 8, height = 5, dpi = 600)

  message("Coefficient plots saved to Plots/")
}
