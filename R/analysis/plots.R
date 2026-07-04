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
# --- Typography --------------------------------------------------------------
# Base font for every figure. A clean, geometric-ish sans in the spirit of the
# Poppins/Inter/Manrope references, using what is actually installed on the
# machine (none of those three are). "Segoe UI" is the modern Windows humanist
# sans and renders correctly through cairo_pdf; the fallbacks cover other boxes.
THESIS_FONT <- local({
  installed <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character(0))
  pick <- c("Poppins", "Inter", "Manrope", "Segoe UI", "Calibri")
  hit  <- pick[pick %in% installed]
  if (length(hit)) hit[1] else ""   # "" => graphics-device default sans
})

# --- Panel "cards" -----------------------------------------------------------
# Soft panel fill on a white plot background, so each panel reads as a distinct
# rounded card (matching the reference dashboards). Rounded corners come from
# elementalist::element_rect_round(); it is GitHub-only
# (remotes::install_github("teunbrand/elementalist")), so if it is not available
# we fall back to a flat soft-fill rect and the pipeline still runs.
PANEL_FILL      <- "#F6F6F9"
PANEL_BORDER    <- "grey90"  # hairline around each card
PANEL_BORDER_LW <- 0.2       # border thickness (bump for a bolder frame)
PANEL_RADIUS    <- 6         # corner radius in pt
.have_elementalist <- requireNamespace("elementalist", quietly = TRUE)
if (.have_elementalist) suppressPackageStartupMessages(library(elementalist))

# The card FILL goes on panel.background and the STROKE on panel.border, on
# purpose. panel.background is drawn clipped to the panel rectangle, so a stroke
# there loses its straight edges (only the inward-curving corners survive); a
# border-only element on color = NA has nothing to clip. panel.border is drawn
# UNclipped, so the rounded stroke shows on all four edges. Splitting them is
# what makes the full hairline appear.
panel_card_bg <- function(radius = PANEL_RADIUS) {
  if (.have_elementalist) {
    elementalist::element_rect_round(radius = unit(radius, "pt"), fill = PANEL_FILL, color = NA)
  } else {
    element_rect(fill = PANEL_FILL, color = NA)
  }
}
panel_card_border <- function(radius = PANEL_RADIUS) {
  if (.have_elementalist) {
    elementalist::element_rect_round(radius = unit(radius, "pt"), fill = NA,
                                     color = PANEL_BORDER, linewidth = PANEL_BORDER_LW)
  } else {
    element_rect(fill = NA, color = PANEL_BORDER, linewidth = PANEL_BORDER_LW)
  }
}

# --- Colour + shape ----------------------------------------------------------
# viridis: perceptually uniform AND distinguishable under all common forms of
# colour-vision deficiency. Kept for 4-level mappings (the 4 outcomes / Regimes),
# where it is already colourblind-verified. begin/end stay inside (0, 1) so the
# lightest yellow / darkest purple (low-contrast on white / on the panel fill)
# are excluded.
THESIS_VIRIDIS_RANGE <- c(0.08, 0.85)

scale_color_thesis_d <- function(...) {
  scale_color_viridis_d(begin = THESIS_VIRIDIS_RANGE[1], end = THESIS_VIRIDIS_RANGE[2], ...)
}
scale_fill_thesis_d <- function(...) {
  scale_fill_viridis_d(begin = THESIS_VIRIDIS_RANGE[1], end = THESIS_VIRIDIS_RANGE[2], ...)
}

# Okabe-Ito colourblind-safe triple for the 3-level Nodes mapping (blue / orange
# / green). Distinct hues at well-separated luminance under deuteranopia /
# protanopia / tritanopia. Legend entries use the actual node counts ("N = 4"),
# mirroring the "M = 1" convention for Regimes.
OKABE_ITO_NODES <- c("4" = "#0072B2", "6" = "#E69F00", "8" = "#009E73")
NODE_LABELS     <- c("4" = "N = 4",   "6" = "N = 6",   "8" = "N = 8")

scale_color_thesis_nodes <- function(...) {
  scale_color_manual(values = OKABE_ITO_NODES, labels = NODE_LABELS, ...)
}
scale_fill_thesis_nodes <- function(...) {
  scale_fill_manual(values = OKABE_ITO_NODES, labels = NODE_LABELS, ...)
}

# Redundant shape coding for Nodes, IN ADDITION to colour: keeps figures legible
# in greyscale printouts/photocopies -- the "don't rely on colour alone"
# recommendation. Fillable variants (circle/triangle/square) so points can be
# drawn as a coloured fill with a white stroke (shape 21/24/22).
NODE_SHAPES <- c("4" = 21, "6" = 24, "8" = 22)

# --- Facet-strip labels ------------------------------------------------------
# Density on rows (right-hand strips), labelled directly as low/mid/high so no
# separate super-label is needed. Regimes on columns as "M = 1" ... "M = 4".
DENSITY_LABELS <- c("0.25" = "niedrige Dichte",
                    "0.5"  = "mittlere Dichte",
                    "0.75" = "hohe Dichte")

density_labeller <- function(x) {
  out <- DENSITY_LABELS[as.character(x)]
  ifelse(is.na(out), paste0("Dichte: ", x), out)
}
regime_labeller <- function(x) paste0("M = ", x)

# --- Shared theme ------------------------------------------------------------
# One theme, used by every figure below. theme_minimal() base (clean, low
# ink-to-data ratio) refreshed to match the reference dashboards: rounded
# soft-fill panel cards on white, plain bold strip text (no grey box), a single
# faint dotted horizontal guide, and the geometric base font.
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = THESIS_FONT) +
    theme(
      text               = element_text(family = THESIS_FONT),
      legend.position    = "top",
      legend.title       = element_text(family = THESIS_FONT, face = "bold", size = rel(0.95)),
      legend.text        = element_text(family = THESIS_FONT, size = rel(0.85)),
      strip.text         = element_text(family = THESIS_FONT, face = "bold", size = rel(0.9)),
      strip.background   = element_blank(),   # plain-text headers, no grey box
      strip.placement    = "outside",
      panel.background   = panel_card_bg(),     # rounded soft-fill "card"
      panel.border       = panel_card_border(), # rounded hairline frame (unclipped)
      plot.background     = element_rect(fill = "white", color = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3, linetype = "dotted"),
      panel.spacing      = unit(0.9, "lines"),
      plot.title         = element_text(family = THESIS_FONT, face = "bold", hjust = 0.5, size = rel(1.15)),
      plot.subtitle      = element_text(family = THESIS_FONT, hjust = 0.5, color = "grey35", size = rel(0.9)),
      plot.caption       = element_text(family = THESIS_FONT, color = "grey45", size = rel(0.75), hjust = 0),
      axis.title         = element_text(family = THESIS_FONT, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1)
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
        aes(group = Nodes), linewidth = 0.6,
        lineend = "round", linejoin = "round"
      ) +
      coord_cartesian(ylim = c(NA, 1)) +
      facet_grid(Density ~ Regimes,
                 labeller = labeller(
                   Density = density_labeller,
                   Regimes = regime_labeller
                 )) +
      scale_color_thesis_nodes() +
      scale_fill_thesis_nodes() +
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

    output_file <- plots_file(paste0(col_name, "_plot_with_labels.pdf"))

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
        position = dodge, lineend = "round", linejoin = "round"
      ) +
      geom_point(
        data = summary_data,
        aes(x = Timesteps, y = mean_val, fill = Nodes, shape = Nodes),
        position = dodge, size = 3, stroke = 0.5, color = "white"
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      facet_grid(Density ~ Regimes,
                 labeller = labeller(
                   Density = density_labeller,
                   Regimes = regime_labeller
                 )) +
      scale_color_thesis_nodes() +
      scale_fill_thesis_nodes() +
      scale_shape_manual(values = NODE_SHAPES, labels = NODE_LABELS) +
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

    output_file <- plots_file(paste0(net, "_senspec_plot_with_labels.pdf"))

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

  ggsave(plots_file("coefficient_plot_full.pdf"),
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

  ggsave(plots_file("coefficient_plot_main_effects.pdf"),
         plot = p_coef_main, device = cairo_pdf, width = 8, height = 5, dpi = 600)

  message("Coefficient plots saved to ", PLOTS_DIR)
}
