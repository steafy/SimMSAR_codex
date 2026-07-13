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
THESIS_VIRIDIS_RANGE <- c(0.15, 0.88)

scale_color_thesis_d <- function(...) {
  scale_color_viridis_d(begin = THESIS_VIRIDIS_RANGE[1], end = THESIS_VIRIDIS_RANGE[2], ...)
}
scale_fill_thesis_d <- function(...) {
  scale_fill_viridis_d(begin = THESIS_VIRIDIS_RANGE[1], end = THESIS_VIRIDIS_RANGE[2], ...)
}

# Mako triple for the 3-level Nodes mapping: an ordered, perceptually-uniform,
# colourblind-safe ramp (deep blue = N4 -> mint = N8). Mako is deliberately a
# DIFFERENT viridis-family palette from the Regime ramp (viridis proper): the
# two ordered dimensions therefore never look identical across figures, even
# though they never co-occur as two colour legends in one figure (Nodes is the
# colour dimension only where Regimes is a facet, and vice versa). Change these
# three hues here to restyle Nodes everywhere (Figures 3, 4 and the line /
# sens-spec plots). Legend entries use the node counts ("N = 4"), mirroring the
# "M = 1" Regime convention.
NODES_COLORS <- c("4" = "#3B5799", "6" = "#359BAA", "8" = "#84D9B1")
NODE_LABELS  <- c("4" = "4", "6" = "6", "8" = "8")
SENSPEC_COLORS <- c(
  "4.Sensitivity" = "#3B5799",
  "6.Sensitivity" = "#359BAA",
  "8.Sensitivity" = "#84D9B1",
  "4.Specificity"   = "#C11754",
  "6.Specificity"   = "#8B1E1E",
  "8.Specificity"   = "#EF5640"
  )

# Correlation line plots ONLY (make_line_plots): their own, wider viridis Nodes
# range so the three heavily overlapping IQR bands separate clearly -- mako's
# near-neighbour blue-greens muddy together there. Kept separate from
# NODES_COLORS so Figures 3/4 stay on mako. Edit begin/end to retune.
LINE_NODES_COLORS <- setNames(
  viridisLite::viridis(3, begin = 0.20, end = 0.80),
  c("4", "6", "8")
)

scale_color_thesis_nodes <- function(...) {
  scale_color_manual(values = NODES_COLORS, labels = NODE_LABELS, ...)
}
scale_fill_thesis_nodes <- function(...) {
  scale_fill_manual(values = NODES_COLORS, labels = NODE_LABELS, ...)
}

# Redundant shape coding for Nodes, IN ADDITION to colour: keeps figures legible
# in greyscale printouts/photocopies -- the "don't rely on colour alone"
# recommendation. Fillable variants (circle/triangle/square) so points can be
# drawn as a coloured fill with a white stroke (shape 21/24/22).
NODE_SHAPES <- c("4" = 21, "6" = 24, "8" = 22)

# --- Outcome palette + labels (third house palette, alongside Regime/Nodes) ---
# The four recovery outcomes, renamed to the thesis's own matrix notation
# (A = temporal Beta, K = contemporaneous Kappa) rather than the internal
# column names. Defined ONCE here and reused in every figure's legends/strips/
# titles via OUTCOME_LABELS, so the naming can never drift between figures.
OUTCOME_LABELS <- c(
  Beta_corr            = "temporal (A)",
  Kappa_corr           = "contemporaneous (K)",
  Beta_ac_corr_pearson = "controllability",
  RQ4                  = "regime sequence"
)

# Qualitative quartet (the "Viridis-System" outcome palette: a warm plasma
# spread, cool viridis being reserved for the ordered Regime ramp) keyed by the
# DISPLAY outcome name. Used only by Figure 1 (the coefficient plot), the single
# figure where outcome itself is the colour dimension -- Regime and Nodes never
# appear there, so the hue reuse across figures carries no within-figure
# ambiguity. Wrapped in scale_*_thesis_outcome() below, mirroring
# scale_*_thesis_d() (Regime) and scale_*_thesis_nodes() (Nodes), so all three
# palettes are edited in one place.
OUTCOME_COLORS <- c(
  "temporal (A)"                  = "#2B0594",  # deep violet
  "contemporaneous (K)"                  = "#99159F",  # magenta
  "controllability" = "#E26561",  # coral
  "regime sequence"      = "#FDC527"   # gold
)

scale_color_thesis_outcome <- function(...) {
  scale_color_manual(values = OUTCOME_COLORS, ...)
}
scale_fill_thesis_outcome <- function(...) {
  scale_fill_manual(values = OUTCOME_COLORS, ...)
}

# --- Shared point/line sizing ------------------------------------------------
# The line/prediction figures (2 and 4) share ONE point size + linewidth so
# their markers read as uniform across figures. The coefficient forest plot
# (Figure 1) is a different chart type -- dot-and-whisker, no data lines -- and
# its marker size is tuned independently (via THESIS_COEF_POINT_SIZE), so
# retuning the line-plot dots never silently resizes the forest-plot dots and
# vice versa. Figure 3 (the dumbbell) overrides both with larger points -- a
# lollipop chart needs prominent markers.
THESIS_POINT_SIZE      <- 1.7   # line/prediction figures (2, 3, 4): uniform dots
THESIS_LINEWIDTH       <- 0.9   # line/prediction figures (2, 3, 4)
THESIS_COEF_POINT_SIZE <- 3     # coefficient forest plot (1), controlled apart

# --- Opacity knobs (all independent) -----------------------------------------
# Figures 2-4 get SEPARATE dot vs. line opacity so each can be tuned on its own:
# translucent lines let overlapping prediction curves show through each other,
# while the data-point markers can stay crisp (or be softened) independently.
# Figure 1 (the forest plot) has its own knobs for its dots and its error bars,
# and both coefficient variants (compact + full) read the SAME error-bar
# linewidth/opacity so they stay visually identical.
THESIS_LINE_ALPHA  <- 0.60  # figures 2-4: prediction lines
THESIS_POINT_ALPHA <- 0.80   # figures 2-4: data-point markers
THESIS_COEF_ALPHA  <- 0.60  # figure 1: forest-plot dots
THESIS_COEF_ERRORBAR_LW    <- 0.7   # figure 1: error-bar line width (both variants)
THESIS_COEF_ERRORBAR_ALPHA <- 0.60  # figure 1: error-bar opacity (both variants)

# White outline width of the shape-21 dots in Figures 3 and 4. This white ring
# is what creates the little gap between a dot and the line running through it:
# a bigger stroke = a wider gap, 0 = no ring and no gap. (Set the dots' outline
# colour to NA in those geom_point() calls instead if you want them fully
# borderless.)
THESIS_POINT_STROKE <- 0.3

# Clean raw model term names into thesis notation, shared by the coefficient
# plot's interaction facets AND get_averaged_main_effects(), so the two never
# disagree on how a term is spelled: "a:b" -> "a x b", drop the _s/_num scaling
# suffixes ("Density_s" -> "Density"). A trailing "*" (RQ4 Regimes-vs-M=2 marker)
# is left untouched.
clean_term <- function(x) {
  x <- gsub(":", " × ", x)
  x <- gsub("_s\\b", "", x)
  x <- gsub("_num\\b", "", x)
  x
}

# --- Facet-strip labels ------------------------------------------------------
# Density on rows (right-hand strips), labelled directly as low/mid/high so no
# separate super-label is needed. Regimes on columns as "M = 1" ... "M = 4".
DENSITY_LABELS <- c("0.25" = "low density",
                    "0.5"  = "mid density",
                    "0.75" = "high density")

density_labeller <- function(x) {
  out <- DENSITY_LABELS[as.character(x)]
  ifelse(is.na(out), paste0("density: ", x), out)
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
      # Legend: pinned to the TOP and right-JUSTIFIED (not centered), with any
      # multiple legends stacked vertically (legend.box = "vertical"). Figures
      # carrying a second legend (e.g. Figure 1's significance shape) render it
      # BELOW the main colour legend. legend.justification can behave differently
      # across ggplot2 versions, so this is verified by eye on the rendered plot.
      legend.position    = "top",
      legend.justification = "right",
      legend.box         = "vertical",
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
      # More breathing room between axis titles and the panel. Set on the theme
      # itself (not per-figure) so every current and future figure inherits it.
      axis.title         = element_text(family = THESIS_FONT, face = "bold"),
      axis.title.x       = element_text(family = THESIS_FONT, face = "bold", margin = margin(t = 10)),
      axis.title.y       = element_text(family = THESIS_FONT, face = "bold", margin = margin(r = 10)),
      axis.text.x        = element_text(angle = 45, hjust = 1),
      # A little extra breathing room on the right of every figure (keeps
      # content/strips off the page edge and gives the top-right legend room).
      plot.margin        = margin(t = 6, r = 22, b = 6, l = 6)
    )
}

# -----------------------------------------------------------------------------
# PART 7: Line plots with facets (one per correlation outcome)
# -----------------------------------------------------------------------------
make_line_plots <- function(corr_results, cols) {

  line_plots <- list()

  for (col_name in cols) {

    title <- switch(
      col_name,
      "Beta_corr"           = "Estimation accuracy of temporal network (A)",
      "Kappa_corr"           = "Estimation accuracy of contemporaneous Netzwerks (K)",
      "Beta_ac_corr_pearson"  = "Estimation accuracy of average controllability (AC)",
      "Beta_ac_corr_spearman" = "Estimation accuracy of average controllability (AC) (Spearman)"
    )

    # Median line + empirical IQR (25th-75th percentile) ribbon per Nodes group,
    # computed directly on the raw per-regime correlations. BOUNDED BY
    # CONSTRUCTION: the ribbon edges are quantiles of the data, so they can never
    # exceed r = 1 the way a symmetric mean +- SD interval does, and the
    # asymmetry near the ceiling is shown faithfully. Non-parametric, and immune
    # to the near-|1| sensitivity of a Fisher-z mean (atanh(r) -> Inf, which can
    # drag a back-transformed mean far above where most of the data sit): a few
    # extreme replicates can only shift which side of the median they fall on,
    # never pull it toward them. (A lighter 10th-90th band can be added the same
    # way if a wider spread is wanted.)
    summ <- corr_results %>%
      mutate(T_num = as.numeric(as.character(Timesteps))) %>%
      group_by(T_num, Density, Nodes, Regimes) %>%
      summarise(
        med = median(.data[[col_name]], na.rm = TRUE),
        q25 = quantile(.data[[col_name]], 0.25, na.rm = TRUE),
        q75 = quantile(.data[[col_name]], 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    t_breaks <- sort(unique(summ$T_num))

    p <- ggplot(summ, aes(x = T_num, color = Nodes, fill = Nodes, group = Nodes)) +
      # Translucent IQR fill with a very thin same-colour outline: the edge stays
      # legible where the three Nodes bands overlap, without darkening the fill.
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.1, linewidth = 0.1, color = NA) +
      geom_line(aes(y = q25), linewidth = 0.1, alpha = 0.2) +
      geom_line(aes(y = q75), linewidth = 0.1, alpha = 0.2) +
      geom_line(aes(y = med), linewidth = THESIS_LINEWIDTH, alpha = THESIS_LINE_ALPHA,
                lineend = "round", linejoin = "round") +
      geom_point(aes(y = med), shape = 21, size = THESIS_POINT_SIZE,
                 alpha = THESIS_POINT_ALPHA, color = "white", stroke = THESIS_POINT_STROKE) +
      coord_cartesian(ylim = c(NA, 1)) +
      facet_grid(Regimes ~ Density,
                 labeller = labeller(
                   Density = density_labeller,
                   Regimes = regime_labeller
                 )) +
      # Log-scaled T-axis, matching the prediction figure (design is log-spaced).
      scale_x_log10(breaks = t_breaks, labels = t_breaks) +
      # Line-plot-specific viridis Nodes palette (LINE_NODES_COLORS): wider range
      # than the shared mako, to keep the overlapping IQR bands distinguishable.
      scale_color_manual(values = LINE_NODES_COLORS, labels = NODE_LABELS) +
      scale_fill_manual(values = LINE_NODES_COLORS, labels = NODE_LABELS) +
      labs(
        title    = title,
        x        = "time steps (T, log-scale)",
        y        = "correlations",
        color    = "nodes",
        fill     = "nodes"
      ) +
      theme_thesis(base_size = 12) +
      theme(
        plot.margin = margin(t = 2, r = 4, b = 0, l = 4),
        panel.spacing = unit(0.15, "in")
      ) +
    theme(aspect.ratio = 1)

    # Output filenames use the thesis's own matrix notation (A/K/AC) instead of
    # the internal column names, matching OUTCOME_LABELS above and the senspec
    # plot filenames below.
    file_prefix <- switch(
      col_name,
      "Beta_corr"              = "A_corr",
      "Kappa_corr"             = "K_corr",
      "Beta_ac_corr_pearson"   = "AC_corr",
      "Beta_ac_corr_spearman"  = "AC_corr_spearman",
      col_name
    )
    output_file <- plots_file(paste0(file_prefix, "_plot.pdf"))

    # cairo_pdf instead of the default pdf() device: the base device's font-
    # metric calculation assumes a single-byte locale and mangles multi-byte
    # UTF-8 characters (umlauts, sharp s, multiplication sign, plus-minus,
    # ...) when R is running under a C/POSIX locale (as in this sandbox) --
    # cairo renders the glyphs directly and isn't affected.
    ggsave(
      filename = output_file,
      plot     = p,
      device   = cairo_pdf,
      width    = 7,
      height   = 10,
      units    = "in",
      dpi      = 1200
    )

    message("saved: ", output_file)
    line_plots[[col_name]] <- p
  }

  invisible(line_plots)
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
        median_val = median(Value, na.rm = TRUE),
        q25        = quantile(Value, 0.25, na.rm = TRUE),
        q75        = quantile(Value, 0.75, na.rm = TRUE),
        .groups    = "drop"
      )
    
    dodge <- position_dodge(width = 0.5)

    p <- ggplot() +
      geom_errorbar(
        data = summary_data,
        aes(x = Timesteps, ymin = q25, ymax = q75,
            color = interaction(Nodes, Metric), group = interaction(Nodes, Metric)),
        width = 0.3, linewidth = 0.4, position = dodge
      ) +
      geom_line(
        data = summary_data,
        aes(x = Timesteps, y = median_val,
            color = interaction(Nodes, Metric), group = interaction(Nodes, Metric),
            linetype = Metric),
        position = dodge, lineend = "round", linejoin = "round"
      ) +
      geom_point(
        data = summary_data,
        aes(x = Timesteps, y = median_val, fill = interaction(Nodes, Metric),
            shape = Nodes, group = interaction(Nodes, Metric)),
        position = dodge, size = 1.8, stroke = 0.4, color = "white"
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      facet_grid(Regimes ~ Density,
                 labeller = labeller(
                   Density = density_labeller,
                   Regimes = regime_labeller
                 )) +
      scale_color_manual(values = SENSPEC_COLORS, guide = "none") +
      scale_fill_manual(values = SENSPEC_COLORS, guide = "none") +
      scale_shape_manual(values = NODE_SHAPES, labels = NODE_LABELS) +
      guides(
        linetype = guide_legend(override.aes = list(
          color = c(unname(SENSPEC_COLORS["6.Sensitivity"]),
                    unname(SENSPEC_COLORS["6.Specificity"]))
        )),
        shape = guide_legend(override.aes = list(fill = "grey35", color = "grey35"))
      ) +
      labs(
        title    = paste0("Sensitivity & Specificity: ",
                          if (net == "Beta") "temporal network (A)" else "contemporaneous network (K)"),
        x        = "time steps (T)",
        y        = "median value",
        color    = "nodes",
        shape    = "nodes",
        linetype = "metric"
      ) +
      theme_thesis(base_size = 12) +
      theme(
        plot.margin = margin(t = 2, r = 4, b = 0, l = 4),
        panel.spacing = unit(0.15, "in"),
        legend.box = "horizontal",
        legend.box.just = "right"
      ) +
      theme(aspect.ratio = 1)
    

    # Same A/K matrix-notation prefix as the line-plot filenames above.
    net_prefix <- switch(net, "Beta" = "A", "Kappa" = "K")
    output_file <- plots_file(paste0(net_prefix, "_senspec_plot.pdf"))

    ggsave(
      filename = output_file,
      plot     = p,
      device   = cairo_pdf,
      width    = 7,
      height   = 10,
      units    = "in",
      dpi      = 1200
    )

    message("saved: ", output_file)
  }
}

# -----------------------------------------------------------------------------
# PART 8.0: Averaged main effects (shared by BOTH coefficient plots)
# -----------------------------------------------------------------------------
# The three design covariates (logT, Density_s, Nodes_s) each interact with
# Regimes in the primary RQ1-3 models, so a raw summary(model) coefficient for
# them is the slope CONDITIONAL on Regimes = 1 (the reference level) -- not a
# fair "main effect". Here they are instead emtrends-averaged across Regimes
# with weights = "equal": Regimes is a fully-crossed, deliberately balanced
# design factor, so equal-weighting across its levels is correct (not the
# proportional/"proper" weighting emmeans would use by default, which would let
# whichever Regimes level happens to have more surviving fits dominate).
#
# The Regimes dummies themselves are taken as RAW coefficients: Density_s /
# Nodes_s / logT are centered, so each Regimes contrast is already evaluated at
# typical conditions and needs no further averaging.
#
# RQ4 (Regimesequenz) is added as a 4th outcome: logT interacts with Regimes
# there too (-> emtrends-averaged), but Density_s / Nodes_s do NOT (-> raw, they
# are already unconditional main effects), and its Regimes contrasts are
# relative to M = 2 (M = 1 is excluded from RQ4 entirely), so they are marked
# with a trailing "*" and flagged in the plot caption.
#
# Computed ONCE and returned tidy (outcome, term, estimate, se, ci_lo, ci_hi,
# p, sig) so the compact main-effects plot and the full plot's "Haupteffekte"
# facet are fed from the exact same numbers and can never silently disagree.

# Shared top-to-bottom / row order for the six (RQ1-3) or five (RQ4) main-
# effect terms. Used by make_coefficient_plots() (y-axis order) AND
# export_main_effects_table() (row order) so the plot and its companion
# table can never silently diverge in ordering.
MAIN_EFFECTS_TERM_ORDER <- c("logT", "Density", "Nodes",
                             "Regimes2", "Regimes3", "Regimes4",
                             "Regimes3*", "Regimes4*")

get_averaged_main_effects <- function(all_results, seq_recovery, corr_cols) {

  # Asymptotic (z-based) inference: matches the 1.96*SE CIs used everywhere else
  # in these plots, and avoids the very slow Kenward-Roger df on a dataset this
  # large. Set locally each call so sourcing order can't leave it unset.
  emmeans::emm_options(lmer.df = "asymptotic")

  # emtrends-averaged slope of numeric predictor `v`, marginalized over Regimes
  # with equal weights. `~ 1` collapses to a single overall trend.
  avg_slope <- function(model, v, data) {
    emt <- emmeans::emtrends(model, ~ 1, var = v, weights = "equal", data = data)
    s   <- as.data.frame(summary(emt, infer = c(TRUE, TRUE)))
    data.frame(
      term      = v,
      estimate  = s[[paste0(v, ".trend")]],
      se        = s[["SE"]],
      df        = s[["df"]],          # wird bei asymptotic = Inf sein
      statistic = s[["z.ratio"]],     # heißt bei df=Inf "z.ratio", nicht "t.ratio"
      ci_lo     = s[["asymp.LCL"]],
      ci_hi     = s[["asymp.UCL"]],
      p         = s[["p.value"]],
      stringsAsFactors = FALSE
    )
  }
  
  # Raw model coefficient for one term (95% CI on the normal approximation).
  raw_coef <- function(cs, term, out_term = term) {
    est <- cs[term, "Estimate"]; se <- cs[term, "Std. Error"]
    data.frame(
      term      = out_term,
      estimate  = est,
      se        = se,
      df        = cs[term, "df"],
      statistic = cs[term, "t value"],
      ci_lo     = est - 1.96 * se,
      ci_hi     = est + 1.96 * se,
      p         = cs[term, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }

  out <- list()

  # --- RQ1-3: all three slopes emtrends-averaged; Regimes dummies raw ---------
  for (outcome in corr_cols) {
    model <- all_results[[outcome]]$primary_model
    d     <- model.frame(model)
    cs    <- summary(model)$coefficients
    rows  <- dplyr::bind_rows(
      avg_slope(model, "logT",      d),
      avg_slope(model, "Density_s", d),
      avg_slope(model, "Nodes_s",   d),
      raw_coef(cs, "Regimes2"),
      raw_coef(cs, "Regimes3"),
      raw_coef(cs, "Regimes4")
    )
    rows$outcome <- OUTCOME_LABELS[[outcome]]
    out[[outcome]] <- rows
  }

  # --- RQ4: logT emtrends-averaged; Density/Nodes raw; Regimes* rel. to M=2 ---
  if (!is.null(seq_recovery) && !is.null(seq_recovery$primary_model)) {
    mr  <- seq_recovery$primary_model
    d4  <- model.frame(mr)
    cs4 <- summary(mr)$coefficients
    rows4 <- dplyr::bind_rows(
      avg_slope(mr, "logT", d4),
      raw_coef(cs4, "Density_s"),
      raw_coef(cs4, "Nodes_s"),
      raw_coef(cs4, "Regimes3", "Regimes3*"),
      raw_coef(cs4, "Regimes4", "Regimes4*")
    )
    rows4$outcome <- OUTCOME_LABELS[["RQ4"]]
    out[["RQ4"]] <- rows4
  }

  res <- dplyr::bind_rows(out)

  # Significance via BH within each outcome's own main-effect set (mirrors the
  # per-model BH the interaction facets keep using). Then clean term names to
  # thesis notation via the shared clean_term() ("Density_s" -> "Density").
  res <- res %>%
    group_by(outcome) %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           sig   = ifelse(p_adj < 0.05, "p < .05", "n.s.")) %>%
    ungroup()
  res$term    <- clean_term(res$term)
  res$outcome <- factor(res$outcome, levels = unname(OUTCOME_LABELS))
  res
}

# -----------------------------------------------------------------------------
# PART 8: Coefficient plot (full + main-effects-only) -- FIGURE 1
# -----------------------------------------------------------------------------
make_coefficient_plots <- function(all_results, corr_cols, seq_recovery = NULL) {

  # Top-to-bottom order of main-effect terms on the y-axis (RQ4's starred
  # Regimes contrasts sit just below the shared M=1-referenced ones).
  main_order <- MAIN_EFFECTS_TERM_ORDER

  # --- Averaged main effects (shared source for both plots) ------------------
  main_eff <- get_averaged_main_effects(all_results, seq_recovery, corr_cols) %>%
    mutate(n_cross = 0L)

  # --- Interaction terms: RAW coefficients, RQ1-3 only, per-model BH (as before)
  # The 2-/3-way facets keep the original behaviour verbatim: interactions ARE
  # the statement of how an effect changes across Regimes, so the raw coefficient
  # is already the correct number -- no Regimes-averaging applies. BH is over the
  # full per-model coefficient set, exactly as the previous version computed it,
  # so these facets are unchanged.
  inter <- lapply(corr_cols, function(outcome) {
    cs <- summary(all_results[[outcome]]$primary_model)$coefficients
    data.frame(
      outcome  = OUTCOME_LABELS[[outcome]],
      term     = rownames(cs),
      estimate = cs[, "Estimate"],
      se       = cs[, "Std. Error"],
      p_adj    = p.adjust(cs[, "Pr(>|t|)"], method = "BH"),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  inter$term <- clean_term(inter$term)
  inter <- inter %>%
    filter(term != "(Intercept)") %>%
    mutate(n_cross = stringr::str_count(term, "\u00d7")) %>%
    filter(n_cross >= 1) %>%
    mutate(
      ci_lo = estimate - 1.96 * se,
      ci_hi = estimate + 1.96 * se,
      sig   = ifelse(p_adj < 0.05, "p < .05", "n.s."),
      outcome = factor(outcome, levels = unname(OUTCOME_LABELS))
    )

  keep <- c("outcome", "term", "estimate", "ci_lo", "ci_hi", "sig", "n_cross")
  full_data <- bind_rows(main_eff[keep], inter[keep])

  # y-axis factor: main effects in explicit order, then interactions by
  # complexity then alphabetically (facets are free_y, but one global factor
  # keeps ordering deterministic).
  me_terms  <- main_order[main_order %in% full_data$term]
  int_terms <- full_data %>% filter(n_cross >= 1) %>%
    arrange(n_cross, term) %>% pull(term) %>% unique()
  term_levels <- c(me_terms, int_terms)
  full_data$term <- factor(full_data$term, levels = rev(term_levels))

  # Colour (Outcome) legend on top, significance shape legend BELOW it.
  coef_guides <- guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

  # --- 8.1: Full coefficient plot (all terms) --------------------------------
  p_coef_full <- ggplot(full_data,
                        aes(x = estimate, y = term,
                            color = outcome, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   linewidth = THESIS_COEF_ERRORBAR_LW, alpha = THESIS_COEF_ERRORBAR_ALPHA,
                   height = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(size = THESIS_COEF_POINT_SIZE, alpha = THESIS_COEF_ALPHA,
               position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1), guide = "none") +
    scale_color_thesis_outcome() +
    facet_wrap(~ n_cross, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(
                 "0" = "Main effects (averaged over regimes)",
                 "1" = "2-way interactions",
                 "2" = "3-way interactions"
               ))) +
    labs(
      title    = "Regression coefficients primary model",
      x        = "Regression coefficients (95%-ci)",
      y        = NULL,
      color    = "Outcome",
      shape    = "significant"
    ) +
    coef_guides +
    guides(shape = "none") +
    theme_thesis(base_size = 15) +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(t = 2, r = 7, b = 0, l = 4)
    )

  ggsave(plots_file("coefficient_plot_full.pdf"),
         plot = p_coef_full, device = cairo_pdf, width = 10, height = 14, dpi = 600)

  # --- 8.2: Main effects only (compact version for main text) ----------------
  main_plot_data <- main_eff %>%
    mutate(term = factor(term, levels = rev(me_terms)))

  p_coef_main <- ggplot(main_plot_data,
                        aes(x = estimate, y = term,
                            color = outcome, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   linewidth = THESIS_COEF_ERRORBAR_LW, alpha = THESIS_COEF_ERRORBAR_ALPHA,
                   height = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(size = THESIS_COEF_POINT_SIZE, alpha = THESIS_COEF_ALPHA,
               position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1), guide = "none") +
    scale_color_thesis_outcome() +
    labs(
      title    = "Main effects (averaged over regimes)",
      x        = "Regression coefficients (95%-ci)",
      y        = NULL,
      color    = "Outcome",
      shape    = "significant"
    ) +
    coef_guides +
    guides(shape = "none") +
    theme_thesis(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(plots_file("coefficient_plot_main_effects.pdf"),
         plot = p_coef_main, device = cairo_pdf, width = 8, height = 5.5, dpi = 800)

  message("Coefficient plots saved to ", PLOTS_DIR)

  invisible(list(full = p_coef_full, main = p_coef_main))
}

# -----------------------------------------------------------------------------
# Shared prediction helpers for the model-prediction figures (2, 3, 4)
# -----------------------------------------------------------------------------
# The primary models are fit on SCALED predictors (logT = scale(log(T)),
# Density_s = scale(Density), Nodes_s = scale(Nodes)); the model frame keeps only
# the scaled columns, not the raw design values. To predict at interpretable raw
# levels (T = 200..1600, Density = 0.25/0.5/0.75, Nodes = 4/8) we recover the
# exact linear map raw -> scaled from a data frame that still carries both. The
# map is exactly linear (an affine centre/scale), so a 2-point lm recovers it
# with no approximation.
.lin_map <- function(from, to) {
  d   <- unique(data.frame(from = as.numeric(from), to = as.numeric(to)))
  fit <- stats::lm(to ~ from, data = d)
  co  <- stats::coef(fit)
  function(x) as.numeric(co[[1]] + co[[2]] * x)
}

# Population-level (re.form = NA) prediction over a grid, holding the untouched
# covariates at their scaled mean of 0. `fixed` supplies the scaled predictor
# columns already on the grid; Regimes is coerced to the model's own factor
# levels. tanh back-transform is applied for the Fisher-z correlation outcomes
# and skipped for RQ4's Cohen's kappa (fit on the raw kappa scale).
.predict_recovery <- function(model, grid, tanh_bt) {
  lv <- levels(model.frame(model)$Regimes)
  grid$Regimes <- factor(as.character(grid$Regimes), levels = lv)
  pr <- predict(model, newdata = grid, re.form = NA)
  grid$value <- if (tanh_bt) tanh(pr) else pr
  grid
}

# The four T design points, shared by Figures 2 (and the log x-axis breaks).
THESIS_T_POINTS <- c(200, 400, 800, 1600)

# -----------------------------------------------------------------------------
# FIGURE 2: predicted recovery vs. T (one facet per outcome, line per Regime)
# -----------------------------------------------------------------------------
# For each outcome's PRIMARY model, the predicted recovery is drawn across a
# continuous logT grid (T = 200..1600) at Density_s = Nodes_s = 0, one line per
# Regime level, with points marking the four actual T design values. The x-axis
# is log-scaled: the model is fit on log(T) and the design is log-spaced (each T
# level doubles), so a linear T-axis would visually exaggerate the early curve.
make_prediction_curves <- function(all_results, seq_recovery, dat_sim) {

  map_logT_corr <- .lin_map(log(dat_sim$Timesteps), dat_sim$logT)
  sr            <- seq_recovery$seq_results
  map_logT_rq4  <- .lin_map(log(sr$timesteps), sr$logT)

  T_line <- exp(seq(log(min(THESIS_T_POINTS)), log(max(THESIS_T_POINTS)), length.out = 200))

  # Build line + point data for one outcome/model.
  build_one <- function(model, label, tanh_bt, map_logT) {
    lv <- levels(model.frame(model)$Regimes)
    mk <- function(Tvals) {
      g <- expand.grid(Timesteps = Tvals, Regimes = lv, stringsAsFactors = FALSE)
      g$logT <- map_logT(log(g$Timesteps)); g$Density_s <- 0; g$Nodes_s <- 0
      out <- .predict_recovery(model, g, tanh_bt)
      out$outcome <- label
      out
    }
    list(line = mk(T_line), point = mk(THESIS_T_POINTS))
  }

  specs <- c(
    lapply(corr_cols, function(oc) list(
      model = all_results[[oc]]$primary_model, label = OUTCOME_LABELS[[oc]],
      tanh_bt = TRUE, map = map_logT_corr)),
    list(list(model = seq_recovery$primary_model, label = OUTCOME_LABELS[["RQ4"]],
              tanh_bt = FALSE, map = map_logT_rq4))
  )

  built     <- lapply(specs, function(s) build_one(s$model, s$label, s$tanh_bt, s$map))
  line_df   <- bind_rows(lapply(built, `[[`, "line"))
  point_df  <- bind_rows(lapply(built, `[[`, "point"))

  # Consistent Regime colouring/ordering and outcome facet order across figures.
  lvl_reg <- c("1", "2", "3", "4")
  line_df$Regimes  <- factor(as.character(line_df$Regimes),  levels = lvl_reg)
  point_df$Regimes <- factor(as.character(point_df$Regimes), levels = lvl_reg)
  line_df$outcome  <- factor(line_df$outcome,  levels = unname(OUTCOME_LABELS))
  point_df$outcome <- factor(point_df$outcome, levels = unname(OUTCOME_LABELS))

  p <- ggplot(line_df, aes(x = Timesteps, y = value, color = Regimes)) +
    geom_line(aes(group = Regimes), linewidth = THESIS_LINEWIDTH,
              alpha = THESIS_LINE_ALPHA, lineend = "round", linejoin = "round") +
    # White-outlined shape-21 dots (fill = Regime), matching Figures 3 and 4.
    geom_point(data = point_df, aes(fill = Regimes), shape = 21,
               size = THESIS_POINT_SIZE, alpha = THESIS_POINT_ALPHA,
               color = "white", stroke = THESIS_POINT_STROKE) +
    # SHARED (fixed) y-axis across all four outcomes so the curves are directly
    # comparable panel-to-panel -- an equal predicted recovery sits at the same
    # height everywhere, rather than each facet zooming to its own range.
    facet_wrap(~ outcome,  nrow = 1, ncol = 4) +
    scale_x_log10(breaks = THESIS_T_POINTS, labels = THESIS_T_POINTS) +
    scale_color_thesis_d(drop = FALSE, limits = lvl_reg, labels = regime_labeller) +
    scale_fill_thesis_d(drop = FALSE, limits = lvl_reg, labels = regime_labeller) +
    labs(
      title    = "Interaction of regimes and time steps",
      x        = "time steps (T, log-scaled)",
      y        = "predicted outcome (r / κ)",
      color    = "regime",
      fill     = "regime"
    ) +
    theme_thesis(base_size = 12) + 
    theme(aspect.ratio = 1)
    

  ggsave(plots_file("regimes_x_timesteps.pdf"),
         plot = p, device = cairo_pdf, width = 9, height = 4, dpi = 800)
  message("saved: ", plots_file("regimes_x_timesteps.pdf"))

  invisible(p)
}

# -----------------------------------------------------------------------------
# FIGURE 3: Nodes effect by Regime -- interaction profile (A and K only)
# -----------------------------------------------------------------------------
# For A (Beta_corr) and K (Kappa_corr), the predicted recovery at Nodes = 4 vs.
# Nodes = 8 (holding logT = 0, Density_s = 0) is drawn as two profile lines
# across the regimes with the between-Nodes gap shaded, on TWO scales side by
# side in one 1x4 layout:
#   - LEFT half  (z-scale, tanh_bt = FALSE): the model's native, linear scale.
#     This is where the Nodes x Regimes INTERACTION COEFFICIENTS actually live,
#     so it shows the genuine modelled effect -- for K, the gap narrows as M
#     rises; for A, it stays negligible throughout.
#   - RIGHT half (r-scale, tanh_bt = TRUE): the practically interpretable
#     correlation scale used elsewhere in the results (e.g. Figure 2). Because
#     tanh() compresses differences near the ceiling and stretches them further
#     away from it, the SAME (shrinking) z-scale gap can translate into a
#     WIDER-looking r-scale gap at high M for K -- a ceiling-compression
#     artifact, not a second, independent effect. Showing both scales side by
#     side keeps that distinction visible instead of implicit.
# An italic "delta = N4 - N8" label quantifies the gap at each regime, on
# whichever scale that panel uses.
make_nodes_regime_profile <- function(all_results, dat_sim) {
  
  map_nodes   <- .lin_map(dat_sim$Nodes, dat_sim$Nodes_s)
  node_levels <- c("4", "8")
  
  # Predictions on both scales, from the same predictor grid (Regimes x Nodes,
  # logT = Density_s = 0) -- only the tanh back-transform differs.
  build_pts <- function(tanh_bt) {
    lapply(c("Beta_corr", "Kappa_corr"), function(oc) {
      model <- all_results[[oc]]$primary_model
      lv    <- levels(model.frame(model)$Regimes)
      g <- expand.grid(Regimes = lv, Nodes = c(4, 8), stringsAsFactors = FALSE)
      g$logT <- 0; g$Density_s <- 0; g$Nodes_s <- map_nodes(g$Nodes)
      out <- .predict_recovery(model, g, tanh_bt = tanh_bt)
      out$outcome <- OUTCOME_LABELS[[oc]]
      out$Nodes   <- factor(as.character(out$Nodes), levels = node_levels)
      out
    }) %>% bind_rows()
  }
  
  pts_z <- build_pts(tanh_bt = FALSE)  # Fisher-z (model-native) scale
  pts_r <- build_pts(tanh_bt = TRUE)   # r (tanh-back-transformed) scale
  
  # Shared panel-builder: one scale's worth of data in, one facetted (A, K)
  # ggplot out. Colour pair, y-axis label, subtitle and delta-label prefix are
  # all parameterised so the z- and r-scale calls below are the only place the
  # two halves differ.
  build_panel <- function(pts, line_cols, gap_fill, y_lab, subtitle, delta_prefix) {
    
    pts$outcome <- factor(pts$outcome, levels = unname(OUTCOME_LABELS))
    pts$Mnum    <- as.integer(as.character(pts$Regimes))
    
    n4 <- pts %>% filter(Nodes == "4") %>% select(outcome, Mnum, N4 = value)
    n8 <- pts %>% filter(Nodes == "8") %>% select(outcome, Mnum, N8 = value)
    wide <- left_join(n4, n8, by = c("outcome", "Mnum")) %>%
      mutate(lo = pmin(N4, N8), hi = pmax(N4, N8),
             # Clamp near-zero deltas to exactly 0 (avoids a stray "-0.00").
             delta  = ifelse(abs(N4 - N8) < 0.005, 0, N4 - N8),
             dlabel = sprintf("%s = %.2f", delta_prefix, delta))
    
    reg_breaks <- sort(unique(pts$Mnum))
    
    ggplot(pts, aes(x = Mnum, y = value)) +
      geom_ribbon(data = wide, inherit.aes = FALSE,
                  aes(x = Mnum, ymin = lo, ymax = hi),
                  fill = gap_fill, alpha = 0.18) +
      geom_line(aes(color = Nodes, group = Nodes),
                linewidth = THESIS_LINEWIDTH, alpha = THESIS_LINE_ALPHA,
                lineend = "round") +
      geom_point(aes(fill = Nodes), shape = 21, alpha = THESIS_POINT_ALPHA,
                 size = THESIS_POINT_SIZE, color = "white", stroke = THESIS_POINT_STROKE) +
      # geom_text(data = wide, inherit.aes = FALSE,
      #           aes(x = Mnum, y = hi, label = dlabel),
      #           vjust = -1, fontface = "italic", size = 5, color = "grey30") +
      # Fixed (shared) y-axis WITHIN this scale's two facets (A, K) -- ggplot's
      # facet_wrap default -- so the A vs. K gap sizes stay directly
      # comparable to each other on this scale. The z- and r-halves each get
      # their own such fixed axis via being separate ggplot objects below.
      facet_wrap(~ outcome, nrow = 1) +
      scale_x_continuous(breaks = reg_breaks, labels = regime_labeller(reg_breaks),
                         expand = expansion(mult = c(0.08, 0.08))) +
      scale_color_manual(values = line_cols, labels = NODE_LABELS[c("4", "8")]) +
      scale_fill_manual(values = line_cols, labels = NODE_LABELS[c("4", "8")]) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
      labs(
        x        = "regime",
        y        = y_lab,
        color    = "nodes",
        fill     = "nodes"
      ) +
      theme_thesis(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(margin = margin(r = 12)))
  }
  
  # r-scale colours: UNCHANGED from the original single-scale figure (mako
  # blue/teal, via the shared NODES_COLORS palette: "4" and "6" hues).
  r_line_cols <- c("4" = unname(NODES_COLORS[["4"]]),
                   "8" = unname(NODES_COLORS[["6"]]))
  r_gap_fill  <- unname(NODES_COLORS[["4"]])
  
  # z-scale colours: a separate pink/rose pair, deliberately NOT reusing any
  # hue from NODES_COLORS or OUTCOME_COLORS, so the two scales are visually
  # distinguishable at a glance even before reading axis labels, with no
  # accidental colour overlap with Nodes elsewhere or the outcome hues in
  # Figure 1.
  z_line_cols <- c("4" = "#B3084C", "8" = "#F2789F")
  z_gap_fill  <- "#B3084C"
  
  p_z <- build_panel(pts_z, z_line_cols, z_gap_fill,
                     y_lab        = "pred. corr. (Fisher-z)",
                     delta_prefix = "\u0394z")
  
  p_r <- build_panel(pts_r, r_line_cols, r_gap_fill,
                     y_lab        = "pred. corr. (r-Scale)",
                     delta_prefix = "\u0394r")
  
  # Shared title spanning the whole combined figure. A labs(title=) on p_z or
  # p_r would only title THAT half -- cowplot has no native concept of one
  # title over several combined plots, so it's built as its own tiny "plot"
  # (draw_label on an empty canvas) and stacked above the z/r row.
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Interaction of nodes and regimes",
      fontface   = "bold",
      size       = 24,
      fontfamily = THESIS_FONT,
      x = 0.5, hjust = 0.5
    )
  
  plot_row <- cowplot::plot_grid(p_z, p_r, nrow = 1, rel_widths = c(1, 1))
  
  p <- cowplot::plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.06, 1))
  
  ggsave(plots_file("node_x_regimes.pdf"),
         plot = p, device = cairo_pdf, width = 15, height = 6.1, dpi = 1200)
  message("saved: ", plots_file("node_x_regimes.pdf"))
  
  invisible(p)
}

# -----------------------------------------------------------------------------
# FIGURE 4: Density x Nodes x Regime grid (A, K, Kontrollierbarkeit)
# -----------------------------------------------------------------------------
# facet_grid(outcome ~ Regimes): outcome as rows (strip on the right), all four
# Regime levels as columns (regime_labeller strips on top). Within each small
# panel, predicted recovery is drawn against Density (0.25/0.5/0.75), one
# line+points per Nodes level (4 vs 8), holding logT = 0. The three correlation
# outcomes are tanh-back-transformed. Nodes uses the same two colours as
# Figure 3 (N=4 / N=6 hues) so the two N=4-vs-N=8 comparisons stay consistent.
make_density_nodes_regime_grid <- function(all_results, dat_sim) {

  map_dens  <- .lin_map(dat_sim$Density, dat_sim$Density_s)
  map_nodes <- .lin_map(dat_sim$Nodes,   dat_sim$Nodes_s)
  dens_vals <- sort(unique(dat_sim$Density))
  outcomes  <- c("Beta_corr", "Kappa_corr", "Beta_ac_corr_pearson")

  dat <- lapply(outcomes, function(oc) {
    model <- all_results[[oc]]$primary_model
    lv    <- levels(model.frame(model)$Regimes)
    g <- expand.grid(Density = dens_vals, Nodes = c(4, 8), Regimes = lv,
                     stringsAsFactors = FALSE)
    g$logT <- 0; g$Density_s <- map_dens(g$Density); g$Nodes_s <- map_nodes(g$Nodes)
    out <- .predict_recovery(model, g, tanh_bt = TRUE)
    out$outcome <- OUTCOME_LABELS[[oc]]
    out$Nodes   <- factor(as.character(out$Nodes), levels = c("4", "8"))
    out
  }) %>% bind_rows()

  dat$outcome <- factor(dat$outcome,
                        levels = unname(OUTCOME_LABELS[c("Beta_corr", "Kappa_corr",
                                                         "Beta_ac_corr_pearson")]))
  dat$Regimes <- factor(as.character(dat$Regimes), levels = c("1", "2", "3", "4"))

  # Figure-4-specific Nodes colours (requested): the two warm hues from the
  # "Sonnenuntergang" Regime ramp -- N=4 takes that scheme's M3 colour, N=8 its
  # M2 colour.
  fig4_node_cols <- c("4" = "#EF5640",   # Sonnenuntergang M3
                      "8" = "#C11754")   # Sonnenuntergang M2

  p <- ggplot(dat, aes(x = Density, y = value, color = Nodes)) +
    geom_line(aes(group = Nodes), linewidth = THESIS_LINEWIDTH,
              alpha = THESIS_LINE_ALPHA, lineend = "round") +
    geom_point(aes(fill = Nodes), shape = 21, size = THESIS_POINT_SIZE,
               alpha = THESIS_POINT_ALPHA, color = "white", stroke = THESIS_POINT_STROKE) +
    # One shared y-axis across ALL panels (not free per row) so every panel is
    # directly comparable, at the cost of compressing the flatter outcomes.
    facet_grid(outcome ~ Regimes,
               labeller = labeller(Regimes = regime_labeller)) +
    scale_x_continuous(breaks = dens_vals,
                       labels = c("low", "mid", "high")) +
    scale_color_manual(values = fig4_node_cols, labels = NODE_LABELS[c("4", "8")]) +
    scale_fill_manual(values = fig4_node_cols, labels = NODE_LABELS[c("4", "8")]) +
    labs(
      title = "Interaction of density, nodes and regimes",
      x     = "density",
      y     = "predicted correlations",
      color = "nodes",
      fill  = "nodes"
    ) +
    theme_thesis(base_size = 15)

  ggsave(plots_file("density_x_nodes_x_regimes.pdf"),
         plot = p, device = cairo_pdf, width = 10, height = 7.5, dpi = 1200)
  message("saved: ", plots_file("density_x_nodes_x_regimes.pdf"))

  invisible(p)
}
