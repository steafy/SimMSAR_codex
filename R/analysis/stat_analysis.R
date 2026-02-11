# =============================================================================
# Statistical Analysis Script for MSAR Models (FINAL VERSION)
# =============================================================================
# -----------------------------------------------------------------------------
# SETUP: Load required packages
# -----------------------------------------------------------------------------

# Core packages
required_packages <- c("dplyr", "ggplot2", "kableExtra", "xtable", "lmtest", "sandwich")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("✓ All required packages loaded\n\n")

# -----------------------------------------------------------------------------
# PART 1: DATA PREPARATION
# -----------------------------------------------------------------------------

cat("=== PART 1: DATA PREPARATION ===\n\n")

# Verify we have an msar_results object
if (!inherits(MSAR_dynamics_list, "msar_results")) {
  stop("MSAR_dynamics_list is not an msar_results object. Please re-run estimate_MSAR().")
}

cat("✓ MSAR_dynamics_list verified\n")
cat("  Dimensions:", nrow(MSAR_dynamics_list), "rows ×", ncol(MSAR_dynamics_list), "columns\n\n")

# -----------------------------------------------------------------------------
# 1.1: Create corr_results with grouping variables
# -----------------------------------------------------------------------------

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
  # Create Condition identifier
  mutate(Condition = interaction(Timesteps, Density, Nodes, Regimes, drop = TRUE)) %>%
  # IMPORTANT: Create SimUID (globally unique simulation identifier)
  mutate(SimUID = interaction(Condition, SimID, drop = TRUE)) %>%
  # Reorder columns for clarity
  select(SimUID, Condition, SimID, Timesteps, Density, Nodes, Regimes, RegimeIndex, N, everything())

cat("✓ corr_results created\n")
cat("  Unique SimUID:", length(unique(corr_results$SimUID)), "\n")
cat("  Unique Conditions:", length(unique(corr_results$Condition)), "\n\n")

# Transform to factors
corr_results <- corr_results %>%
  mutate(
    Timesteps = factor(Timesteps),
    Density = factor(Density),
    Nodes = factor(Nodes),
    Regimes = factor(Regimes),
    RegimeIndex = factor(RegimeIndex)
  )

cat("✓ Variables converted to factors\n\n")

# -----------------------------------------------------------------------------
# 1.2: Calculate descriptive statistics (optional - for reporting)
# -----------------------------------------------------------------------------


descriptive_stats_agg <- data.frame(
  Metric = character(),
  Mean = numeric(),
  SD = numeric(),
  Median = numeric(),
  Min = numeric(),
  Max = numeric(),
  stringsAsFactors = FALSE
)

# Metrics to summarize (all outcome variables, not just correlations)
all_metrics <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")

# Add other metrics if they exist in your data
other_metrics <- c("MAE_Wtemp", "Wtemp_sen", "Wtemp_spec",
                   "MAE_Wcont", "Wcont_sen", "Wcont_spec")

# Check which metrics are available in corr_results
available_other_metrics <- other_metrics[other_metrics %in% names(corr_results)]

if (length(available_other_metrics) > 0) {
  all_metrics <- c(all_metrics, available_other_metrics)
}

# Aggregate OTHER metrics (not correlations, which are already aggregated as _z)
if (length(available_other_metrics) > 0) {
  
  cat("  Aggregating other metrics (MAE, sensitivity, specificity)...\n")
  
  # Aggregate by taking mean across regimes for each simulation
  dat_sim_other <- corr_results %>%
    group_by(SimUID, Condition, SimID, 
             Timesteps, Density, Nodes, Regimes) %>%
    summarise(
      across(all_of(available_other_metrics), mean, .names = "{.col}"),
      .groups = "drop"
    )
  
  # Merge with existing dat_sim
  dat_sim <- dat_sim %>%
    left_join(dat_sim_other, by = c("SimUID", "Condition", "SimID"))
  
  cat("  ✓ Other metrics aggregated\n\n")
}

# Now calculate descriptive stats from aggregated data
for (metric in all_metrics) {
  
  if (metric %in% names(dat_sim)) {
    # For already aggregated data
    values <- dat_sim[[metric]]
  } else if (paste0(metric, "_z") %in% names(dat_sim)) {
    # For z-transformed correlations, back-transform first
    values <- tanh(dat_sim[[paste0(metric, "_z")]])
  } else {
    cat("  ⚠ Metric", metric, "not found, skipping\n")
    next
  }
  
  descriptive_stats_agg <- rbind(
    descriptive_stats_agg,
    data.frame(
      Metric = metric,
      Mean = mean(values, na.rm = TRUE),
      SD = sd(values, na.rm = TRUE),
      Median = median(values, na.rm = TRUE),
      Min = min(values, na.rm = TRUE),
      Max = max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}






# Note: get_stats() is assumed to be available from your package
# If not available, this section can be skipped or replaced

if (exists("get_stats") && is.function(get_stats)) {
  cat("Calculating descriptive statistics...\n")
  descript_stats_data <- get_stats(MSAR_dynamics_list)
  
  # Convert to old format for compatibility
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
    
    descript_stats[[old_name]] <- descript_stats_data %>%
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
  
  cat("✓ Descriptive statistics calculated\n\n")
} else {
  cat("⚠ get_stats() not available - skipping descriptive statistics\n\n")
  descript_stats <- NULL
}


descriptive_stats_report <- as.data.frame(matrix(nrow = 5, ncol = 10,
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
  descriptive_stats_report[ ,i] <- val
}

# -----------------------------------------------------------------------------
# 1.3: Calculate estimation process statistics
# -----------------------------------------------------------------------------

cat("Calculating estimation process statistics...\n")

aggr_factors <- corr_results %>%
  select(Timesteps, Density, Nodes, Regimes, N) %>%
  distinct()

n_ts_per_condition <- n_distinct(corr_results$SimID)
n_conditions <- nrow(aggr_factors)
expected_total <- n_ts_per_condition * n_conditions

omissions <- expected_total - sum(aggr_factors$N)
omissions_per <- (omissions / expected_total) * 100

cat("✓ Estimation statistics:\n")
cat("  Expected total:", expected_total, "\n")
cat("  Actual total:", sum(aggr_factors$N), "\n")
cat("  Omissions:", omissions, "(", round(omissions_per, 2), "%)\n\n")

# -----------------------------------------------------------------------------
# PART 2: FISHER-Z TRANSFORMATION
# -----------------------------------------------------------------------------

cat("=== PART 2: FISHER-Z TRANSFORMATION ===\n\n")

corr_cols <- c("Wtemp_corr", "Wtemp_ac_corr", "Wcont_corr", "Wcont_ac_corr")
eps <- 1e-6

for (cc in corr_cols) {
  r <- corr_results[[cc]]
  # Clip extreme values to avoid Inf/-Inf
  r <- pmin(pmax(r, -1 + eps), 1 - eps)
  # Fisher-z transformation: z = atanh(r)
  corr_results[[paste0(cc, "_z")]] <- atanh(r)
}

cat("✓ Fisher-z transformation completed\n")
cat("  Example - Wtemp_corr:\n")
cat("    Original range:", round(range(corr_results$Wtemp_corr), 3), "\n")
cat("    Fisher-z range:", round(range(corr_results$Wtemp_corr_z), 3), "\n")

# Check for extreme values
z_cols <- paste0(corr_cols, "_z")
for (zc in z_cols) {
  n_extreme <- sum(abs(corr_results[[zc]]) > 5)
  if (n_extreme > 0) {
    cat("  ⚠", zc, "has", n_extreme, "extreme values (|z| > 5)\n")
  }
}
cat("\n")

# -----------------------------------------------------------------------------
# PART 3: AGGREGATION TO SIMULATION-LEVEL
# -----------------------------------------------------------------------------

cat("=== PART 3: AGGREGATION TO SIMULATION-LEVEL ===\n\n")

# Convert factors to numeric for aggregation
corr_results_for_agg <- corr_results %>%
  mutate(
    Timesteps_num = as.numeric(as.character(Timesteps)),
    Density_num = as.numeric(as.character(Density)),
    Nodes_num = as.numeric(as.character(Nodes)),
    Regimes_num = as.numeric(as.character(Regimes))
  )

# Aggregate: Mean across RegimeIndex for each SimUID
dat_sim <- corr_results_for_agg %>%
  group_by(SimUID, Condition, SimID, 
           Timesteps_num, Density_num, Nodes_num, Regimes_num) %>%
  summarise(
    across(all_of(z_cols), mean, .names = "{.col}"),
    .groups = "drop"
  ) %>%
  rename(
    Timesteps = Timesteps_num,
    Density = Density_num,
    Nodes = Nodes_num,
    Regimes_fac = Regimes_num  # Rename temporarily
  )

cat("✓ Aggregation completed\n")
cat("  Rows before:", nrow(corr_results), "\n")
cat("  Rows after:", nrow(dat_sim), "\n")
cat("  → Each simulation now has equal weight!\n\n")

# Check: Gewichtungsproblem gelöst?
original_regime_dist <- table(corr_results_for_agg$Regimes_num)
cat("Original distribution (row-level):\n")
print(original_regime_dist)
cat("\nPercentages:\n")
print(round(100 * original_regime_dist / sum(original_regime_dist), 1))

sim_regime_dist <- table(dat_sim$Regimes_fac)
cat("\nAggregated distribution (simulation-level):\n")
print(sim_regime_dist)
cat("\nPercentages:\n")
print(round(100 * sim_regime_dist / sum(sim_regime_dist), 1))
cat("\n")

# -----------------------------------------------------------------------------
# PART 4: TRANSFORM & SCALE PREDICTORS
# -----------------------------------------------------------------------------

cat("=== PART 4: PREDICTOR TRANSFORMATION ===\n\n")

dat_sim <- dat_sim %>%
  mutate(
    # log(Timesteps) for diminishing returns
    logT = as.numeric(scale(log(Timesteps))),
    
    # Density and Nodes: centered and scaled
    Density_s = as.numeric(scale(Density)),
    Nodes_s = as.numeric(scale(Nodes)),
    
    # Regimes as FACTOR (not scaled!)
    Regimes = factor(Regimes_fac)
  ) %>%
  select(-Regimes_fac)  # Remove temporary variable

cat("✓ Predictors transformed\n")
cat("  • log(Timesteps) scaled: Mean =", round(mean(dat_sim$logT), 3), 
    ", SD =", round(sd(dat_sim$logT), 3), "\n")
cat("  • Density scaled: Mean =", round(mean(dat_sim$Density_s), 3), 
    ", SD =", round(sd(dat_sim$Density_s), 3), "\n")
cat("  • Nodes scaled: Mean =", round(mean(dat_sim$Nodes_s), 3), 
    ", SD =", round(sd(dat_sim$Nodes_s), 3), "\n")
cat("  • Regimes as factor with", nlevels(dat_sim$Regimes), "levels:", 
    paste(levels(dat_sim$Regimes), collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# PART 4b: RESIDUAL DIAGNOSTICS (before inference)
# -----------------------------------------------------------------------------

cat("=== PART 4b: RESIDUAL DIAGNOSTICS ===\n\n")

diagnostics_results <- list()

for (outcome in corr_cols) {
  
  outcome_z <- paste0(outcome, "_z")
  
  cat("Diagnostics for:", outcome, "\n")
  cat(strrep("-", 60), "\n")
  
  # Fit preliminary 2-way model for diagnostics
  formula_2way <- as.formula(
    paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2")
  )
  model_prelim <- lm(formula_2way, data = dat_sim)
  
  # --- 1. Breusch-Pagan Test (Heteroskedastizität) ---
  bp_test <- lmtest::bptest(model_prelim)
  
  cat("\n  Breusch-Pagan Test:\n")
  cat("    BP =", round(bp_test$statistic, 2), 
      ", df =", bp_test$parameter,
      ", p", ifelse(bp_test$p.value < 0.001, "< .001", 
                    paste("=", round(bp_test$p.value, 3))), "\n")
  
  if (bp_test$p.value < 0.05) {
    cat("    → Heteroskedastizität vorhanden → HC3-Korrektur gerechtfertigt\n")
  } else {
    cat("    → Keine signifikante Heteroskedastizität\n")
  }
  
  # --- 2. Shapiro-Wilk Test (Normalverteilung der Residuen) ---
  resids <- residuals(model_prelim)
  
  # Shapiro-Wilk hat ein Limit von n = 5000
  if (length(resids) <= 5000) {
    sw_test <- shapiro.test(resids)
  } else {
    # Bei großen Stichproben: Zufallsstichprobe von 5000
    set.seed(42)
    sw_test <- shapiro.test(sample(resids, 5000))
    cat("    (Shapiro-Wilk basiert auf Zufallsstichprobe von n = 5000)\n")
  }
  
  cat("\n  Shapiro-Wilk Test:\n")
  cat("    W =", round(sw_test$statistic, 4),
      ", p", ifelse(sw_test$p.value < 0.001, "< .001",
                    paste("=", round(sw_test$p.value, 3))), "\n")
  
  if (sw_test$p.value < 0.05) {
    cat("    → Abweichung von Normalverteilung\n")
    cat("    → Bei großem n durch ZGS abgefedert; Q-Q-Plot prüfen\n")
  } else {
    cat("    → Normalverteilung nicht verworfen\n")
  }
  
  # --- 3. Skewness & Kurtosis (deskriptiv) ---
  skew <- mean(((resids - mean(resids)) / sd(resids))^3)
  kurt <- mean(((resids - mean(resids)) / sd(resids))^4) - 3  # Excess Kurtosis
  
  cat("\n  Verteilungskennwerte der Residuen:\n")
  cat("    Skewness:", round(skew, 3), 
      ifelse(abs(skew) > 1, " ⚠ stark schief", ""), "\n")
  cat("    Excess Kurtosis:", round(kurt, 3),
      ifelse(abs(kurt) > 2, " ⚠ starke Abweichung", ""), "\n")
  
  # --- Ergebnisse speichern ---
  diagnostics_results[[outcome]] <- list(
    bp_test = bp_test,
    sw_test = sw_test,
    skewness = skew,
    kurtosis = kurt,
    heteroskedastic = bp_test$p.value < 0.05,
    non_normal = sw_test$p.value < 0.05
  )
  
  cat("\n")
}

# --- Zusammenfassung ---
cat("=== DIAGNOSTICS SUMMARY ===\n\n")
cat(sprintf("%-20s %-15s %-15s %-10s %-10s\n", 
            "Outcome", "BP p-value", "SW p-value", "Skewness", "Kurtosis"))
cat(strrep("-", 70), "\n")

for (outcome in corr_cols) {
  d <- diagnostics_results[[outcome]]
  cat(sprintf("%-20s %-15s %-15s %-10s %-10s\n",
              outcome,
              ifelse(d$bp_test$p.value < 0.001, "< .001", round(d$bp_test$p.value, 3)),
              ifelse(d$sw_test$p.value < 0.001, "< .001", round(d$sw_test$p.value, 3)),
              round(d$skewness, 3),
              round(d$kurtosis, 3)))
}

# HC3-Empfehlung
any_hetero <- any(sapply(diagnostics_results, `[[`, "heteroskedastic"))
cat("\n")
if (any_hetero) {
  cat("→ Mindestens ein Outcome zeigt signifikante Heteroskedastizität.\n")
  cat("  HC3-robuste Standardfehler werden für alle Modelle verwendet.\n")
} else {
  cat("→ Keine signifikante Heteroskedastizität gefunden.\n")
  cat("  HC3-robuste SE werden dennoch als Vorsichtsmaßnahme beibehalten.\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# PART 5: FIT MODELS (MAIN, 2-WAY, 3-WAY)
# -----------------------------------------------------------------------------

cat("=== PART 5: MODEL FITTING ===\n\n")

all_results <- list()

for (outcome in corr_cols) {
  
  outcome_z <- paste0(outcome, "_z")
  
  cat("Analyzing:", outcome, "\n")
  cat(strrep("-", 60), "\n")
  
  # 5.1: Main effects model (baseline)
  formula_main <- as.formula(
    paste0(outcome_z, " ~ logT + Density_s + Nodes_s + Regimes")
  )
  
  model_main <- lm(formula_main, data = dat_sim)
  
  # 5.2: 2-way interactions (PRIMARY MODEL)
  formula_2way <- as.formula(
    paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^2")
  )
  
  model_2way <- lm(formula_2way, data = dat_sim)
  
  # 5.3: 3-way interactions (EXPLORATORY)
  formula_3way <- as.formula(
    paste0(outcome_z, " ~ (logT + Density_s + Nodes_s + Regimes)^3")
  )
  
  model_3way <- lm(formula_3way, data = dat_sim)
  
  # 5.4: Model comparison
  comparison <- data.frame(
    Model = c("Main Effects", "2-Way (primary)", "3-Way (exploratory)"),
    df = c(
      length(coef(model_main)),
      length(coef(model_2way)),
      length(coef(model_3way))
    ),
    AIC = c(AIC(model_main), AIC(model_2way), AIC(model_3way)),
    BIC = c(BIC(model_main), BIC(model_2way), BIC(model_3way)),
    R2 = c(
      summary(model_main)$r.squared,
      summary(model_2way)$r.squared,
      summary(model_3way)$r.squared
    ),
    Adj_R2 = c(
      summary(model_main)$adj.r.squared,
      summary(model_2way)$adj.r.squared,
      summary(model_3way)$adj.r.squared
    )
  )
  
  comparison$Delta_AIC <- comparison$AIC - min(comparison$AIC)
  comparison$Delta_BIC <- comparison$BIC - min(comparison$BIC)
  
  # ANOVA F-Tests (for nested linear models)
  # Note: We use F-test, not LRT, because these are standard lm() models
  # LRT would be used for glm() or lmer() models
  anova_2v1 <- waldtest(model_main, model_2way, 
                        vcov = vcovHC(model_2way, type = "HC3"))
  anova_3v2 <- waldtest(model_2way, model_3way, 
                        vcov = vcovHC(model_3way, type = "HC3"))

  cat("\nModel Comparison:\n")
  print(comparison, digits = 3)
  
  cat("\nANOVA F-Test (2-Way vs Main Effects):\n")
  cat("  F(", anova_2v1$Df[2], ", ", anova_2v1$Res.Df[2], ") = ", 
      round(anova_2v1$F[2], 2),
      ", p ", ifelse(anova_2v1$`Pr(>F)`[2] < 0.001, "< .001", 
                     paste("=", round(anova_2v1$`Pr(>F)`[2], 3))), "\n", sep = "")
  
  # Decision about 3-way (direct AIC comparison: 2-way vs 3-way)
  delta_aic_3way <- comparison$AIC[2] - comparison$AIC[3]  # positive => 3-way better
  
  if (delta_aic_3way > 10) {
    cat("\n⚠ 3-Way model is substantially better (ΔAIC =", round(delta_aic_3way, 1), ")\n")
    cat("   → Mention in supplement, but focus on 2-Way for interpretation\n")
    include_3way <- TRUE
  } else {
    cat("\n✓ 2-Way model is sufficient (ΔAIC =", round(delta_aic_3way, 1), ")\n")
    include_3way <- FALSE
  }
  
  cat("\n")
  
  # Store results
  all_results[[outcome]] <- list(
    main = model_main,
    two_way = model_2way,
    three_way = model_3way,
    comparison = comparison,
    anova = list(
      two_vs_main = anova_2v1,
      three_vs_two = anova_3v2
    ),
    include_3way = include_3way,
    primary_model = model_2way  # ALWAYS use 2-way as primary
  )
}

# -----------------------------------------------------------------------------
# PART 6: EFFECT INTERPRETATION (SIMPLIFIED - Z-SCALE ONLY)
# -----------------------------------------------------------------------------

cat("\n=== PART 6: EFFECT SUMMARY ON Z-SCALE ===\n\n")

# Function for selective back-transformation (only for context)
z_to_r <- function(z) tanh(z)

interpretation_summary <- list()

for (outcome in corr_cols) {
  
  outcome_z <- paste0(outcome, "_z")
  model <- all_results[[outcome]]$primary_model
  coefs <- coef(model)
  
  cat("Outcome:", outcome, "\n")
  cat(strrep("=", 60), "\n\n")
  
  # --- 1. OVERALL CONTEXT (for text only) ---
  intercept_z <- coefs["(Intercept)"]
  
  # Calculate marginal mean across regime levels
  regime_coefs <- coefs[grepl("^Regimes", names(coefs))]
  regime_z_values <- c(
    intercept_z,
    intercept_z + regime_coefs["Regimes2"],
    intercept_z + regime_coefs["Regimes3"],
    intercept_z + regime_coefs["Regimes4"]
  )
  
  marginal_mean_z <- mean(regime_z_values, na.rm = TRUE)
  marginal_mean_r <- z_to_r(marginal_mean_z)
  
  # Min and max for context
  min_r <- z_to_r(min(regime_z_values, na.rm = TRUE))
  max_r <- z_to_r(max(regime_z_values, na.rm = TRUE))
  
  cat("CONTEXT (for text interpretation):\n")
  cat("  Overall mean (marginal): z =", round(marginal_mean_z, 3), 
      "→ r =", round(marginal_mean_r, 3), "\n")
  cat("  Range across regimes: r =", round(min_r, 3), "to", round(max_r, 3), "\n\n")
  
  # --- 2. MAIN EFFECTS (Z-SCALE ONLY) ---
  cat("MAIN EFFECTS (z-scale, for tables):\n")
  
  main_predictors <- c("logT", "Density_s", "Nodes_s")
  for (pred in main_predictors) {
    if (pred %in% names(coefs)) {
      beta <- coefs[pred]
      cat("  ", pred, ": β =", round(beta, 3), "\n")
    }
  }
  
  cat("\n  Regime effects (vs. Regimes1):\n")
  for (reg_name in names(regime_coefs)) {
    cat("  ", reg_name, ": β =", round(regime_coefs[reg_name], 3), "\n")
  }
  
  # --- 3. TOP INTERACTIONS (Z-SCALE ONLY) ---
  interaction_names <- grep(":", names(coefs), value = TRUE)
  
  if (length(interaction_names) > 0) {
    interaction_coefs <- coefs[interaction_names]
    interaction_sorted <- sort(abs(interaction_coefs), decreasing = TRUE)
    
    cat("\n  Top 3 interactions:\n")
    for (i in 1:min(3, length(interaction_sorted))) {
      name <- names(interaction_sorted)[i]
      value <- interaction_coefs[name]
      cat("  ", gsub(":", " × ", name), ": β =", round(value, 3), "\n")
    }
  }
  
  cat("\n\n")
  
  # Store minimal summary
  interpretation_summary[[outcome]] <- list(
    marginal_mean_z = marginal_mean_z,
    marginal_mean_r = marginal_mean_r,
    range_r = c(min_r, max_r),
    main_effects_z = coefs[main_predictors[main_predictors %in% names(coefs)]],
    regime_effects_z = regime_coefs,
    top_interactions = head(sort(abs(interaction_coefs), decreasing = TRUE), 3)
  )
}

# --- Summary table for quick reference ---
cat("=== QUICK REFERENCE: MARGINAL MEANS ===\n\n")

marginal_means_df <- data.frame(
  Outcome = corr_cols,
  Mean_z = sapply(interpretation_summary, function(x) x$marginal_mean_z),
  Mean_r = sapply(interpretation_summary, function(x) x$marginal_mean_r),
  Min_r = sapply(interpretation_summary, function(x) x$range_r[1]),
  Max_r = sapply(interpretation_summary, function(x) x$range_r[2])
)

print(marginal_means_df, digits = 3)

# Export this simple summary
write.csv(marginal_means_df, 
          file = "interpretation_marginal_means.csv", 
          row.names = FALSE)

cat("\n✓ Exported: interpretation_marginal_means.csv\n")
cat("  Use these r-values for contextualizing results in text.\n")
cat("  All other effects should be reported on z-scale.\n\n")

# -----------------------------------------------------------------------------
# PART 7: EXPORT RESULTS
# -----------------------------------------------------------------------------

cat("=== PART 7: EXPORT RESULTS ===\n\n")

# 7.1: Coefficients for primary models (2-way)
for (outcome in corr_cols) {
  
  model <- all_results[[outcome]]$primary_model
  
  # Coefficient table with robust (HC3) standard errors
  robust_ct <- lmtest::coeftest(model, vcov. = sandwich::vcovHC(model, type = "HC3"))
  coef_df <- data.frame(
    Predictor = rownames(robust_ct),
    Estimate = robust_ct[, 1],
    `Std. Error` = robust_ct[, 2],
    `df` = model$df.residual,
    `t value` = robust_ct[, 3],
    `p_raw` = robust_ct[, 4],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  coef_df$df <- model$df.residual
  
  # FDR correction within each outcome model
  coef_df$`p (BH-adj.)` <- p.adjust(coef_df$p_raw, method = "BH")
  
  # Significance stars based on FDR-adjusted p-values
  coef_df$`sig.` <- ifelse(coef_df$`p (BH-adj.)` < 0.001, "***",
                        ifelse(coef_df$`p (BH-adj.)` < 0.01, "**",
                               ifelse(coef_df$`p (BH-adj.)` < 0.05, "*", "")))
  
  # Final table: only adjusted p-values
  coef_df <- coef_df[, c("Predictor", "Estimate", "Std. Error", "df", "t value",
                         "p (BH-adj.)", "sig.")]
  
  # Clean predictor names
  coef_df$Predictor <- gsub(":", " × ", coef_df$Predictor)
  
  # Export to HTML
  kable(coef_df, caption = paste("Results_2way_", outcome), digits =3, format = "html", row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    save_kable(file = paste0("Results_2way_", outcome, ".html"))
  
  # Export to LaTeX (digits=3)
  print(
    xtable(coef_df, caption = paste("Results_2way_", outcome), digits = 3, row.names = FALSE),
    file = paste0("Results_2way_", outcome, ".tex"),
    include.rownames = FALSE
  )
}

# 7.2: Model comparison table
comparison_all <- do.call(rbind, lapply(names(all_results), function(outcome) {
  comp <- all_results[[outcome]]$comparison
  comp$Outcome <- outcome
  comp
}))

write.csv(comparison_all,
          file = "model_comparison_all_outcomes.csv",
          row.names = FALSE)

cat("✓ Exported: model_comparison_all_outcomes.csv\n")

# 7.3: Diagnostics (PDF)
pdf("diagnostics_all_models.pdf", width = 10, height = 8)

for (outcome in corr_cols) {
  model <- all_results[[outcome]]$primary_model
  
  par(mfrow = c(2, 2))
  plot(model, main = paste(outcome, "- 2-Way Model"))
}

dev.off()


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


