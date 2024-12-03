

t_mean <- list()
for (i in seq_along(T)) {
  timesteps <- T[i]
  t_subset <- MSAR_results %>% filter(Timesteps == timesteps)
  t_mean[[paste0(timesteps, "_Timesteps")]] <- mean(t_subset$N)
}

t_kruskal <- kruskal.test(N ~ Timesteps, MSAR_results)
d_kruskal <- kruskal.test(N ~ Density, MSAR_results)
n_kruskal <- kruskal.test(N ~ Nodes, MSAR_results)


library(dunn.test)
MSAR_results$Timesteps <- as.factor(MSAR_results$Timesteps)
t_dunn <- dunn.test(MSAR_results$N, MSAR_results$Timesteps, method = "bonferroni")

MSAR_results$Nodes <- as.factor(MSAR_results$Nodes)
d_dunn <- dunn.test(MSAR_results$N, MSAR_results$Nodes, method = "bonferroni")
n_dunn <- dunn.test(N ~ Nodes, MSAR_results, method = "bonferroni")



corr_list <- list()
corr_results <- list()

# Extract stats data from MSAR results list
for (t in seq_along(T)) {
  for (density in seq_along(Density)) {
    for (nodes in seq_along(N)) {
      for (regimes in seq_along(M)) {
        for (ts in 1:length(MSAR_models[[t]][[density]][[nodes]][[regimes]][["MSAR_models"]])) {
          for (r in 1:regimes) {
          result <- MSAR_models[[t]][[density]][[nodes]][[regimes]][["MSAR_models"]][[ts]][[r]]
          temp <- data.frame(
          Timesteps = T[t],
          Density = Density[density],
          Nodes = N[nodes],
          Regimes = M[regimes],
          Wtemp_corr = result[["corr. Wtemp"]],
          Wtemp_ac_corr = result[["corr. Wtemp ac"]],
          Wcont_corr = result[["corr. Wcont"]],
          result[["corr. Wcont ac"]]
        )
        corr_list <- append(corr_list, list(temp))
          }
        }
      }
    }
  }
}

# Combine all dataframes in list in one dataframe
corr_results <- do.call(rbind, corr_list)

# Omit na values
corr_results <- na.omit(corr_results)



####################################
corr_list <- list()

# Extract stats data from MSAR results list
for (t in seq_along(T)) {
  for (density in seq_along(Density)) {
    for (nodes in seq_along(N)) {
      for (regimes in seq_along(M)) {
        # Fetch models
        models <- MSAR_models[[t]][[density]][[nodes]][[regimes]][["MSAR_models"]]
        
        for (ts in seq_along(models)) {
          for (r in seq_along(models[[ts]])) {
            result <- models[[ts]][[r]]
            
            # Create dataframe for current parameters
            temp <- data.frame(
              Timesteps = T[t],
              Density = Density[density],
              Nodes = N[nodes],
              Regimes = M[regimes],
              Wtemp_corr = result[["corr. Wtemp"]],
              Wtemp_ac_corr = result[["corr. Wtemp ac"]],
              Wcont_corr = result[["corr. Wcont"]],
              Wcont_ac_corr = result[["corr. Wcont ac"]]
            )
            
            # Add dataframe to list
            corr_list <- append(corr_list, list(temp))
          }
        }
      }
    }
  }
}

# Combine all dataframes
corr_results <- do.call(rbind, corr_list)

# Remove NA
corr_results <- na.omit(corr_results)
###############################################




ANOVA <- aov(Wtemp_corr ~ Timesteps * Density * Nodes * Regimes, corr_results)
residuals <- residuals(ANOVA)
shapiro.test(residuals)

qqnorm(residuals)
qqline(residuals, col = "red")

hist(residuals, breaks = 50, main = "Histogram of Residuals", xlab = "Residuals")


plot(density(residuals), main = "Density Plot of Residuals", xlab = "Residuals")
curve(dnorm(x, mean=mean(residuals), sd=sd(residuals)), add=TRUE, col="red")


corr_results$Timesteps <- as.factor(corr_results$Timesteps)
corr_results$Density <- as.factor(corr_results$Density)
corr_results$Nodes <- as.factor(corr_results$Nodes)
corr_results$Regimes <- as.factor(corr_results$Regimes)



# Robuste ANOVA mit trimmed means
library(ARTool)
Wtemp_art <- art(Wtemp_corr ~ Timesteps * Density * Nodes * Regimes, data = corr_results)
anova_results <- anova(Wtemp_art)
print(anova_results)
##########################################

library(lmPerm)

# Specify dependent and independent variables
dependent_vars <- colnames(corr_results)[5:8]
independent_vars <- colnames(corr_results)[1:4]

permanova_results <- list()

# Loop over dependent variables
for (dep_var in dependent_vars) {
  # Create formula for ANOVA
  formula <- as.formula(paste(dep_var, "~", paste(independent_vars, collapse = " + "), "+",
                              paste(combn(independent_vars, 2, FUN = paste, collapse = ":"), collapse = " + "), "+",
                              paste(combn(independent_vars, 3, FUN = paste, collapse = ":"), collapse = " + "), "+",
                              paste(independent_vars, collapse = ":"), collapse = " "))
  
  # Calculate permutations ANOVA
  model <- aovp(formula, data = corr_results, perm = "Prob")
  
  # Extract summary
  summary_model <- summary(model)
  
  # Extract effects, degrees of freedom, and mean squares.n
  effects <- rownames(summary_model[[1]])
  df <- summary_model[[1]][, "Df"]
  mean_sq <- summary_model[[1]][, "R Mean Sq"]
  
  # Extract esidual Mean Square and p-values
  effects <- trimws(effects) 
  residual_mean_sq <- mean_sq[effects == "Residuals"]
  if (length(residual_mean_sq) == 0 || is.na(residual_mean_sq)) {
    stop("Residual Mean Square konnte nicht berechnet werden.")
  }
  
  p_values <- summary_model[[1]][, "Pr(Prob)"]
  
  # Calculate F-values without residuals
  f_values <- mean_sq / residual_mean_sq
  f_values[effects == "Residuals"] <- NA  
  
  # Store results in dataframe
  results_df <- data.frame(
    Effect = effects,
    Df = df,
    Mean_Sq = mean_sq,
    F_value = f_values,
    P_value = p_values
  )
  
  # Store results in list
  permanova_results[[dep_var]] <- results_df
}

# Show results
for (dep_var in names(permanova_results)) {
  cat("\nResults for dependent variable:", dep_var, "\n")
  print(permanova_results[[dep_var]])
}





##########################################
library(ggplot2)
library(dplyr)

# Daten zusammenfassen
summary_data <- corr_results %>%
  group_by(Timesteps, Density, Nodes, Regimes) %>%
  summarise(mean_Wtemp_corr = mean(Wtemp_corr, na.rm = TRUE))

# Plot erstellen

# Plot erstellen
Wtemp_plot <- ggplot(summary_data, aes(x = Timesteps, y = mean_Wtemp_corr, color = Nodes, group = Nodes)) +
  geom_line() +
  geom_point() +
  facet_grid(Density ~ Regimes, labeller = label_both) +
  labs(title = "Wtemp mean correlations",
       x = "Timesteps",
       y = "Mean correlations",
       color = "Nodes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # X-Achsenbeschriftungen schräg darstellen
  )

Wtemp_plotly <- ggplotly(Wtemp_plot)
htmlwidgets::saveWidget(Wtemp_plotly, "Plots/Wtemp_plot.html")

ggsave(filename = "Plots/Wtemp_plot.svg",
       plot = Wtemp_plot,
       width = 30,        # Breite des Plots
       height = 22.5,       # Höhe des Plots
       units = "cm"       # Einheit der Abmessungen
)






#######################################


library(ggplot2)
library(dplyr)
library(cowplot)

# Daten zusammenfassen
summary_data <- corr_results %>%
  group_by(Timesteps, Density, Nodes, Regimes) %>%
  summarise(mean_Wtemp_corr = mean(Wtemp_corr, na.rm = TRUE),
            sd_Wtemp_corr = sd(Wtemp_corr, na.rm = TRUE))

# HoverInfo-Spalte erstellen
summary_data <- summary_data %>% mutate(HoverInfo = paste("Timesteps:", Timesteps,
                                                          "<br>Nodes:", Nodes,
                                                          "<br>Mean:", round(mean_Wtemp_corr, 2),
                                                          "<br>Sd:", round(sd_Wtemp_corr, 2)
))

# Position für versetzte Punkte
dodge <- position_dodge(width = 0.5)

# Basisplot erstellen
Wtemp_plot <- ggplot() +
  # Einzelwerte als Scatterplot
  geom_jitter(data = corr_results,
              aes(x = Timesteps, y = Wtemp_corr, color = Nodes),
              position = dodge, alpha = 0.2, size = 0.1) +
  # Aggregierte Daten als Linie und Punkte
  geom_line(data = summary_data,
            aes(x = as.numeric(Timesteps), y = mean_Wtemp_corr, color = Nodes, group = Nodes),
            position = dodge) +
  geom_point(data = summary_data,
             aes(x = as.numeric(Timesteps), y = mean_Wtemp_corr, color = Nodes, text = HoverInfo),
             position = dodge, size = 1) +
  # Facetierung
  facet_grid(Density ~ Regimes, labeller = label_value) +
  # Labels und Design
  labs(title = "Wtemp mean correlations",
       x = "Timesteps",
       y = "Mean correlations",
       color = "Nodes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               size = 8),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    strip.text.y = element_text(vjust = -0.25),
    panel.spacing = unit(0.2, "in")
  )

# Zusätzliche Labels erstellen
label_x_right <- ggplot() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Density", angle = -90, size = 4, hjust = 0)

label_y_top <- ggplot() +
  theme_void() +
  xlim(0, 1) +       # Definiere den horizontalen Bereich
  ylim(0, 1) +       # Definiere den vertikalen Bereich
  annotate("text", x = 0.4725, y = 0.5, label = "Regimes", size = 4, hjust = 0.5)


# Kombination der Plots mit `cowplot`
final_plot <- ggdraw() +
  draw_plot(Wtemp_plot, 0, 0, 1, 1) +                      # Hauptplot
  draw_plot(label_x_right, 0.96, 0.08, 0.03, 0.8) +        # X-Beschriftung rechts
  draw_plot(label_y_top, 0.1, 0.875, 0.84, 0.05)            # Y-Beschriftung oben

# Speichern als PDF
ggsave(filename = "Plots/Wtemp_plot_with_labels.pdf",
       plot = final_plot,
       width = 10,       # Breite des Plots
       height = 8,       # Höhe des Plots
       units = "in",
       dpi = 600         # Auflösung in DPI
)








