library(plotly)


## Extract stats from results list to dataframe
# Initialize list to store stats data for each combination from Density*T*N*M
data_list <- list()
MSAR_results <- list()

# Extract stats data from MSAR results list
for (t in names(MSAR_dynamics_list)) {
  for (density in names(MSAR_dynamics_list[[t]])) {
    for (nodes in names(MSAR_dynamics_list[[t]][[density]])) {
      for (regimes in names(MSAR_dynamics_list[[t]][[density]][[nodes]])) {
        stats <- MSAR_dynamics_list[[t]][[density]][[nodes]][[regimes]]$Stats
        temp <- data.frame(
          Timesteps = as.numeric(gsub("_Timesteps", "", t)),
          Density = as.numeric(gsub("%", "", gsub("Density_", "", density))),
          Nodes = as.numeric(gsub("_Nodes", "", nodes)),
          Regimes = MSAR_results$Regimes <- as.numeric(gsub("_Regimes", "", regimes)),
          N = stats$Wtemp_corr$N,
          Wtemp_corr_mean = stats$Wtemp_corr$Mean,
          Wtemp_corr_sd = stats$Wtemp_corr$Sd,
          Wtemp_ac_corr_mean = stats$Wtemp_ac_corr$Mean,
          Wtemp_ac_corr_sd = stats$Wtemp_ac_corr$Sd,
          Wcont_corr_mean = stats$Wcont_corr$Mean,
          Wcont_corr_sd = stats$Wcont_corr$Sd,
          Wcont_ac_corr_mean = stats$Wcont_ac_corr$Mean,
          Wcont_ac_corr_sd = stats$Wcont_ac_corr$Sd
        )
        data_list <- append(data_list, list(temp))
      }
    }
  }
}

# Combine all dataframes in list in one dataframe
MSAR_results <- do.call(rbind, data_list)

# Omit na values
MSAR_results <- na.omit(MSAR_results)
MSAR_results$N <- MSAR_results$N / MSAR_results$Regimes


## Get z-values for plot
# Initialize list for z-values (means)
z_matrices <- list()

global_min <- Inf
global_max <- -Inf

# Create matrix for z-values for each regime
for (regime in M) {
  df_subset <- subset(MSAR_results, Regimes == regime)

  # Initialize z-matrix
  z_matrix <- matrix(NA, nrow = length(N), ncol = length(T))
  
  # Fill z-matrix with mean values from Wtemp_corr_mean
  for (i in seq_along(N)) {
    for (j in seq_along(T)) {
      value <- df_subset$Wtemp_corr_mean[df_subset$Nodes == N[i] & df_subset$Timesteps == T[j]]
      if (length(value) > 0) {
        z_matrix[i, j] <- value
        global_min <- min(global_min, value, na.rm = TRUE)
        global_max <- max(global_max, value, na.rm = TRUE)
      }
    }
  }
  z_matrices[[paste0(regime, "_Regimes")]] <- z_matrix
}


## Create 3D surface plots
# Initialize plot
fig <- plot_ly()

# Set plot colours
colorscale <- list(c(0, "blue"),
                   c(0.33, "cyan"),
                   c(0.66, "yellow"),
                   c(1, "red")
                   )

# Füge Oberflächen für jeden Wert in Regimes hinzu
annotations <- list()
shapes <- list()

for (regime in seq_along(M)) {

  # Add surface for each regime to plot
  fig <- fig %>% add_surface(
    x = T,
    y = N,
    z = z_matrices[[regime]],
    name = paste("Regimes", M[regime]),
    opacity = 1,
    colorscale = colorscale,
    zmin = global_min,
    zmax = global_max,
    showscale = regime == 1
    )
  
  
  # Finde die Koordinaten für die Anmerkung
  x_annotation <- mean(T)
  y_annotation <- mean(N)
  z_annotation <- global_max
  
  # Füge die Anmerkung hinzu
  annotations <- append(annotations, list(
    list(
      x = x_annotation,
      y = y_annotation,
      z = z_annotation,
      text = paste("Regimes", regime),
      xanchor = "center",
      yanchor = "bottom",
      showarrow = TRUE,
      arrowhead = 2,
      ax = 0,
      ay = -40
    )
  ))
  
  # Füge die Linie hinzu
  shapes <- append(shapes, list(
    list(
      type = "line",
      x0 = x_annotation,
      y0 = y_annotation,
      z0 = z_annotation,
      x1 = x_annotation,
      y1 = y_annotation,
      z1 = z_annotation - 0.1 * (global_max - global_min),
      line = list(color = "black", width = 2)
    )
  ))
  }

# Set layout for plot
fig <- fig %>% layout(
  title = "3D Surface Plot für verschiedene Regimes",
  scene = list(
    xaxis = list(title = "Timesteps", tickvals = T),
    yaxis = list(title = "Nodes", tickvals = N),
    zaxis = list(title = "Mean")
  )
)


# Save combined plot as html file
htmlwidgets::saveWidget(fig, "Plots/combined_3d_surface_plot.html")






# Load libraries
library(ggplot2)
library(plotly)

# Erstelle benutzerdefinierte Hover-Informationen
MSAR_results <- MSAR_results %>% mutate(
  hover_text = paste(
    "Nodes:", Nodes, "<br>",
    "Timesteps:", Timesteps, "<br>",
    "Mean Corr:", round(Wtemp_corr_mean, 2), "<br>",
    "SD:", round(Wtemp_corr_sd, 2), "<br>",
    "Regimes:", Regimes
  )
)


# Erstelle den 3D Bubbleplot
fig <- plot_ly(
  MSAR_results, 
  x = ~Nodes, 
  y = ~Wtemp_corr_mean, 
  z = ~Timesteps,
  #size = ~Wtemp_corr_sd * 0.5,
  color = ~factor(Regimes),
  colors = c('#0d0887', '#9c179e', '#ed7953', '#f0f921'),
  marker = list(
   size = ~Wtemp_corr_sd * 50,
    sizemode = 'diameter',
    line = list(width = 0.1, color = 'black')
    ),
  text = ~hover_text,
  hoverinfo = 'text',
  type = 'scatter3d', 
  mode = 'markers'
) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Nodes'),
      yaxis = list(title = 'Wtemp mean correlations'),
      zaxis = list(title = 'Timesteps')
    ),
    title = '3D Bubbleplot für Timesteps, Nodes und Means',
    showlegend = TRUE,
    legend = list(
      title = list(text = 'Regimes'),
      itemsizing = 'constant'
    )
  )

# Zeige den Plot an
htmlwidgets::saveWidget(fig, "3DBubblePlot.html")





















# Definiere die Plasma-Farbpalette
plasma_colors <- c('#0d0887', '#9c179e', '#ed7953', '#f0f921')
jet_colors <- c("blue", "cyan", "yellow", "red")

# Erstelle benutzerdefinierte Hover-Informationen
MSAR_results <- MSAR_results %>% mutate(
  hover_text = paste(
    "Nodes:", Nodes, "<br>",
    "Timesteps:", Timesteps, "<br>",
    "Mean Corr:", round(Wtemp_corr_mean, 2), "<br>",
    "SD:", round(Wtemp_corr_sd, 2), "<br>",
    "Regimes:", Regimes
  )
)

# Erstelle den 3D-Plot mit Linien und Bubbles
fig <- plot_ly()

# Fügen Sie die Linien und Bubbles für jede Kombination von Regimes und Nodes hinzu
for (regime in unique(MSAR_results$Regimes)) {
  for (node in unique(MSAR_results$Nodes)) {
    subset_data <- MSAR_results %>% filter(Regimes == regime, Nodes == node)
    
    # Sortiere die Daten nach Timesteps
    subset_data <- subset_data %>% arrange(Timesteps)
    
    # Füge die Linien hinzu
    for (i in 1:(nrow(subset_data) - 1)) {
      fig <- fig %>%
        add_trace(
          x = subset_data$Nodes[i:(i + 1)],
          y = subset_data$Wtemp_corr_mean[i:(i + 1)],
          z = subset_data$Timesteps[i:(i + 1)],
          type = 'scatter3d',
          mode = 'lines',
          color = factor(subset_data$Regimes[i:(i + 1)]),
          line = list(color = 'black',
                      width = subset_data$Wtemp_corr_sd[i]*55
                      #opacity = 0.70
          ),
          hoverinfo = 'none',
          showlegend = FALSE
        )
    }
    
    # Füge die Linien hinzu
    for (i in 1:(nrow(subset_data) - 1)) {
      fig <- fig %>%
        add_trace(
          x = subset_data$Nodes[i:(i + 1)],
          y = subset_data$Wtemp_corr_mean[i:(i + 1)],
          z = subset_data$Timesteps[i:(i + 1)],
          type = 'scatter3d',
          mode = 'lines',
          color = factor(subset_data$Regimes[i:(i + 1)]),
          line = list(color = jet_colors[as.numeric(regime)],
                      width = subset_data$Wtemp_corr_sd[i]*50
                      #opacity = 0.70
                      ),
          hoverinfo = 'none',
          showlegend = FALSE
        )
    }
    
    # Füge die Bubbles hinzu
    fig <- fig %>%
      add_trace(
        x = subset_data$Nodes, 
        y = subset_data$Wtemp_corr_mean, 
        z = subset_data$Timesteps,
        type = 'scatter3d',
        mode = 'markers',
        color = factor(subset_data$Regimes),
                marker = list(
                  size = subset_data$Wtemp_corr_sd * 50,
                  color = jet_colors[as.numeric(regime)],
                  sizemode = 'diameter',
                  line = list(width = 0.1, color = 'black')
                  #opacity = 0.7
                  ),
        text = subset_data$hover_text,
        hoverinfo = 'text',
        showlegend = FALSE
      )
  }
}



# Eigene Legende hinzufügen
for (regime in unique(MSAR_results$Regimes)) {
  fig <- fig %>%
    add_trace(
      x = c(3), y = c(0), z = c(0),
      type = 'scatter3d',
      mode = 'markers',
      marker = list(size = 1, color = jet_colors[as.numeric(regime)]),
      showlegend = TRUE,
      name = paste("Regime", regime),
      hoverinfo = 'none'
    )
}

fig <- fig %>%
  layout(
    scene = list(
      xaxis = list(title = 'Nodes'),
      yaxis = list(title = 'Wtemp mean correlations'),
      zaxis = list(title = 'Timesteps')
    ),
    title = '3D Plot für Timesteps, Nodes und Means mit Bubbles und Linien',
    showlegend = TRUE,
    legend = list(
      title = list(text = 'Regimes'),
      itemsizing = 'constant'
    )
  )

# Zeige den Plot an
htmlwidgets::saveWidget(fig, "Plots/3DPlot_with_Custom_Hover_Info.html")








library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(colorspace)

scale_factor <- 0.2

# Funktion zur Verschiebung der Timesteps basierend auf der Anzahl der Nodes
shift_timesteps <- function(timesteps, nodes) {
  shift_amount <- 50 * (nodes - mean(nodes))  # Verschiebungsgröße anpassen
  return(timesteps + shift_amount)
}

# Erstelle die Plots
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime)
  
  # Füge eine kleine Verschiebung zu den Timesteps basierend auf Nodes hinzu
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = shift_timesteps(Timesteps, Nodes))
  
  p <- ggplot(regime_data, aes(x = shifted_timesteps, y = Wtemp_corr_mean, color = factor(Nodes))) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = Wtemp_corr_mean - (Wtemp_corr_sd * scale_factor), 
                      ymax = Wtemp_corr_mean + (Wtemp_corr_sd * scale_factor)), 
                  width = 0.2, position = position_dodge(width = 0.2)) +
    labs(title = paste("Regime", regime, ": Mean Correlation over Timesteps"),
         x = "Timesteps", y = "Wtemp Corr Mean", color = "Nodes") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_viridis_d()
  
  print(p)
}

# Optional: Speichere die Plots als PDF
pdf("Plots/Regimes_Plots.pdf")
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime)
  
  # Füge eine kleine Verschiebung zu den Timesteps basierend auf Nodes hinzu
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = shift_timesteps(Timesteps, Nodes))
  
  p <- ggplot(regime_data, aes(x = shifted_timesteps, y = Wtemp_corr_mean, color = factor(Nodes))) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = Wtemp_corr_mean - (Wtemp_corr_sd * scale_factor), 
                      ymax = Wtemp_corr_mean + (Wtemp_corr_sd * scale_factor)), 
                  width = 0.2, position = position_dodge(width = 0.2)) +
    labs(title = paste("Regime", regime, ": Mean Correlation over Timesteps"),
         x = "Timesteps", y = "Wtemp Corr Mean", color = "Nodes") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_viridis_d()
  
  print(p)
}
dev.off()





# Skalierungsfaktor für die Linienbreite
scale_factor <- 10

# Erstelle die Plots
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime)
  

  
  p <- ggplot(regime_data, aes(x = Timesteps,
                               y = Wtemp_corr_mean,
                               color = factor(Nodes),
                               fill = factor(Nodes))) +
    geom_point() +
    geom_ribbon(aes(ymin = Wtemp_corr_mean - Wtemp_corr_sd,
                    ymax = Wtemp_corr_mean + Wtemp_corr_sd),
                alpha = 0.2,
                color = NA,
                show.legend = FALSE) +
    geom_errorbar(aes(ymin = Wtemp_corr_mean - Wtemp_corr_sd,
                      ymax = Wtemp_corr_mean + Wtemp_corr_sd,
                      color = factor(Nodes)), # Ändere die Farbe der Errorbars
                  width = 1) +
    labs(title = paste("Regimes", regime, ": Mean Correlation over Timesteps"),
         x = "Timesteps", y = "Wtemp mean correlation", color = "Nodes") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_discrete_qualitative() + 
    scale_fill_discrete_qualitative(guide = "none")
    
  
  print(p)
}


scaling <- 1

# Optional: Speichere die Plots als PDF
pdf("Plots/Regimes_Plots.pdf")
scaling <- 0.07
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime)

  
  p <- ggplot(regime_data, aes(x = Timesteps,
                               y = Wtemp_corr_mean,
                               color = factor(Nodes),
                               fill = factor(Nodes))) +
    #geom_line(aes(size = 0, color = factor(Nodes)), alpha = 1, show.legend = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin = Wtemp_corr_mean - Wtemp_corr_sd * scaling,
                    ymax = Wtemp_corr_mean + Wtemp_corr_sd * scaling),
                alpha = 0.2,
                color = NA,
                show.legend = FALSE) +
    labs(title = paste("Regimes", regime, ": Mean Correlation over Timesteps"),
         x = "Timesteps", y = "Wtemp mean correlations", color = "Nodes") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_discrete_qualitative() + 
    scale_fill_discrete_qualitative(guide = "none")
   
  
  print(p)
}
dev.off()



library(colorspace)


# Skalierungsfaktor für die Breite der Ribbons
scale_factor <- 0.5 # Passe diesen Wert nach Bedarf an

# Erstelle die Plots
pdf("Plots/Regimes_Plots.pdf")
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime) 

  
  # Berechne die Ribbon-Werte unabhängig vom Mittelwert
  regime_data <- regime_data %>%
    mutate(ymin = Wtemp_corr_mean - Wtemp_corr_sd * scale_factor,
           ymax = Wtemp_corr_mean + Wtemp_corr_sd * scale_factor)
  
  p <- ggplot(regime_data, aes(x = Timesteps,
                               y = Wtemp_corr_mean,
                               color = factor(Nodes),
                               fill = factor(Nodes))) +
    geom_point(size = 4) +  # Vergrößere die Marker
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
                alpha = 0.2,
                color = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 0.7) +  # Linie für den Mittelwert
    geom_text(aes(label = round(Wtemp_corr_sd, 2)), 
              color = "black",  # Setze die Schriftfarbe auf schwarz
              vjust = 0.5, hjust = 0.5, size = 2.5,  # Zentriere den Text in den Markern
              fontface = "bold") +  # Setze den Text in Fettdruck
    labs(title = paste("Regime", regime, ": Mean Correlation over Timesteps"),
         x = "Timesteps", y = "Wtemp mean correlation", color = "Nodes") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_discrete_qualitative() +
    scale_fill_discrete_qualitative(guide = "none")
  
  print(p)
}
dev.off()





# # Funktion zur Verschiebung der Timesteps basierend auf der Anzahl der Nodes
# shift_timesteps <- function(timesteps, nodes) {
#   shift_amount <- 50 * nodes - mean(nodes)  # Verschiebungsgröße anpassen
#   return(timesteps + shift_amount)
# }

# Skalierungsfaktor für die Breite der Ribbons
scale_factor <- 0.2  # Passe diesen Wert nach Bedarf an
x_offset <- 50  # Betrag für die Verschiebung der x-Werte

# Erstelle die Plots
pdf("Plots/Regimes_Plots.pdf")
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime) 
  
  # Füge eine kleine Verschiebung zu den Timesteps basierend auf Nodes hinzu
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = Timesteps + (as.numeric(factor(Nodes)) - 1) * x_offset)  # Verschiebung ohne Abhängigkeit vom Mittelwert
  
  p <- ggplot(regime_data, aes(x = shifted_timesteps,
                               y = Wcont_corr_mean,
                               color = factor(Nodes),
                               fill = factor(Nodes))) +
    geom_point(size = 1.5,
               show.legend = FALSE) +  
    geom_ribbon(aes(ymin = Wcont_corr_mean - Wcont_corr_sd * scale_factor, 
                    ymax = Wcont_corr_mean + Wcont_corr_sd * scale_factor),
                alpha = 0.2,
                color = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 0.7) +  # Linie für den Mittelwert
    # geom_text(aes(label = round(Wtemp_corr_sd, 2)),
    #           color = "black",  # Setze die Schriftfarbe auf schwarz
    #           vjust = 0.5, hjust = 0.5, size = 2.5 # Zentriere den Text in den Markern
    #           ) +  # Setze den Text in Fettdruck
    geom_segment(aes(x = shifted_timesteps, xend = shifted_timesteps, 
                     y = Wcont_corr_mean - Wcont_corr_sd * scale_factor,
                     yend = Wcont_corr_mean + Wcont_corr_sd * scale_factor, 
                     color = factor(Nodes)), 
                 linetype = "solid",
                 linewidth = 0.4,
                 alpha = 0.85) +  # Linie in der Farbe von Nodes
    labs(title = paste("Regime", regime, ": Mean Correlation over Timesteps"),
         x = "Timesteps", y = "Wcont mean correlation", color = "Nodes") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_discrete_qualitative() +
    scale_fill_discrete_qualitative(guide = "none")
  
  print(p)
}
dev.off()







library(plotly)
library(dplyr)
library(RColorBrewer)

# Skalierungsfaktor für die Breite der Ribbons
scale_factor <- 0.2  # Passe diesen Wert nach Bedarf an
x_offset <- 50  # Betrag für die Verschiebung der x-Werte

# Erstelle den Plot für jedes Regime
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime)
  
  # Füge eine kleine Verschiebung zu den Timesteps basierend auf Nodes hinzu
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = Timesteps + (as.numeric(factor(Nodes)) - 1) * x_offset)
  
  # Definiere eine Farbpalette basierend auf der Anzahl der Nodes
  node_colors <- RColorBrewer::brewer.pal(length(unique(regime_data$Nodes)), "Set1")
  
  # Erstelle das Plotly-Objekt
  fig <- plot_ly()
  
  # Füge die Daten für jede Node-Gruppe hinzu
  for (i in seq_along(unique(regime_data$Nodes))) {
    node <- unique(regime_data$Nodes)[i]
    node_data <- regime_data %>% filter(Nodes == node)
    
    # Bestimme die Farbe für den aktuellen Node
    node_color <- node_colors[i]
    
    # Ribbon darstellen
    fig <- fig %>%
      add_trace(x = c(node_data$shifted_timesteps, rev(node_data$shifted_timesteps)), 
                y = c(node_data$Wcont_corr_mean - node_data$Wcont_corr_sd * scale_factor, 
                      rev(node_data$Wcont_corr_mean + node_data$Wcont_corr_sd * scale_factor)),
                type = 'scattergl', mode = 'lines',
                fill = 'toself', fillcolor = sprintf('rgba(%s,0.2)', paste(col2rgb(node_color), collapse = ',')),
                line = list(color = 'transparent', width = 0),
                name = paste("Node", node), legendgroup = paste("Node", node), showlegend = FALSE,
                hoverinfo = 'none')
    
    # Mittelwertslinie hinzufügen
    fig <- fig %>%
      add_trace(x = node_data$shifted_timesteps, y = node_data$Wcont_corr_mean, type = 'scatter', mode = 'lines+markers',
                name = paste("Node", node), legendgroup = paste("Node", node),
                line = list(width = 2, color = node_color),
                marker = list(size = 6, color = node_color),
                hoverinfo = 'text',
                text = paste("Nodes:", node, "<br>Mean:", round(node_data$Wcont_corr_mean, 2), "<br>SD:", round(node_data$Wcont_corr_sd, 2), "<br>Timesteps:", node_data$Timesteps))
    
    # SD-Linien hinzufügen
    fig <- fig %>%
      add_trace(x = node_data$shifted_timesteps, 
                y = node_data$Wcont_corr_mean - node_data$Wcont_corr_sd * scale_factor,
                type = 'scatter', mode = 'lines', line = list(dash = 'dot', width = 1, color = node_color),
                name = paste("Node", node, "SD"), legendgroup = paste("Node", node), showlegend = FALSE, hoverinfo = 'none') %>%
      add_trace(x = node_data$shifted_timesteps, 
                y = node_data$Wcont_corr_mean + node_data$Wcont_corr_sd * scale_factor,
                type = 'scatter', mode = 'lines', line = list(dash = 'dot', width = 1, color = node_color),
                name = paste("Node", node, "SD"), legendgroup = paste("Node", node), showlegend = FALSE, hoverinfo = 'none')
  }
  
  # Layout anpassen
  fig <- fig %>%
    layout(title = paste("Regime", regime, ": Mean Correlation over Timesteps"),
           xaxis = list(title = "Timesteps"),
           yaxis = list(title = "Wcont mean correlation"),
           legend = list(title = list(text = "Nodes")))
  
  # Plot anzeigen
  fig
}

htmlwidgets::saveWidget(fig, "Plots/Wcont_cor_lineplots.html")





library(plotly)
library(dplyr)
library(RColorBrewer)

# Skalierungsfaktor für die Breite der Ribbons
scale_factor <- 0.2  # Passe diesen Wert nach Bedarf an
x_offset <- 20  # Betrag für die Verschiebung der x-Werte

# Schleife über jeden Wert von Regimes
for (regime in unique(MSAR_results$Regimes)) {
  regime_data <- MSAR_results %>% filter(Regimes == regime)
  
  # Füge eine kleine Verschiebung zu den Timesteps basierend auf Nodes hinzu
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = Timesteps + (as.numeric(factor(Nodes)) - 1) * x_offset)
  
  # Definiere eine Farbpalette basierend auf der Anzahl der Nodes
  node_colors <- RColorBrewer::brewer.pal(length(unique(regime_data$Nodes)), "Set1")
  
  # Erstelle das Plotly-Objekt
  fig <- plot_ly()
  
  # Füge die Daten für jede Node-Gruppe hinzu
  for (i in seq_along(unique(regime_data$Nodes))) {
    node <- unique(regime_data$Nodes)[i]
    node_data <- regime_data %>% filter(Nodes == node)
    
    # Bestimme die Farbe für den aktuellen Node
    node_color <- node_colors[i]
    
    # Ribbon darstellen
    fig <- fig %>%
      add_trace(x = c(node_data$shifted_timesteps, rev(node_data$shifted_timesteps)), 
                y = c(node_data$Wcont_corr_mean - node_data$Wcont_corr_sd * scale_factor, 
                      rev(node_data$Wcont_corr_mean + node_data$Wcont_corr_sd * scale_factor)),
                type = 'scattergl', mode = 'lines',
                fill = 'toself', fillcolor = sprintf('rgba(%s,0.2)', paste(col2rgb(node_color), collapse = ',')),
                line = list(color = 'transparent', width = 0),
                name = paste("Node", node), legendgroup = paste("Node", node), showlegend = FALSE,
                hoverinfo = 'none')
    
    # Mittelwertslinie hinzufügen
    fig <- fig %>%
      add_trace(x = node_data$shifted_timesteps, y = node_data$Wcont_corr_mean, type = 'scattergl', mode = 'lines+markers',
                name = paste("Node", node), legendgroup = paste("Node", node),
                line = list(width = 2, color = node_color),
                marker = list(size = 6, color = node_color),
                hoverinfo = 'text',
                text = paste("Nodes:", node, "<br>Mean:", round(node_data$Wcont_corr_mean, 2), "<br>SD:", round(node_data$Wcont_corr_sd, 2), "<br>Timesteps:", node_data$Timesteps))
    
    # SD-Linien hinzufügen
    fig <- fig %>%
      add_segments(x = node_data$shifted_timesteps, 
                   xend = node_data$shifted_timesteps, 
                   y = node_data$Wcont_corr_mean - node_data$Wcont_corr_sd * scale_factor,
                   yend = node_data$Wcont_corr_mean + node_data$Wcont_corr_sd * scale_factor,
                   type = 'scattergl', mode = 'lines',
                   line = list(dash = 'dot', width = 1, color = node_color),
                   name = paste("Node", node, "SD"), legendgroup = paste("Node", node), showlegend = FALSE,
                   hoverinfo = 'none')
  }
  
  # Layout anpassen und Titel fett machen
  fig <- fig %>%
    layout(title = list(text = paste(regime, "Regimes: Mean Correlation over Timesteps"),
                        font = list(family = "Arial", size = 24, color = "black")),
           xaxis = list(title = "Timesteps"),
           yaxis = list(title = "Wcont mean correlation"),
           legend = list(title = list(text = "Nodes")))
  
  # Speichere das Plot für das aktuelle Regime als HTML-Datei
  file_name <- paste0("Plots/Wcont_cor_lineplot_", regime, "_regimes.html")
  htmlwidgets::saveWidget(fig, file_name)
  
 }






## Make lineplots for Wtemp_corr and Wcont_corr for each set of regimes
library(plotly)
library(dplyr)
library(RColorBrewer)

# Function to create lineplots with ribbons, markers and errorbars
create_ribbon_plot <- function(df, mean_col, sd_col, regime, scale_factor, x_offset) {
  regime_data <- df %>% filter(Regimes == regime)
  
  # Offset of x-values for nodes
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = Timesteps + (as.numeric(factor(Nodes)) - 1) * x_offset)
  
  # Set color for each node
  node_colors <- RColorBrewer::brewer.pal(length(unique(regime_data$Nodes)), "Set1")
  
  # Initialize plotly object
  fig <- plot_ly()
  
  # Iterate through set of nodes
  for (i in seq_along(N)) {
    node <- N[i]
    node_data <- regime_data %>% filter(Nodes == node)
    
    # Set node color
    node_color <- node_colors[i]
    
    # Plot ribbon
    fig <- fig %>%
      add_trace(x = c(node_data$shifted_timesteps, rev(node_data$shifted_timesteps)), 
                y = c(node_data[[mean_col]] - node_data[[sd_col]] * scale_factor, 
                      rev(node_data[[mean_col]] + node_data[[sd_col]] * scale_factor)),
                type = 'scattergl', mode = 'lines',
                fill = 'toself', fillcolor = sprintf('rgba(%s,0.2)', paste(col2rgb(node_color), collapse = ',')),
                line = list(color = 'transparent', width = 0),
                name = paste(node, "Nodes" ), legendgroup = paste("Node", node), showlegend = FALSE,
                hoverinfo = 'none')
  
    # Plot meanline
    fig <- fig %>%
      add_trace(x = node_data$shifted_timesteps, y = node_data[[mean_col]], type = 'scattergl', mode = 'lines+markers',
                name = paste(node, "Nodes"), legendgroup = paste("Node", node),
                line = list(width = 2, color = node_color),
                marker = list(size = 6, color = node_color),
                hoverinfo = 'text',
                text = paste("Nodes:", node, "<br>Mean:", round(node_data[[mean_col]], 2), "<br>SD:", round(node_data[[sd_col]], 2), "<br>Timesteps:", node_data$Timesteps))
  
    # Plot errorbars
    fig <- fig %>%
      add_segments(x = node_data$shifted_timesteps, 
                   xend = node_data$shifted_timesteps, 
                   y = node_data[[mean_col]] - node_data[[sd_col]] * scale_factor,
                   yend = node_data[[mean_col]] + node_data[[sd_col]] * scale_factor,
                   type = 'scattergl', mode = 'lines',
                   line = list(dash = 'dot', width = 1, color = node_color),
                   name = paste("Node", node, "SD"), legendgroup = paste("Node", node), showlegend = FALSE,
                   hoverinfo = 'none')
  }
  
  
  # Define layout
  fig <- fig %>%
    layout(title = list(text = paste(sub("_corr_mean", " - mean correlations per timesteps: ", mean_col), regime, " regimes"),
                        font = list(family = "Arial", size = 18, color = "black")),
           xaxis = list(title = "Timesteps"),
           yaxis = list(title = paste(mean_col, "correlation")))#,
          # legend = list(title = list(text = "Nodes")))
  
  return(fig)
}



# Set scalingfactor and nodewise offset
scale_factor <- 0.2
x_offset <- 50

# Generate Wtemp plots for each no. of regimes 
for (regime in unique(MSAR_results$Regimes)) {
  fig_temp <- create_ribbon_plot(MSAR_results, mean_col = "Wtemp_corr_mean", sd_col = "Wtemp_corr_sd", regime, scale_factor, x_offset)
  file_name_temp <- paste0("Plots/Wtemp_corr_lineplot_", regime, "_Regimes.html")
  htmlwidgets::saveWidget(fig_temp, file_name_temp)
}

# Generate Wcont plots for each no. of regimes 
for (regime in unique(MSAR_results$Regimes)) {
  fig_cont <- create_ribbon_plot(MSAR_results, mean_col = "Wcont_corr_mean", sd_col = "Wcont_corr_sd", regime, scale_factor, x_offset)
  file_name_cont <- paste0("Plots/Wcont_corr_lineplot_lineplot_", regime, "_Regimes.html")
  htmlwidgets::saveWidget(fig_cont, file_name_cont)
}