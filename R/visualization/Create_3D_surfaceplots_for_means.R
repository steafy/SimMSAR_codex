library(plotly)


## Extract stats from results list to dataframe
# Initialize list to store stats data for each combination from Density*T*N*M
data_list <- list()
MSAR_results <- list()

# Extract stats data from MSAR results list
for (t in names(MSAR_models)) {
  for (density in names(MSAR_models[[t]])) {
    for (nodes in names(MSAR_models[[t]][[density]])) {
      for (regimes in names(MSAR_models[[t]][[density]][[nodes]])) {
        stats <- MSAR_models[[t]][[density]][[nodes]][[regimes]][["Stats"]]
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


### Erstellen Sie separate Linienplot-Panels für jede Dichte mit 4 individuellen Plots für jede Anzahl von Regimen
# Skalierungsfaktor und Offset für die Bänder festlegen
scale_factor <- 1
x_offset <- 20

## Funktion zum Erstellen von Linienplots
library(plotly)
library(dplyr)
library(RColorBrewer)

# Funktion zum Erstellen von Linienplots
line_plot <- function(df, mean_col, sd_col, regime, scale_factor, x_offset, y_limits, show_legend) {
  regime_data <- df %>% filter(Regimes == regime)
  
  # x-Offset für jede Anzahl von Knoten festlegen
  regime_data <- regime_data %>%
    mutate(shifted_timesteps = Timesteps + (as.numeric(factor(Nodes)) - 1) * x_offset)
  
  # Farbschema für Knoten festlegen
  node_colors <- RColorBrewer::brewer.pal(length(unique(regime_data$Nodes)), "Set1")
  
  # Plotly-Objekt initialisieren
  fig <- plot_ly()
  
  # Schleife über Knoten
  N <- unique(regime_data$Nodes)
  for (i in seq_along(N)) {
    node <- N[i]
    node_data <- regime_data %>% filter(Nodes == node)
    
    # Farbe für Knoten festlegen
    node_color <- node_colors[i]
    
    # Liniengraph zeichnen
    fig <- fig %>%
      # Band für Standardabweichung hinzufügen
      add_trace(
        x = c(node_data$shifted_timesteps, rev(node_data$shifted_timesteps)), 
        y = c(node_data[[mean_col]] - node_data[[sd_col]] * scale_factor, 
              rev(node_data[[mean_col]] + node_data[[sd_col]] * scale_factor)),
        type = 'scatter',
        mode = 'lines',
        fill = 'toself',
        fillcolor = sprintf('rgba(%s,0.2)', paste(col2rgb(node_color), collapse = ',')),
        line = list(color = 'transparent'),
        name = paste(node, "Nodes"),
        legendgroup = paste("Node", node),
        showlegend = FALSE,
        hoverinfo = 'none'
      )
    
    # Mittelwertlinie hinzufügen
    fig <- fig %>%
      add_trace(
        x = node_data$shifted_timesteps,
        y = node_data[[mean_col]],
        type = 'scatter',
        mode = 'lines+markers',
        name = paste(node, "Nodes"),
        legendgroup = paste("Node", node),
        showlegend = show_legend,
        line = list(width = 2, color = node_color),
        marker = list(size = 6, color = node_color),
        hoverinfo = 'text',
        text = paste("Nodes:", node,
                     "<br>Mean:", round(node_data[[mean_col]], 2),
                     "<br>SD:", round(node_data[[sd_col]], 2),
                     "<br>Timesteps:", node_data$Timesteps)
      )
    
    # Fehlerbalken hinzufügen
    fig <- fig %>%
      add_segments(
        x = node_data$shifted_timesteps,
        xend = node_data$shifted_timesteps,
        y = node_data[[mean_col]] - node_data[[sd_col]] * scale_factor,
        yend = node_data[[mean_col]] + node_data[[sd_col]] * scale_factor,
        type = 'scatter',
        mode = 'lines',
        line = list(dash = 'dot', width = 1, color = node_color),
        name = paste("Node", node, "SD"),
        legendgroup = paste("Node", node),
        showlegend = FALSE,
        hoverinfo = 'none'
      )
  }
  
  # Layout definieren ohne Titel
  fig <- fig %>%
    layout(
      xaxis = list(
        title = list(
          text = "Timesteps",
          standoff = 8,
          font = list(size = 11)
        ),
        tickfont = list(size = 10)#,
        # showline = FALSE,
        # linewidth = 1,
        # linecolor = "black"
      ),
      yaxis = list(
        title = list(
          text = paste("Mean correlation"),
          standoff = 8,
          font = list(size = 11)
        ),
        range = y_limits,
        tickfont = list(size = 10)#,
        # showline = TRUE,
        # linewidth = 1,
        # linecolor = "black"
      ),
      margin = list(
        l = 10,
        r = 10,
        b = 40,
        t = 50  # Oberer Rand für Subplot-Titel
      )
    )
  
  return(fig)
}

dir.create("Plots", showWarnings = FALSE)

for (result in seq(6, ncol(MSAR_results) - 1, by = 2)) {
  if (result + 1 <= ncol(MSAR_results)) {
    results <- cbind(
      MSAR_results[, 1:5],
      setNames(MSAR_results[, result, drop = FALSE], colnames(MSAR_results)[result]),
      setNames(MSAR_results[, result + 1, drop = FALSE], colnames(MSAR_results)[result + 1])
    )
    
    # Get names for the respective results
    mean_col <- colnames(MSAR_results)[result]
    sd_col <- colnames(MSAR_results)[result + 1]
    
    # Set limits for y-axis
    means <- results[[mean_col]]
    sd <- results[[sd_col]]
    y_limits <- range(
      c(means - sd * scale_factor,
        means + sd * scale_factor
        ),
      na.rm = TRUE
    )
    
    # Loop over density
    for (density in unique(results$Density)) {
      density_data <- results %>% filter(Density == density)
      
      # Get values for each no. of regimes
      regime_values <- unique(density_data$Regimes)
      subplot_list <- list()
      
      for (regime in regime_values) {
        # Add legend only once
        show_legend <- regime == min(regime_values)
        
        fig_temp <- line_plot(
          density_data, 
          mean_col, 
          sd_col, 
          regime, 
          scale_factor, 
          x_offset, 
          y_limits, 
          show_legend = show_legend
        )
        subplot_list <- append(subplot_list, list(fig_temp))
      }
      
      # Combine subplots
      fig_panel <- subplot(
        subplot_list, 
        nrows = 2, 
        titleX = TRUE, 
        titleY = TRUE,
        margin = 0.075,
        #margin = c(0.02, 0.02, 0.02, 0.02),
        shareY = FALSE,
        shareX = FALSE
      )
      
      # Extract subplot domains
      layout_fig <- fig_panel$x$layout
      x_domains <- list()
      y_domains <- list()
      for (i in seq_along(regime_values)) {
        x_axis_name <- paste0("xaxis", ifelse(i == 1, "", i))
        y_axis_name <- paste0("yaxis", ifelse(i == 1, "", i))
        x_domains[[i]] <- layout_fig[[x_axis_name]]$domain
        y_domains[[i]] <- layout_fig[[y_axis_name]]$domain
      }
      
      # Add annotations
      annotations <- list()
      for (i in seq_along(regime_values)) {
        x_domain <- x_domains[[i]]
        y_domain <- y_domains[[i]]
        x_mid <- mean(x_domain)
        y_top <- y_domain[2]
        
        annotations[[i]] <- list(
          text = paste(regime_values[i], "Regimes"),
          x = x_mid,
          y = y_top + 0.02, 
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 14)
        )
      }
      if (mean_col == "Wtemp_corr_mean") {header <- "Mean correlations for Wtemp"}
      if (mean_col == "Wcont_corr_mean") {header <- "Mean correlations for Wcont"}
      if (mean_col == "Wtemp_ac_corr_mean") {header <- "Mean correlations for average controlablilty of Wtemp"}
      if (mean_col == "Wcont_ac_corr_mean") {header <- "Mean correlations for average controlablilty of Wtemp"}
      
      
      # Adjust panel layout
      fig_panel <- fig_panel %>%
        layout(
          annotations = c(
            annotations,
            list(
              list(
                text = paste0(header, "<br>Density ", density, "%"),
                x = 0.5,
                y = 1.15,  
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE,
                font = list(size = 16)
              )
            )
          ),
          margin = list(
            l = 50,
            r = 30,
            b = 50,
            t = 150
          )
        )
      
      # Store panel
      file_name_panel <- paste0("Plots/Panel_Density_", density, "_", mean_col, ".html")
      htmlwidgets::saveWidget(fig_panel, file_name_panel)
    }
  }
}
