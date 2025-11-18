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




# #### Make 3D surface plots for each density
# #Get z-values for plot
# global_min <- Inf
# global_max <- -Inf
# reg_val <- unique(MSAR_results$Regimes)
# dens_val <- unique(MSAR_results$Density)
# T_val <- unique(MSAR_results$Timesteps)
# N_val <- unique(MSAR_results$Nodes)
# 
# z_matrices <- list()
# # Create matrix for z-values for each regime
# for (regime in seq_along(reg_val)) {
#   z_mat1 <- list()
#   df_r_subset <- subset(MSAR_results, Regimes == reg_val[regime])
#   for (density in seq_along(dens_val)) {
#     df_dens_subset <- subset(df_r_subset, Density == dens_val[density])
#     # Initialize z-matrix
#     z_matrix <- matrix(NA, nrow = length(N_val), ncol = length(T_val))
#     rownames(z_matrix) <- N_val
#     colnames(z_matrix) <- T_val
#   
# 
#   # Fill z-matrix with mean values from Wtemp_corr_mean
# 
#   for (nodes in seq_along(N_val)) {
#     for (timesteps in seq_along(T_val)) {
#       value <- df_dens_subset$Wtemp_corr_mean[df_dens_subset$Nodes == N_val[nodes] & df_dens_subset$Timesteps == T_val[timesteps]]
#       if (length(value) > 0) {
#         z_matrix[nodes,timesteps] <- value
#         global_min <- min(global_min, value, na.rm = TRUE)
#         global_max <- max(global_max, value, na.rm = TRUE)
#         }
#       }
#     }
#   z_mat1[[paste0(dens_val[density], "_Density")]] <- z_matrix
#   }
#   z_matrices[[paste0(reg_val[regime], "_Regimes")]] <- z_mat1
# }
# 
# ## Create 3D surface plots
# # Initialize plot
# fig <- plot_ly()
# 
# # Set plot colours
# colorscale <- list(c(0, "blue"),
#                    c(0.33, "cyan"),
#                    c(0.66, "yellow"),
#                    c(1, "red")
#                    )
# 
# # Füge Oberflächen für jeden Wert in Regimes hinzu
# for(density in seq_along(dens_val)) {
#   fig <- plot_ly()
#   annotations <- list()
#   shapes <- list()
# 
#   for (regime in seq_along(reg_val)) {
#     z <- z_matrices[[regime]][[density]]
#     # Add surface for each regime to plot
#     fig <- fig %>% add_surface(
#       x = T_val,
#       y = N_val,
#       z = z,
#       name = paste("Regimes", reg_val[regime]),
#       opacity = 1,
#       colorscale = colorscale,
#       zmin = global_min,
#       zmax = global_max,
#       showscale = regime == 1
#       )
#     
#     
#     # Set coordinates for annotations
#     x_annotation <- mean(T_val)
#     y_annotation <- mean(N_val)
#     z_annotation <- global_max
#     
#     # Add annotations
#     annotations <- append(annotations, list(
#       list(
#         x = x_annotation,
#         y = y_annotation,
#         z = z_annotation,
#         text = paste("Regimes", regime),
#         xanchor = "center",
#         yanchor = "bottom",
#         showarrow = TRUE,
#         arrowhead = 2,
#         ax = 0,
#         ay = -40
#       )
#     ))
#     
#     # Add outlines
#     shapes <- append(shapes, list(
#       list(
#         type = "line",
#         x0 = x_annotation,
#         y0 = y_annotation,
#         z0 = z_annotation,
#         x1 = x_annotation,
#         y1 = y_annotation,
#         z1 = z_annotation - 0.1 * (global_max - global_min),
#         line = list(color = "black", width = 2)
#       )
#     ))
#     }
#   
#   # Set layout for plot
#   fig <- fig %>% layout(
#     title = paste0("3D surface plot for density", dens_val[density]),
#     scene = list(
#       xaxis = list(title = "Timesteps", tickvals = T_val),
#       yaxis = list(title = "Nodes", tickvals = N_val),
#       zaxis = list(title = "Mean")
#     ),
#     annotations = annotations,
#     shapes = shapes
#   )
#   # Save combined plot as html file
#   filename <- paste0("Plots/surface_plot_density", dens_val[density], ".html")
#   htmlwidgets::saveWidget(fig, filename)
# }
# 

# 
# ### Make separate line plot panels for each density, with 4 individual plots for each no. regimes 
# # Specify offset for nodes and scaling factor for ribbons
# scale_factor <- 0.2
# x_offset <- 10
# 
# ## Function to make line plots
# library(plotly)
# library(dplyr)
# library(RColorBrewer)
# 
# # Function to make lineplots
# line_plot <- function(df, mean_col, sd_col, regime, scale_factor, x_offset, y_limits, show_legend) {
#   regime_data <- df %>% filter(Regimes == regime)
#   
#   # Set x-offset for each no. of nodes
#   regime_data <- regime_data %>%
#     mutate(shifted_timesteps = Timesteps + (as.numeric(factor(Nodes)) - 1) * x_offset)
#   
#   # Set nodes color palette
#   node_colors <- RColorBrewer::brewer.pal(length(unique(regime_data$Nodes)), "Set1")
#   
#   # Initialise plotly object
#   fig <- plot_ly()
#   
#   # Loop over nodes
#   N <- unique(regime_data$Nodes)
#    for (i in seq_along(N)) {
#     node <- N[i]
#     node_data <- regime_data %>% filter(Nodes == node)
#     
#     # Set node color
#     node_color <- node_colors[node]
#     
#     # Plot line graph
#     fig <- fig %>%
#       # Add ribbon for sd
#       add_trace(x = c(node_data$shifted_timesteps,
#                       rev(node_data$shifted_timesteps)), 
#                 y = c(node_data[[mean_col]] - node_data[[sd_col]] * scale_factor, 
#                       rev(node_data[[mean_col]] + node_data[[sd_col]] * scale_factor)),
#                 type = 'scatter',
#                 mode = 'lines',
#                 fill = 'toself',
#                 fillcolor = sprintf('rgba(%s,0.2)', paste(col2rgb(node_color), collapse = ',')),
#                 line = list(color = 'transparent'),
#                 name = paste(node, "Nodes"),
#                 legendgroup = paste("Node", node),
#                 showlegend = FALSE,
#                 hoverinfo = 'none'
#                 )
#     
#                 
#     
#     # Add meanline
#     fig <- fig %>%
#       add_trace(x = node_data$shifted_timesteps,
#                 y = node_data[[mean_col]],
#                 type = 'scatter',
#                 mode = 'lines+markers',
#                 name = paste(node, "Nodes"),
#                 legendgroup = paste("Node", node),
#                 showlegend = show_legend,
#                 line = list(width = 2, color = node_color),
#                 marker = list(size = 6, color = node_color),
#                 hoverinfo = 'text',
#                 text = paste("Nodes:", node,
#                              "<br>Mean:", round(node_data[[mean_col]], 2),
#                              "<br>SD:", round(node_data[[sd_col]], 2),
#                              "<br>Timesteps:", node_data$Timesteps)
#                 )
#     
#     # Plot errorbars
#     fig <- fig %>%
#       add_segments(x = node_data$shifted_timesteps,
#                    xend = node_data$shifted_timesteps,
#                    y = node_data[[mean_col]] - node_data[[sd_col]] * scale_factor,
#                    yend = node_data[[mean_col]] + node_data[[sd_col]] * scale_factor,
#                    type = 'scatter',
#                    mode = 'lines',
#                    line = list(dash = 'dot', width = 1, color = node_color),
#                    name = paste("Node", node, "SD"),
#                    legendgroup = paste("Node", node),
#                    showlegend = FALSE,
#                    hoverinfo = 'none')
#   }
#   
#   # Define layout without title
#   fig <- fig %>%
#     layout(
#       xaxis = list(title = "Timesteps"),
#       yaxis = list(title = paste(mean_col, "Correlation"), range = y_limits)
#     )
#   
#   return(fig)
# }
# 
# 
# # Create dir 
# dir.create("Plots", showWarnings = FALSE)
# 
# for (result in seq(6, ncol(MSAR_results) -1, by = 2)) {
#   results <- data.frame()
#   if (result+1 <= ncol(MSAR_results)) {
#     results <- cbind(
#       MSAR_results[ ,1:5],
#       setNames(MSAR_results[, result, drop = FALSE], colnames(MSAR_results)[result]),
#       setNames(MSAR_results[, result + 1, drop = FALSE], colnames(MSAR_results)[result + 1])
#                )
#     }
#   mean_col <- colnames(results)[ncol(results) - 1]
#  sd_col <- colnames(results)[ncol(results)]
# 
#   
# # Set y-axis limits
#   means <- results[ ,ncol(results) -1]
#   sd <- results[ ,ncol(results)]
#   y_limits <- range(c(
#     means - sd * scale_factor,
#     means + sd * scale_factor
#     ),
#     na.rm = TRUE)
# 
# # Loop over density
# for (density in unique(results$Density)) {
#   density_data <- results %>% filter(Density == density)
#   
#   # Get values for each no. of regimes
#   regime_values <- unique(density_data$Regimes)
#   subplot_list <- list()
# 
#   for (regime in regime_values) {
#     # Show legend only once
#     show_legend <- regime == min(regime_values)
#     
#     fig_temp <- line_plot(
#       density_data, 
#       mean_col, 
#       sd_col, 
#       regime, 
#       scale_factor, 
#       x_offset, 
#       y_limits, 
#       show_legend = show_legend
#     )
#     subplot_list <- append(subplot_list, list(fig_temp))
#   }
# 
#   
#   # Combine plots in panel
#   fig_panel <- subplot(
#     subplot_list, 
#     nrows = 2, 
#     titleX = TRUE, 
#     titleY = TRUE, 
#     margin = 0.05,
#     shareY = TRUE,
#     shareX = TRUE
#   )
#   
#   ## Add subplot titles
#   # Set position of subplots
#   annotation_positions <- list(
#     list(x = 0.225, y = 1.0),
#     list(x = 0.775, y = 1.0),
#     list(x = 0.225, y = 0.5),
#     list(x = 0.775, y = 0.5)
#   )
#   
#   # Make annotations for subplot titles
#   annotations <- list()
#   for (i in seq_along(regime_values)) {
#     annotations[[i]] <- list(
#       text = paste(regime_values[i], " Regimes"),
#       x = annotation_positions[[i]]$x,
#       y = annotation_positions[[i]]$y,
#       xref = "paper",
#       yref = "paper",
#       xanchor = "center",
#       yanchor = "bottom",
#       showarrow = FALSE,
#       font = list(size = 14)
#     )
#   }
#   
#   # Adjust panel layout
#   fig_panel <- fig_panel %>%
#     layout(
#       #title = paste("Density", density, " %"),
#       annotations = c(
#         annotations,
#         list(
#           list(
#             text = paste("Density", density, " %"),
#             x = 0.5,
#             y = 0,
#             xref = "paper",
#             yref = "paper",
#             showarrow = FALSE,
#             font = list(size = 16)
#           )
#         )
#       )
#     )
#   
#   # Storing panel
#   file_name_panel <- paste0("Plots/Panel_Density_", density, ".html")
#   htmlwidgets::saveWidget(fig_panel, file_name_panel)
#   }
# }
# 




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
