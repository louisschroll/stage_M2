#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/plot_prediction.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-03
#'
#' Script description:
#' This function makes two plots of the study area representing the mean intensity
#' of space use and the standard deviation associated with the predicted values, and 
#' return them in a list.
#' Mandatory arguments:
#' @new_grid: a sf object, with at least two columns called mean_psi and sd_psi
#' Facultative arguments:
#' @add_colony: a boolean indicating whether colony stored in countLL file should 
#' be plot on the map
#' @species_colony: if add_colony is true, species_colony is the name of the species
#' for which we want to plot colony
#' -------------------------------------------------------------------------------

load("~/stage_M2/1.data/contour_golfe_du_lion.rdata")
load("~/stage_M2/1.data/countLL.rdata")

plot_prediction <- function(new_grid, add_colonies = F, species_colony=NULL, 
                            legend_position="bottom", plot_title = ""){
  # Plot the intensity of space use
  mean_psi_plot <- ggplot() +
    geom_sf(data = new_grid, aes(fill = mean_psi), color = NA) +
    scale_fill_distiller(palette = "Spectral", 
                         breaks = c(min(new_grid$mean_psi), max(new_grid$mean_psi)),
                         labels = c("Low", "High")) +
    geom_sf(data = contour_golfe,
            color = "black",
            fill = "antiquewhite") +
    labs(title = plot_title) +
    theme_bw() +
    theme(
      legend.position = legend_position,
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.title.position = 'plot',
      ) +
    guides(
      fill = guide_colourbar(
        title = "Relative space-use",
        title.theme = element_text(
          #family = "Helvetica",
          face = "bold",
          size = 16,
          hjust = 0.5),
        title.position = 'top',
        ticks = FALSE,
        title.hjust = .5,
        barwidth = unit(10, 'lines'),
        barheight = unit(.5, 'lines')),
      colour = "none")
  
  if (add_colonies) {
    mean_psi_plot <- mean_psi_plot +
      geom_sf(data = countLL %>% filter(Species == species_colony) %>% st_crop(st_bbox(contour_golfe)),
              pch = 16, size = 2)
  }
  
  # Plot uncertainty in the prediction
  sd_psi_plot <- ggplot() + 
    geom_sf(data = new_grid, aes(fill = sd_psi), color = NA) +
    scale_fill_viridis_c(option = "inferno",
                         breaks = c(min(new_grid$sd_psi), max(new_grid$sd_psi)),
                         labels = c("Low", "High")
                         ) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = plot_title) +
    theme_bw() +
    theme(legend.position = legend_position, 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    guides(
      fill = guide_colourbar(
        title = "Standard deviation",
        title.theme = element_text(
          #family = "Helvetica",
          face = "bold",
          size = 16,
          hjust = 0.5),
        title.position = 'top',
        ticks = FALSE,
        title.hjust = .5,
        barwidth = unit(10, 'lines'),
        barheight = unit(.5, 'lines')),
      colour = "none")
  
  return(list(mean_psi_plot = mean_psi_plot, sd_psi_plot = sd_psi_plot))
}

library(extrafont)
blue <-  "#3d405b"
orange <- "#f2cc8f"
red <- "#e07a5f"
purple <-  "#81b29a"

pred_precision <- tibble(value = new_grid_nmix$sd_psi, model = "Poisson GLM") %>%
  bind_rows(tibble(value = new_grid_RSF$sd_psi, model = "RSF")) %>%
  bind_rows(tibble(value = new_grid_int$sd_psi, model = "Integrated model")) %>%
  ggplot(aes(x = value, y = model, color = as.factor(model))) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)
  ) +
  geom_boxplot(width = .25,
               outlier.shape = NA) +
  # geom_point(
  #   size = 1.3,
  #   alpha = .01,
  #   aes(color = model),
  #   position = position_jitter(
  #     seed = 1, width = .1
  #   )) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Space-use precision",
    subtitle = "Predicted CV among grid-cells",
    x = "Coefficient of variation (CV)",
    y = ""
  ) +
  scale_color_manual(
    values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "Poisson GLM" = blue
    )
  ) +
  scale_fill_manual(
    values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "Poisson GLM" = blue
    )
  ) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.text.y = element_text(colour = 'black', size = 8),
    plot.title = element_text(
      size = 16,
      face = "bold",
      family  = "Helvetica",
      hjust = 0.5
    ),
    plot.subtitle = element_text(family  = "Helvetica", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  guides(color = "none", 
         fill = "none")

