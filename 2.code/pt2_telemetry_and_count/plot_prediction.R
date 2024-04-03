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

plot_prediction <- function(new_grid, add_colonies = F, species_colony=NULL){
  # Plot the intensity of space use
  mean_psi_plot <- ggplot() + 
    geom_sf(data = new_grid, aes(fill = mean_psi), color = NA) +
    scale_fill_distiller(palette = "Spectral",
                         guide = guide_colorbar(ticks = FALSE,
                                                barwidth = 9,
                                                barheight = 0.7)) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = "Intensity of space use") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) 
  
  if (add_colonies){
    mean_psi_plot <- mean_psi_plot +
      geom_sf(data = countLL %>% filter(Species == species_colony) %>% st_crop(st_bbox(contour_golfe)),
              pch = 16, size = 2)
  }
  
  # Plot uncertainty in the prediction
  sd_psi_plot <- ggplot() + geom_sf(data = new_grid, aes(fill = sd_psi), color = NA) +
    scale_fill_viridis_c(option = "inferno", 
                         guide = guide_colorbar(ticks = FALSE,
                                                barwidth = 9,
                                                barheight = 0.7)) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = "SD intensity") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(list(mean_psi_plot = mean_psi_plot, sd_psi_plot = sd_psi_plot))
}

