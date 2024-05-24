
plot_maps_by_methods <- function(species){
  # Load predictions
  load(file = paste0(local_path, "3.results/prediction_grid/grid_", species, ".rdata"))
  
  # Plot all the maps
  # 0/3 spOccupancy 
  grid_occ <- spOccupancy_res_grid %>% 
    select(all_of("grid_c", paste0(c("mean_psi_", "sd_psi_"), species)))
  names(grid_occ) <- c("mean_psi", "sd_psi", "grid_c")
  plots_occ <- plot_prediction(grid_occ, 
                                       plot_title = "Occupancy")
  # 1/3 - N-mixture
  plots_nmix <- plot_prediction(grid_nmix, 
                                plot_title = "N-mixture")
  # 2/3 - RSF
  plots_rsf <- plot_prediction(grid_RSF, 
                               plot_title = "RSF")
  # 3/3 - RSF & N-mixture
  plots_int <- plot_prediction(grid_int, 
                               plot_title = "Integrated model")
  
  predictive_maps <- 
    (plots_occ$mean_psi_plot + plots_occ$sd_psi_plot) /
    (plots_rsf$mean_psi_plot + plots_rsf$sd_psi_plot) /
    (plots_nmix$mean_psi_plot + plots_nmix$sd_psi_plot) /
    (plots_int$mean_psi_plot + plots_int$sd_psi_plot) +
    plot_layout(guides = "collect", nrow = 4) +
    plot_annotation(
      theme = theme(
        legend.position = "bottom", 
        legend.spacing.x = unit(3, "cm")
      )
    )
  
  return(predictive_maps)
}