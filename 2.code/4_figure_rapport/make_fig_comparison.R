
make_fig_comparison <- function(species, fig_title = "", best_cov){
  grid_occ <- spOccupancy_res_grid %>% 
    select(all_of(c("grid_c", paste0(c("mean_psi_", "sd_psi_"), species))))
  names(grid_occ) <- c("mean_psi", "sd_psi", "grid_c")
  
  predictive_maps <- wrap_elements(plot_maps_by_models(species = species, 
                                         grid_occ = grid_occ, 
                                         grid_nmix = grid_nmix, 
                                         grid_RSF = grid_RSF, 
                                         grid_int = grid_int))

  coefficient_plot <- gather_coeff_values(sampleNmixture = sampleNmixture, 
                                          samples_spOcc = samples_spOcc, 
                                          samplesRSF = samplesRSF, 
                                          samplesint = samplesint,
                                          selected_cov = best_cov[[species]]) %>%  
    change_covariate_names() %>% 
    plot_coeff_by_model2(ncol = 2) 
    wrap_elements()
  
  precision_plot <- wrap_elements(plot_space_use_precision(grid_nmix = grid_nmix, 
                                             grid_RSF = grid_RSF, 
                                             grid_int = grid_int, 
                                             grid_spOcc = grid_occ))
  
  plot_left <- (coefficient_plot / precision_plot) +
    plot_layout(heights = c(7, 4))
  
  plot_final <- (predictive_maps | plot_left) +
    plot_layout(
      guide = "keep"
    ) +
    plot_annotation(
      title = fig_title,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.position = 'top'
      )
    )
  
  return(plot_final)
}



plot_maps_by_models <- function(species, grid_occ, grid_nmix, grid_RSF, grid_int){
  # Plot all the maps
  # 0/3 spOccupancy 
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
                               plot_title = "RSF + N-mixture")
  
  predictive_maps <- 
    (plots_occ$mean_psi_plot + plots_occ$sd_psi_plot) /
    (plots_nmix$mean_psi_plot + plots_nmix$sd_psi_plot) /
    (plots_rsf$mean_psi_plot + plots_rsf$sd_psi_plot) /
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



plot_maps_by_models_2col <- function(species, grid_occ, grid_nmix, grid_RSF, grid_int, title=""){
  # Plot all the maps
  # 0/3 spOccupancy 
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
                               plot_title = "RSF + N-mixture")
  
  predictive_maps <- 
    ((plots_occ$mean_psi_plot + plots_occ$sd_psi_plot) |
    (plots_nmix$mean_psi_plot + plots_nmix$sd_psi_plot)) /
    ((plots_rsf$mean_psi_plot + plots_rsf$sd_psi_plot) |
    (plots_int$mean_psi_plot + plots_int$sd_psi_plot)) +
    plot_layout(guides = "collect", nrow = 2) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          family = "Calibri"),
        legend.position = "bottom", 
        legend.spacing.x = unit(3, "cm")
      )
    ) 
  
  return(predictive_maps)
}
