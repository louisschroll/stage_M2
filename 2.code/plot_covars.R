
contour_golfe2 <- st_crop(contour_golfe, xmin=691000, ymin=6129601, xmax=915000, ymax=6280000)

plot_list <- list(
  sst_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = mean_SST), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = "Température de l'eau en surface (SST)") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe2),
  
  depth_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = bathymetry), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = 'Bathymétrie') +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe2),
  
  dist_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = dist_to_shore), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = 'Distance à la côte') +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe2),
  
  chla_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = mean_CHL), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = 'Cholophylle A') +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe2),
  
  ssh_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = mean_SSH), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    labs(title = "Anomalie d'hauteur de l'eau (SSH)") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe2),
  
  conc_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = concavity), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    labs(title = "Concavité du fond marin") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe2)
)


library(cowplot)
plot_grid(plotlist = plot_list, ncol = 2)
