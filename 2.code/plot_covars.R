
library(sf)
library(tidyverse)

# load data
load("1.data/all_seabirds_counts.rdata")
load("1.data/covariates_data.rdata")
load("1.data/contour_golfe_du_lion.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

plot_list <- list(
  sst_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = mean_SST), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = "Température de l'eau en surface (SST)") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe),
  
  depth_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = bathymetry), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = 'Bathymétrie') +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe),
  
  dist_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = dist_to_shore), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = 'Distance à la côte') +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe),
  
  chla_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = mean_CHL), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    #scale_fill_viridis_c(option = "B") + 
    labs(title = 'Cholophylle A') +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe),
  
  ssh_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = mean_SSH), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    labs(title = "Anomalie d'hauteur de l'eau (SSH)") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe),
  
  conc_plot = ggplot() + 
    geom_sf(data = grid, aes(fill = concavity), lwd = 0.1) +
    scale_fill_distiller(palette = "Spectral") +
    labs(title = "Concavité du fond marin") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_sf(data = contour_golfe)
)


library(cowplot)
plot_grid(plotlist = plot_list, ncol = 2, 
          label = "Figure 2 : Variables environnementales utilisées")
