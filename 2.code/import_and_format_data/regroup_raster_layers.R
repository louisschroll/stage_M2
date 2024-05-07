# HEADER ------------------------------------------------------------------------
#
# Script name:  regroup_raster_layers.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-15
#
# Script description:
# The purpose of this script is to group all the data downloaded from different
# sources in one .tif file, with the same extent and precision for all raster
# layers. 
# The dynamic covariables are averaged across all years. 
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Import packages
library(terra)
library(tidyverse)

# Define variables used in the script
study_area <- ext(3, 5.6, 42.25, 43.7) # (xmin, xmax, ymin, ymax)


combine_ncfile <- function(path_non_int, path_interim){
  # Combine non-interim and interim periods
  ## Open the .nc files
  raster <- rast(path_non_int)
  raster_int <- rast(path_interim)
  
  ## Give a number to each month to rename the columns
  i_max <- nlyr(raster)
  i_max_int <- nlyr(raster_int)
  new_names <- as.character(1:i_max)
  new_names_int <- as.character(1:i_max_int + i_max)
  
  names(raster) <- new_names
  names(raster_int) <- new_names_int
  ## Combine the two rasters 
  raster_combine <- c(raster, raster_int)
  return(raster_combine)
}



# ---- Import Sea Surface Temperature (SST) and average by season -----
raster_sst <- combine_ncfile("1.data/copernicus_data/SST.nc", "1.data/copernicus_data/SST_int.nc") %>% 
  terra::crop(study_area) %>% 
  terra::extend(study_area) 

# Average seasonal SST across all years
mean_season <- raster_sst %>% as_tibble() %>% 
  bind_cols(crds(raster_sst)) %>% 
  pivot_longer(-c(x, y), names_to = "date_nb", values_to = "value") %>% 
  mutate(month_nb = (as.numeric(date_nb)-1)%%12 + 1) %>% 
  mutate(season_nb =  ceiling(month_nb/3)) %>% 
  group_by(x, y, season_nb) %>%  
  summarize(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE)
            ) %>% 
  ungroup()

# Plot the raster
if (FALSE) {
  mean_season %>% mutate(mean_value = scale(mean_value)) %>%
    ggplot(aes(x=x, y=y, fill=mean_value)) +
      geom_raster() +
      facet_wrap(~season_nb) +
      scale_fill_distiller(palette = "Spectral") 
}
# Create one col for each season and compute all season mean SST
raster_mean_SST <- mean_season %>% 
  pivot_wider(names_from = season_nb, values_from = c(mean_value, sd_value)) %>% 
  rename(mean_winter_SST = mean_value_1,
         mean_spring_SST = mean_value_2,
         mean_summer_SST = mean_value_3,
         mean_autumn_SST = mean_value_4,
         sd_winter_SST = sd_value_1,
         sd_spring_SST = sd_value_2,
         sd_summer_SST = sd_value_3,
         sd_autumn_SST = sd_value_4) %>% 
  #select(c(x, y) | starts_with("mean_")) %>% 
  terra::rast(type="xyz", crs="", digits=6, extent=NULL) 
 
all_SST <- c(app(raster_sst, mean),app(raster_sst, sd)) %>% 
  resample(raster_mean_SST)
names(all_SST) <- c("mean_SST", "sd_SST")
raster_mean_SST <- c(raster_mean_SST, all_SST)
plot(raster_mean_SST)
corrplot::corrplot.mixed(raster_mean_SST %>% as_tibble() %>% cor(method = "pearson", use = "na.or.complete"))

# ---- Import and crop static covariates data ----

raster_static_cov <- rast("1.data/static_covariates_raster.tif")
crs(raster_mean_SST) <- crs(raster_static_cov)

# Use resample to have same extent and resolution
raster_static_cov <- raster_static_cov %>% 
  terra::resample(raster_mean_SST) 

plot(raster_mean_SST[[1:4]])
plot(raster_static_cov)
# ---- Import and average the other dynamic covariates ----
# Salinity (SAL), euphotic zone depth (EZD) and Net primary productivity are
# excluded because their correlation whith CHL is too high
# Sea surface heigth (SSH) is too correlated with dist_to_shore and bathymetry
# And sd(chl) to correlated with mean(chl)

CHL <- combine_ncfile("1.data/copernicus_data/CHL.nc", "1.data/copernicus_data/CHL_int.nc") 
SAL <- combine_ncfile("1.data/copernicus_data/SAL.nc", "1.data/copernicus_data/SAL_int.nc") 
SSH <- combine_ncfile("1.data/copernicus_data/SSH.nc", "1.data/copernicus_data/SSH_int.nc") 

raster_dynamic_cov <- c(app(CHL, mean), app(CHL, sd),
  app(SAL, mean), app(SAL, sd),
  app(SSH, mean), app(SSH, sd)
  ) %>% 
  terra::crop(raster_mean_SST) %>% 
  terra::extend(raster_mean_SST) 

VEL <- combine_ncfile("1.data/copernicus_data/velocity.nc", "1.data/copernicus_data/velocity_int.nc") %>% 
  terra::resample(raster_mean_SST) 

# EZD <- rast("1.data/copernicus_data/euphotic_depth.nc")%>%
#   terra::resample(raster_mean_SST)
# NPP <- rast("1.data/copernicus_data/NPP.nc") %>%
#   terra::resample(raster_mean_SST)
raster_dynamic_cov <- c(raster_dynamic_cov, 
                        app(VEL, mean), app(VEL, sd)
                        # app(EZD, mean), app(EZD, sd),
                        # app(NPP, mean), app(NPP, sd)
                        )

names(raster_dynamic_cov) <- c("mean_CHL", "sd_CHL",
                               "mean_SAL", "sd_SAL",
                               "mean_SSH", "sd_SSH",
                               "mean_VEL", "sd_VEL"
                               # "mean_EZD", "sd_EZD",
                               # "mean_NPP", "sd_NPP"
                               )
plot(raster_dynamic_cov)
corrplot::corrplot.mixed(raster_dynamic_cov %>% as_tibble() %>% cor(method = "pearson", use = "na.or.complete"))

# ---- Import and average fishing effort data ----
# load(file = "fishing_effort.rdata")
# 
# hist(ais_tibble$apparent_fishing_hours, breaks = 1000)
# ggplot(data = ais_tibble, aes(x=lon, y=lat, fill=log(apparent_fishing_hours))) +
#   geom_raster() +
#   scale_fill_distiller(palette = "Spectral") +
#   facet_wrap(~year)
# ais_tibble$apparent_fishing_hours[ais_tibble$apparent_fishing_hours > 5] <- 5
# ais_tibble$apparent_fishing_hours[ais_tibble$apparent_fishing_hours == 0] <- NA

# ---- Combine all the rasters ----
combined_rasters <- c(raster_mean_SST, raster_static_cov, raster_dynamic_cov) 
plot(combined_rasters)

# cut the area too far offshore
triangle <- vect(matrix(c(3, 5.7, 3, 3, 42.1, 43.1, 46, 42.1), ncol = 2), type="polygons", atts=NULL, crs=crs(combined_rasters))
combined_rasters2 <- mask(x=combined_rasters, mask=triangle)
if (FALSE){
  load("1.data/contour_golfe_du_lion.rdata")
  combined_rasters2 %>% as_tibble() %>% 
    drop_na() %>% 
    bind_cols(crds(combined_rasters2)) %>% 
    ggplot() +
    geom_raster(aes(x=x, y=y, fill = bathymetry)) +
    #geom_line(aes(x=x, y = 40.8+0.4*x), linewidth = 1) +
    scale_fill_distiller(palette = "Spectral") +
    geom_sf(data = contour_golfe %>% st_transform(crs=crs(combined_rasters))) +
    geom_sf(data = migralion_eff)
}

# Exclude covars that have a correlation coeff > 0.8
corrplot::corrplot.mixed(combined_rasters2 %>% as_tibble() %>% cor(method = "pearson", use = "na.or.complete"))
# exclude: mean_automn_SST (cor w/ mean_winter_SST)
# slope, mean_VEL (cor w/ bathymetry)
# mean_autumn_SST cor w/ mean_winter_SST
layers_to_keep <- c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", #"mean_autumn_SST",
                    #"sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST",
                    "mean_SST", "sd_SST", #"slope",
                    "dist_to_shore", "bathymetry", "mean_CHL", 
                    #"mean_SAL", "mean_VEL", 
                    "sd_SAL",
                    "mean_SSH", "sd_SSH",  "sd_VEL")

combined_rasters3 <- combined_rasters2[[c(layers_to_keep)]]
corrplot::corrplot.mixed(combined_rasters3 %>% as_tibble() %>% cor(method = "pearson", use = "na.or.complete"))

# Add the log values
log_raster <- log(combined_rasters3)
log_raster$concavity <- combined_rasters3$concavity
log_raster$mean_SSH <- combined_rasters3$mean_SSH
log_raster$bathymetry <- log(-combined_rasters3$bathymetry)
names(log_raster) <- paste0("log_", names(log_raster))


final_raster <- c(log_raster,combined_rasters3) %>% scale()
plot(final_raster)

# save 
terra::writeRaster(final_raster, filename = "1.data/all_covariates.tif", overwrite=TRUE)
