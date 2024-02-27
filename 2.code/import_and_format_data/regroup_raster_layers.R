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
setwd("/1.data")

# Import packages
library(terra)
library(tidyverse)

# Define variables used in the script
study_area <- ext(3, 6, 42, 43.7) # (xmin, xmax, ymin, ymax)


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
raster_sst <- combine_ncfile("copernicus_data/SST.nc", "copernicus_data/SST_int.nc") %>% 
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
            sd_value = sd(value, na.rm = TRUE)) %>% 
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
  terra::rast(type="xyz", crs="", digits=6, extent=NULL) 
 
all_SST <- c(app(raster_sst, mean),app(raster_sst, sd)) %>% 
  resample(raster_mean_SST)
names(all_SST) <- c("mean_SST", "sd_SST")
raster_mean_SST <- c(raster_mean_SST, all_SST)
plot(raster_mean_SST)
corrplot::corrplot.mixed(raster_mean_SST %>% as_tibble() %>% cor(method = "pearson", use = "na.or.complete"))

# ---- Import and crop static covariates data ----

raster_static_cov <- rast("static_covariates_raster.tif")
crs(raster_mean_SST) <- crs(raster_static_cov)

# Use resample to have same extent and resolution
raster_static_cov <- raster_static_cov %>% 
  terra::resample(raster_mean_SST) 

# ---- Import and average the other dynamic covariates ----
# Salinity (SAL), euphotic zone depth (EZD) and Net primary productivity are
# excluded because their correlation whith CHL is too high
# Sea surface heigth (SSH) is too correlated with dist_to_shore and bathymetry
# And sd(chl) to correlated with mean(chl)

CHL <- combine_ncfile("copernicus_data/CHL.nc", "copernicus_data/CHL_int.nc") 
#SAL <- combine_ncfile("copernicus_data/SAL.nc", "copernicus_data/SAL_int.nc") 
#SSH <- combine_ncfile("copernicus_data/SSH.nc", "copernicus_data/SSH_int.nc") 

raster_dynamic_cov <- c(app(CHL, mean)#, app(CHL, sd),
  #app(SAL, mean), app(SAL, sd),
  #app(SSH, mean), app(SSH, sd)
  ) %>% 
  terra::crop(raster_mean_SST) %>% 
  terra::extend(raster_mean_SST) 

VEL <- combine_ncfile("copernicus_data/velocity.nc", "copernicus_data/velocity_int.nc") %>% 
  terra::resample(raster_mean_SST) 

# EZD <- rast("copernicus_data/euphotic_depth.nc")%>% 
#   terra::resample(raster_mean_SST) 
# NPP <- rast("copernicus_data/NPP.nc") %>% 
#   terra::resample(raster_mean_SST) 
raster_dynamic_cov <- c(raster_dynamic_cov, 
                        app(VEL, mean), app(VEL, sd)
                        #app(EZD, mean), app(EZD, sd),
                        #app(NPP, mean), app(NPP, sd)
                        )

names(raster_dynamic_cov) <- c("mean_CHL", #"sd_CHL",
                               #"mean_SAL", "sd_SAL", 
                               #"mean_SSH", "sd_SSH", 
                               "mean_VEL", "sd_VEL"
                               #"mean_EZD", "sd_EZD", 
                               #"mean_NPP", "sd_NPP"
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

# Add the log values
log_raster <- log(combined_rasters)
log_raster$concavity <- combined_rasters$concavity
log_raster$mean_VEL <- combined_rasters$mean_VEL
log_raster$bathymetry <- log(-combined_rasters$bathymetry)
names(log_raster) <- paste0("log_", names(log_raster))


final_raster <- c(log_raster,combined_rasters) %>% scale()
plot(final_raster)
#combined_rasters[combined_rasters$bathymetry<(-1200)] <- NA
plot(combined_rasters)

corrplot::corrplot.mixed(combined_rasters %>% as_tibble() %>% cor(method = "pearson", use = "na.or.complete"))

terra::writeRaster(final_raster, filename = "all_covariates.tif", overwrite=TRUE)
