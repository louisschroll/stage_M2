# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/import_and_format_data/cov_data_in_polygon_grid.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-21
#
# Script description:
# Transform the raster into a polygon grid
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load package
library(stars)
library(sf)
library(tidyverse)
raster_cov = read_stars("1.data/all_covariates.tif")

grid_c <- st_make_grid(raster_cov, 
                       cellsize = 0.04,
                       what = "polygons",
                       square = FALSE)

polygon_values <- split(raster_cov, "band") %>% 
  st_extract(at = st_centroid(grid_c), mean()) %>% 
  split("band")

covariates_data <- st_sf(grid_c) %>%
  mutate(polygon_values$band) %>% 
  filter(!if_any(mean_winter_SST:sd_VEL, is.na))

plot(covariates_data[10:11])

save(covariates_data, file = "1.data/covariates_data.rdata")

