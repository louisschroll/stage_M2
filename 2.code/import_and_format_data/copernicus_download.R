

# ---- Download environmental data from Copernicus ----
# https://data.marine.copernicus.eu/products?facets=tempResolutions%7EMonthly--areas%7EMediterranean+Sea

library(CopernicusMarine) 
library(terra)
# options(CopernicusMarine_uid= "")
# options(CopernicusMarine_pwd="")

x_min = 0
x_max = 6
y_min = 41
y_max = 45
study_area = c(x_min, y_min, x_max, y_max)

path_save <-  "~/stage_M2/1.data/copernicus_data/"

t_min <- "2011-01-01"
t_max <- "2022-08-31"
t_min_interim <- "2022-08-01"
t_max_interim <- "2023-12-31"

# Dataset: Mediterranean Sea Physics Reanalysis
product_name <- "MEDSEA_MULTIYEAR_PHY_006_004"
## Sea Surface Temperature (SST)
  # Non interim period
cms_download_subset(
  destination   = paste0(path_save, "SST.nc"),
  product       = product_name,
  layer         = "med-cmcc-tem-rean-m",
  variable      = "thetao",
  region        = study_area,
  timerange     = c(t_min, t_max),
  verticalrange = c(0, -1),
  overwrite = T
)

  # Interim period
cms_download_subset(
  destination   = paste0(path_save, "SST_int.nc"),
  product       = product_name,
  layer         = "med-cmcc-tem-int-m",
  variable      = "thetao",
  region        = study_area,
  timerange     = c(t_min_interim, t_max_interim),
  verticalrange = c(0, -1),
  overwrite = T
)

## Salinity (Sal)
# Non interim period
cms_download_subset(
  destination   = paste0(path_save, "SAL.nc"),
  product       = product_name,
  layer         = "med-cmcc-sal-rean-m",
  variable      = "so",
  region        = study_area,
  timerange     = c(t_min, t_max),
  verticalrange = c(0, -1),
  overwrite = T
)

# Interim period
cms_download_subset(
  destination   = paste0(path_save, "SAL_int.nc"),
  product       = product_name,
  layer         = "med-cmcc-sal-int-m",
  variable      = "so",
  region        = study_area,
  timerange     = c(t_min, t_max),
  verticalrange = c(0, -1),
  overwrite = T
)


## Sea Water velocity
# Non interim period
cms_download_subset(
  destination   = paste0(path_save, "velocity.nc"),
  product       = product_name,
  layer         = "med-cmcc-cur-rean-m",
  variable      = "sea_water_velocity",
  region        = study_area,
  timerange     = c(t_min, t_max),
  verticalrange = c(0, -1),
  overwrite = T
)

# Interim period
cms_download_subset(
  destination   = paste0(path_save, "SST_int.nc"),
  product       = product_name,
  layer         = "med-cmcc-tem-int-m",
  variable      = "thetao",
  region        = study_area,
  timerange     = c(t_min, t_max),
  verticalrange = c(0, -1),
  overwrite = T
)


## Sea Surface Height (SSH)
cms_download_subset(
  destination   = paste0(path_save, "SSH.nc"),
  product       = product_name,
  layer         = "med-cmcc-ssh-rean-m",
  variable      = "zos",
  region        = study_area,
  timerange     = c(t_min, t_max),
  overwrite = T
)

cms_download_subset(
  destination   = paste0(path_save, "SSH_int.nc"),
  product       = product_name,
  layer         = "med-cmcc-ssh-int-m",
  variable      = "zos",
  region        = study_area,
  timerange     = c(t_min_interim, t_max_interim),
  overwrite = T
)

# Dataset: Mediterranean Sea Biogechemistry Reanalysis
product_name <- "MEDSEA_MULTIYEAR_BGC_006_008"
t_max <- "2021-06-30"
t_min_interim <- "2021-07-01"
cms_download_subset(
  destination   = paste0(path_save, "CHL.nc"),
  product       = product_name,
  layer         = "med-ogs-pft-rean-m",
  variable      = "chl",
  region        = study_area,
  timerange     = c(t_min, t_max),
  verticalrange = c(0, -1),
  overwrite = T
)

cms_download_subset(
  destination   = paste0(path_save, "CHL_int.nc"),
  product       = product_name,
  layer         = "cmems_mod_med_bgc-pft_myint_4.2km_P1M-m",
  variable      = "chl",
  region        = study_area,
  timerange     = c(t_min_interim, t_max_interim),
  verticalrange = c(0, -1),
  overwrite = T
)


# Dataset:  Global Ocean Low and Mid Trophic Levels Biomass Content Hindcast Product
## euphotic zone depth
cms_download_subset(
  destination   = paste0(path_save, "euphotic_depth.nc"),
  product       = "GLOBAL_MULTIYEAR_BGC_001_033",
  layer         = "cmems_mod_glo_bgc_my_0.083deg-lmtl_PT1D-i",
  variable      = "zeu",
  region        = study_area,
  timerange     = c(t_min, "2022-12-31"),
  overwrite = T
)

## NPP: Net primary productivity
cms_download_subset(
  destination   = paste0(path_save, "NPP.nc"),
  product       = "GLOBAL_MULTIYEAR_BGC_001_033",
  layer         = "cmems_mod_glo_bgc_my_0.083deg-lmtl_PT1D-i",
  variable      = "npp",
  region        = study_area,
  timerange     = c(t_min, "2022-12-31"),
  overwrite = T
)



# ---- Download data from Global Fishing Watch ----
# 1st create an account and get a token on GFW website
# https://globalfishingwatch.org/
# remotes::install_github("GlobalFishingWatch/gfwr")
# https://github.com/GlobalFishingWatch/gfwr

library(gfwr)

key <- gfw_auth()

get_region_id(region_name = "France", region_source = "eez", key = key)
get_region_id(region_name = "Spain", region_source = "eez", key = key)

start_year = 2013 # no data before
end_year = 2023
ais_tibble <- tibble()
for (year in start_year:end_year){
  date_range <- paste0(year,"-01-01,",year,"-12-31")
  # Data in French EEZ
  ais_french_eez <- get_raster(
    spatial_resolution = "high",
    temporal_resolution = "monthly",
    group_by = "flag",
    date_range = date_range,
    region = 5677,
    region_source = "eez",
    key = key
  ) 
  
  # Data in Spanish EEZ
  ais_spain_eez <- get_raster(
    spatial_resolution = "high",
    temporal_resolution = "monthly",
    group_by = "flag",
    date_range = date_range,
    region = 5693,
    region_source = "eez",
    key = key
  ) 
  
  # Bind data & filter with coordinates
  new_data <- bind_rows(ais_french_eez, ais_spain_eez) %>% 
    janitor::clean_names() %>% 
    filter(lon >= 2 & lon <= 6 &
             lat >= 42 & lat <= 44)
  ais_tibble <- bind_rows(ais_tibble, new_data)
  
  print(paste("data for year", year, "collected"))
}

ais_tibble <- ais_tibble %>% mutate(year =substr(ais_tibble$time_range, 1, 4))

if(FALSE){
  ais_tibble %>% 
    filter(apparent_fishing_hours < 50) %>% 
    mutate(apparent_fishing_hours = log(apparent_fishing_hours)) %>% 
    ggplot(aes(x = lon, y = lat, fill = apparent_fishing_hours)) +
      geom_tile() +
      facet_wrap(~year) +
      labs(title = "Apparent Fishing Hours Map",
         x = "Longitude",
         y = "Latitude",
         fill = "Fishing Hours") +
      theme_minimal()
  
}

save(ais_tibble, file = "~/stage_M2/1.data/fishing_effort.rdata")



# ---- Download data from MARSPEC ----
# Download compressed file here: 
# https://www.esapubs.org/archive/ecol/E094/086/#data

# Open .7z files with archive
library(archive)
archive_extract(
  "~/biogeo01_07_30s.7z",
  dir = "MARSPEC_data",
  files = NULL
)

archive_extract(
  "~/bathymetry_30s.7z",
  dir = "MARSPEC_data",
  files = NULL
)

library(terra)
library(sf)

## Zoom in the map
ext <- ext(2, 42, 6, 44)

# retrieve bathymetry in m
depth_raster <- rast("~/MARSPEC_data/bathymetry_30s/bathy_30s") %>% 
  raster::crop(ext) %>% 
  as.numeric()

plot(depth_raster)

# retrieve distance to shore in km (biogeo_05), bathymetric slope in degrees (biogeo_06), concavity in degrees (biogeo_07)
file_nb <- c(5, 6, 7)
for (i in file_nb){
  depth_raster <- rast(paste0("~/MARSPEC_data/biogeo01_07_30s/biogeo0",i,"_30s")) %>% 
    crop(ext) %>% 
    as.numeric() %>% 
    c(depth_raster)
}


names(depth_raster) <- c("concavity", "slope", "dist_to_shore", "bathymetry")
varnames(depth_raster) <- c("concavity", "slope", "dist_to_shore", "bathymetry")

plot(depth_raster)

terra::writeRaster(depth_raster, filename = "~/1.data/static_covariates_raster.tif")

