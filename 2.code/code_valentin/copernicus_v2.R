# use CopernicusMarine package
library(CopernicusMarine)
library(ggplot2)
library(stars)
load("~/stage_M2/1.data/gdlmap.rdata")
# destination <- tempfile("copernicus", fileext = ".nc")

info <- cms_products_list(freeText = "GLOBAL_MULTIYEAR_PHY_001_030")
cms_product_details(product ="GLOBAL_MULTIYEAR_PHY_001_030",
                      layer    = "cmems_mod_glo_phy-thetao_myint_0.083deg_P1M-m")
 

# cms_download_subset(
#   destination   = "Data/test.nc",
#   product       = "GLOBAL_MULTIYEAR_PHY_001_030",
#   layer         = "cmems_mod_glo_phy_myint_0.083deg_P1M-m",
#   variable      = "thetao",
#   region        = c(2, 42, 6, 44),
#   timerange     = c("2022-01-01", "2022-02-01"),
#   verticalrange = c(0, -1),
#   overwrite = T
# )

# options(CopernicusMarine_uid= "")
# options(CopernicusMarine_pwd="")

# Non interim period
cms_download_subset(
  destination   = "~/stage_M2/1.data/copernicus_data/SST.nc",
  product       = "MEDSEA_MULTIYEAR_PHY_006_004",
  layer         = "med-cmcc-tem-rean-m",
  variable      = "thetao",
  region        = c(2, 42, 6, 44),
  timerange     = c("2000-01-01", "2023-06-01"),
  verticalrange = c(0, -1),
  overwrite = T
)

# interim period 
cms_download_subset(
  destination   = "~/stage_M2/1.data/copernicus_data/SST_int.nc",
  product       = "MEDSEA_MULTIYEAR_PHY_006_004",
  layer         = "med-cmcc-tem-int-m",
  variable      = "thetao",
  region        = c(2, 42, 6, 44),
  timerange     = c("2021-07-01", "2023-12-01"),
  verticalrange = c(0, -1),
  overwrite = F
)

#> Preparing job...
#> Waiting for job to finish...
#> Downloading file...
#> Done

mydata <- stars::read_ncdf("~/stage_M2/1.data/copernicus_data/SST_int.nc",
                           make_time = TRUE,
                           var = c("thetao"))

class(mydata)
dim(mydata)
#> vo, uo,

# plot 
ggplot() + 
  geom_stars(data = mydata[ , , , , 1], aes(color = thetao), sf = T) +
  geom_sf(data = gdlmap %>% st_transform(crs = st_crs(mydata)))

ggplot() + 
  geom_stars(data = mydata[ , , , , 29], aes(color = thetao), sf = T) +
  geom_sf(data = gdlmap %>% st_transform(crs = st_crs(mydata)))

