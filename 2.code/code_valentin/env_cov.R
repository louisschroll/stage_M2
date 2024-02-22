# extract environmental covariates from a grid and a set of dates

# load packages
library(tidyverse)
library(sf)

# load the grid and the boundaries
load("gdlhex.rdata")
load("gdlmap.rdata")

if(FALSE){
  ggplot() +
    geom_sf(data = sea.gdl2, lwd = 0.1) 
  ggplot()+ 
    geom_sf(data = gdlmap) 
}

# ---- dist to coast ----
st_crs(sea.gdl2$x) == st_crs(gdlmap$geometry)

dist_to_coast <- unlist(map2(sea.gdl2$x, st_union(gdlmap), st_distance))

if(FALSE){
  ggplot() +
    geom_sf(data = sea.gdl2, aes(fill = dist_to_coast), lwd = 0.1) + 
    scale_fill_viridis_c(direction = -1 ) +
    geom_sf(data = gdlmap) 
}

# ---- dist to colony (the closest) based on GISOM field count ----

species <- "Sterne caugek"

load("countLL.rdata")

colony <- countLL %>% 
  filter(Species == "Sterne caugek") %>% 
  group_by(Year) %>% 
  nest() 

head(colony$data[1])
size_col <- dist_to_col <- matrix(NA, nrow = nrow(sea.gdl2), ncol = nrow(colony))
for(y in 1:nrow(colony)){
  dist_to_col[,y] <- unlist(map2(sea.gdl2$x, st_union(colony$data[[y]]), st_distance))
  # size_col[,y] <- unlist(map2(sea.gdl2$x, st_union(colony$data[[y]]), st_distance))
}

colnames(dist_to_col) <- colony$Year

dist_col <- sea.gdl2 %>% 
  bind_cols(dist_to_col)
save(dist_col, file = here("Data/Fond de carte/dist_to_col.rdata"))
if(FALSE){
  ggplot() +
    geom_sf(data = sea.gdl2, aes(fill = dist_to_col[,1]), lwd = 0.1) + 
    scale_fill_viridis_c(direction = -1 , option = "B") +
    geom_sf(data = gdlmap) 
}

# ---- Extract dynamic data from Copernicus ----
library(ncdf4)
library(lubridate)


## Set your Copernicus Marine credentials
USERNAME <- readline("Enter your username: ")                  
PASSWORD <- readline("Enter your password: ")

copernicus_extract <- function(datasetID = datasetID, variable = "analysed_sst",
                               depth = logical(), 
                               time_interval = c("2020-01-01", "2021-01-01")){
  
  ## Open connection
  url <- paste("https://", USERNAME, ":", PASSWORD,"@my.cmems-du.eu/thredds/dodsC/", datasetID, sep ="") 
  
  
  ## Open file and check metadata
  ds <- nc_open(url)
  print(ds)
  
  # Info about the SST variable
  ncatt_get(ds,variable)
  
  # check columns names 
  if("latitude" %in% names(ds$dim)){
    # Longitude
    lon <- ncvar_get(ds, "longitude")
    nlon <- dim(lon)
    
    # Latitude
    lat <- ncvar_get(ds, "latitude")
    nlat <- dim(lat)
  } else{
    # Longitude
    lon <- ncvar_get(ds, "lon")
    nlon <- dim(lon)
    
    # Latitude
    lat <- ncvar_get(ds, "lat")
    nlat <- dim(lat)
  }
  
  # Check dimensions
  print(c(nlon,nlat))
  
  # Time
  time<-ncvar_get(ds,"time")
  nt <- dim(time)
  t_units <- ncatt_get(ds, "time", "units")
  t_units$value
  
  # convert time -- split the time units string into fields
  t_ustr <- strsplit(t_units$value, " ")
  t_dstr <- strsplit(unlist(t_ustr)[3], "-")
  
  if(word(t_units$value) == "minutes"){
    date <- lubridate::ymd(t_dstr) + lubridate::dminutes(time)
  } 
  if(word(t_units$value) == "seconds"){
    date <- lubridate::ymd(t_dstr) + lubridate::dseconds(time)
    range(date)
  }
  
  ## Define the parameters and ranges for subset
  #Bounding box
  bbox <- st_bbox(st_transform(sea.gdl2, crs = 4326))
  x <- c(bbox["xmin"], bbox["xmax"])                # longitude
  y <- c(bbox["ymin"], bbox["ymax"])                # latitude
  t <- time_interval   # time
  
  
  # Function to get the indices from the ranges
  btw <- function(data, num){
    c(min(which(num<=data)), max(which(num>=data)))
  }
  # Starting indices
  lon_indices <- btw(data = lon, num = x)
  lat_indices <- btw(data = lat, num = y)
  time_indices <- btw(data = date, num = t)
  
  # Count number of indices to extract along each dimension
  lon_range <- lon_indices[-1] - lon_indices[1]+1
  lat_range <- lat_indices[-1] - lat_indices[1]+1
  time_range <- time_indices[-1] - time_indices[1]+1
  
  
  if(depth == T){
    # Depth
    # z <- c(1)                       # depth
    # depth<-ncvar_get(ds,"depth")
    # dim(depth)
    # range(depth)
    depth_indices <- 1 #btw(data = depth, num = z)
    depth_range <- 1 # depth_indices[-1] - depth_indices[1]+1
    
    offset <- c(lon_indices[1], lat_indices[1], depth_indices[1], time_indices[1])    #lon,lat,depth,time
    count <- c(lon_range, lat_range, depth_range, time_range)
  }else{
    # Start and Count vectors
    offset <- c(lon_indices[1], lat_indices[1], time_indices[1])    #lon,lat,depth,time
    count <- c(lon_range, lat_range,  time_range)
  }
  
  # Get subsetted variable   
  sstx <- ncvar_get(ds,variable, start = offset, count = count)
  dim(sstx)
  # lon / lat / (depth) / time
  
  xlon <- lon[lon_indices[1]:lon_indices[2]]
  xlat <- lat[lat_indices[1]:lat_indices[2]]
  xtime <- date[time_indices[1]:time_indices[2]]
  
  range(xlon)
  range(xlat)
  # make the tibble
  sst_tib <- matrix(NA, nrow = length(xlon) * length(xlat), ncol = (dim(sstx)[3]+2))
  colnames(sst_tib) <- c("lon","lat",as.character(xtime))
  id <- 0
  
  for(i in 1:length(xlon)){
    for(j in 1:length(xlat)){
      id <- id+1
      sst_tib[id,1] <- xlon[i] # longitude
      sst_tib[id,2] <- xlat[j] # latitude
      sst_tib[id,3:(dim(sstx)[3]+2)] <- sstx[i,j,]
    }
  }
  
  apply(sst_tib[,1:2],2,range)
  
  which(is.na(sst_tib[,1]))
  # spatialize
  sstt <- sst_tib %>% 
    as.data.frame() %>%
    mutate(var_mean = rowMeans(sst_tib[,5:10], na.rm=T)) %>% 
    st_as_sf(coords = c("lon","lat"), crs = 4326) %>% 
    st_transform(crs = st_crs(gdlmap))
  
  return(sstt)
  
} # end function

time_check <- function(url = url){
  
  ds <- nc_open(url)
  
  # Time
  time<-ncvar_get(ds,"time")
  nt <- dim(time)
  t_units <- ncatt_get(ds, "time", "units")
  t_units$value
  
  # convert time -- split the time units string into fields
  t_ustr <- strsplit(t_units$value, " ")
  t_dstr <- strsplit(unlist(t_ustr)[3], "-")
  
  if(word(t_units$value) == "minutes"){
    date <- lubridate::ymd(t_dstr) + lubridate::dminutes(time)
  } 
  if(word(t_units$value) == "seconds"){
    date <- lubridate::ymd(t_dstr) + lubridate::dseconds(time)
    range(date)
  }
  return(range(date))
}
# ---- SST ----

# datasetID
datasetID <- "cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021"

variable = "analysed_sst"

## Open connection
url <- paste("https://",USERNAME, ":", PASSWORD,"@my.cmems-du.eu/thredds/dodsC/",datasetID, sep ="") 

time_check(url)

sstt <- copernicus_extract(datasetID = datasetID, variable = "analysed_sst",
                           time_interval = c("2019-01-01", "2022-10-31"), depth = F)

st_crs(sstt$geometry)
names(sstt)
# 
if(FALSE){
  ggplot() + 
    geom_sf(data = sstt %>% 
              filter(var_mean > 0), aes(color = (var_mean-273))) +
    geom_sf(data = gdlmap) +
    scale_color_viridis_c()
  
}

# ---- SSS----
datasetID <- "med-cmcc-sal-rean-m" #"cmems_mod_med_phy-sal_anfc_4.2km_P1D-m"

## Open connection
url <- paste("https://",USERNAME, ":", PASSWORD,"@my.cmems-du.eu/thredds/dodsC/",datasetID, sep ="") 

## Open file and check metadata
ds <- nc_open(url)
print(ds)

variable = "so"

time_check(url)

sssm <- copernicus_extract(datasetID = datasetID, variable = "so", depth = T,
                           time_interval = c("2019-01-01", "2021-06-30"))

names(sss)
if(FALSE){
  p_sal <-   ggplot() + 
    geom_sf(data = sss , aes(color = (var_mean))) +
    geom_sf(data = gdlmap) + 
    scale_color_viridis_c()
  p_sal
}

# ---- ChlA ----

datasetID <-"cmems_obs-oc_med_bgc-plankton_my_l3-multi-1km_P1D"  #"med-ogs-bio-rean-d" # 

## Open connection
url <- paste("https://",USERNAME, ":", PASSWORD,"@my.cmems-du.eu/thredds/dodsC/",datasetID, sep ="") 

## Open file and check metadata
ds <- nc_open(url)
print(ds)

variable = "CHL" # "nppv" # 

time_check(url)

chl <- copernicus_extract(datasetID = datasetID, variable = variable,
                          depth = F,time_interval = c("2019-01-01", "2022-10-31"))

names(chl)
if(FALSE){
  p_chl <-   ggplot() + 
    geom_sf(data = chl %>% 
              filter(var_mean > 0), aes(color = (var_mean))) +
    geom_sf(data = gdlmap) + 
    scale_color_viridis_c()
  
  p_chl #| p_sal
}

# ---- summary and export cov----


# length(unique(month(ymd(names(sstt))))) * length(unique(year(ymd(names(sstt)))))

dealthecov <- function(sstt = sstt){
  
  gridcov <- sea.gdl2
  # sst 
  sstt2 <- sstt %>% 
    filter(is.na(var_mean)== F)
  
  
  head(sstt[,1:10])
  
  tmp <- tibble()
  for(i in 1:nrow(sea.gdl2)){
    ssti <- st_intersection(sea.gdl2[i,],sstt2)
    if(nrow(ssti)> 0){
      sstii <- tibble(date = names(sstt2),sst = c(as.numeric(ssti[1,3:ncol(ssti)])))
      
      tmpi <- sstii %>% 
        mutate(year = year(ymd(date)),
               month = month(ymd(date))) %>% 
        group_by(year, month) %>% 
        summarise(sstmean= mean(sst)) %>% 
        ungroup() %>% 
        pivot_wider(names_from = c("year","month"), values_from = "sstmean",names_sep = ".") %>% 
        mutate(id = i) %>% 
        as.matrix()
      
      tmp <- rbind(tmp, tmpi)
    }
  }
  
  
  gridcov <- sea.gdl2 %>% 
    mutate(id = 1:nrow(sea.gdl2)) %>% 
    left_join(tmp, by = "id")
  
  #  deal with NAs # pfiou on peut faire mieux mais je suis fatigué ce soir
  nas <- which(is.na(gridcov$`2019.1`))
  
  if(length(nas)>0){
    gridnas <- gridcov %>% 
      filter(is.na(gridcov$`2019.1`))
    gridval <- gridcov %>% 
      filter(!is.na(gridcov$`2019.1`))
    
    #sf::st_drop_geometry(data)
    for(i in 1:nrow(gridnas)){
      id <- gridnas$id[i]
      ptna <- st_centroid(gridnas[i,])
      t <- which(st_distance(ptna,gridval) == min(st_distance(ptna,gridval)))
      gridcov[id,3:(ncol(gridcov)-1)] <- gridval[t,3:(ncol(gridval)-1)]
    }
  }
  
  return(gridcov)
}

## sst

sstcov <- dealthecov(sstt)

## sssm
head(sssm)
dates <- colnames(sssm)
dates2 <- dates %>% 
  ymd_hms() %>% 
  floor_date(unit = "day") %>% 
  as.character()

dates2

colnames(sssm)[1:(ncol(sssm)-2)] <- dates2[1:(ncol(sssm)-2)]
head(sssm)
subsssm <- sssm[2222:2322,] # test
ssscov <- dealthecov(sst= sssm)

ggplot( ) + geom_sf(data = ssscov, aes(fill = `2019.1`), lwd = 0.1)

## chla
head(chl)

chlcov <- dealthecov(sst= chl)

ggplot( ) + geom_sf(data = ssscov, aes(fill = `2019.1`), lwd = 0.1)

save(ssscov, sstcov, chlcov, file = "chl_sss_sst_monthly2019-22.rdata")

library(CopernicusMarine)
copernicus_product_details("GLOBAL_MULTIYEAR_PHY_001_030")

