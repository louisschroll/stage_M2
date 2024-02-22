
library(tidyverse)
library(sf)


format_pelmed_data <- function(pelmed_obs=pelmed, pelmed_eff=pelmed_eff, sitesocc=sitesocc){
  # Pelmed
  pelmed_obs <- pelmed_obs %>% 
    mutate(year = lubridate::year(date)) 
  
  pelmed_formated = format_data(obs_data = pelmed_obs,
                   eff_data = pelmed_eff,
                   gridocc = sitesocc$grid,
                   sitesocc_id = sitesocc$pelmed)
  
  return(pelmed_formated)
}


format_samm_data <- function(samm_obs = samm, samm_eff = effsamm, sitesocc=sitesocc){
  # samm
  samm_obs <- samm_obs %>% 
    mutate(year = lubridate::year(date)) %>% 
    st_transform(crs = st_crs(sitesocc$grid))
  
  samm_eff <- samm_eff %>% 
    mutate(year = lubridate::year(date)) %>% 
    st_transform(crs = st_crs(sitesocc$grid))
  
  samm_formated = format_data(obs_data = samm_obs,
                   eff_data = samm_eff,
                   gridocc = sitesocc$grid,
                   sitesocc_id = sitesocc$samm)
  
  return(samm_formated)
}


format_pnm_data <- function(pnm_obs = pnm, pnm_eff = transect, sitesocc=sitesocc){
  # PNM
  pnm_obs <- pnm_obs %>% 
    mutate(year = lubridate::year(dmy(date)),
           month = lubridate::month(dmy(date))) %>% 
    st_transform(crs = st_crs(sitesocc$grid))
  
   pnm_formated = format_data(obs_data = pnm_obs,
                   eff_data = pnm_eff,
                   gridocc = sitesocc$grid,
                   sitesocc_id = sitesocc$pnm)
   
   return(pnm_formated)
}

format_migralion_data <- function(migralion_obs = migralion, 
                                  migralion_eff = postnup2022_eff, sitesocc=sitesocc){
  # Migralion
  migralion_obs <- migralion_obs %>% 
    mutate(year = lubridate::year(Date_UTC),
           month = lubridate::month(Date_UTC)) %>% 
    mutate(year = ifelse(month<7, 1, 2))
  
  migralion_eff <- st_transform(migralion_eff, crs = st_crs(sitesocc$grid))
  
  migralion_formated = format_data(obs_data = migralion_obs,
                   eff_data = migralion_eff,
                   gridocc = sitesocc$grid,
                   sitesocc_id = sitesocc$migralion,
                   is_migralion = T)
  
  return(migralion_formated)
}


format_data <- function(obs_data, eff_data, gridocc, sitesocc_id, is_migralion=F){
  formated_data = get_y_and_detcov(obs_data, eff_data, gridocc, sitesocc_id, is_migralion)
  formated_data$occ.covs = get_occ_cov(gridocc, sitesocc_id = sitesocc_id)
  formated_data$coords = get_coord(gridocc, sitesocc_id=sitesocc_id)
  return(formated_data)
}


get_y_and_detcov <- function(obs_data, eff_data, gridocc, sitesocc_id, is_migralion=F){

  year_list = unique(obs_data$year) %>% sort(decreasing = F)
  print(year_list)
  
  occurence_df <- gridocc[sitesocc_id,] %>% 
    mutate(id_data = 1:length(sitesocc_id))
  
  intersect_obs <- st_intersection(obs_data, occurence_df) 
  intersect_eff <- st_intersection(eff_data, occurence_df) 
  intersect_eff$length <- st_length(intersect_eff)
  
  occurence_list = c()
  effort_list = c()
  for (YEAR in year_list){
    ### fill seff
    occurence_df$effort <- 0
    for(i in 1:nrow(occurence_df)){
      id <- which(intersect_eff$ido == occurence_df$ido[i])
      if (is_migralion){
        occurence_df$effort[i] <- sum(intersect_eff$length[id])
      }
      else {
        occurence_df$effort[i] <- sum(intersect_eff$length[id][year(intersect_eff$date)[id] == YEAR])
      }
    }
    
    ### fill det/ non-det
    occurence_df$y <- NA
    #### fill w/ 0
    occurence_df$y[occurence_df$effort > 0] <- 0
    #### fill w/ 1
    occurence_df$y[unique(intersect_obs$id_data[intersect_obs$year==YEAR])] <- 1
    
    effort_list = c(effort_list, occurence_df$effort)
    occurence_list = c(occurence_list, occurence_df$y)
  }
  
  y <- matrix(occurence_list, ncol = length(year_list))
  
  detcov <- matrix(effort_list, ncol = length(year_list))
  detcov[detcov == 0] <- NA
  detcov <- scale(detcov)
  
  # check dimensions 
  diff.na <-  length(which(is.na(y)==T)) - length(which(is.na(detcov)==T))
  
  if( !(length(which(is.na(y))) == length(which(is.na(detcov))) ) ){
    print("Warnings: inconsistent dimensions in pelmed data")
  }
  if( diff.na < 0){
    id <- setdiff(which(is.na(detcov)), which(is.na(y)))
    detcov[id] <-  0
  }
  if( diff.na > 0){
    id <- setdiff(which(is.na(y)) ,which(is.na(detcov)))
    y[id] <-  0
  }
  if( length(which(is.na(y))) == length(which(is.na(detcov)) ) ){
    print("Message: dimensions corrected")
  }
  
  # return ...
  data <- list()
  data$y <- y
  data$det.covs <- list()
  data$det.covs$transect_length <- detcov
  return(data)
}


get_coord <- function(gridocc, sitesocc_id){
  coord = st_coordinates(st_centroid(gridocc[sitesocc_id,]))
  return(coord)
}


get_occ_cov <- function(gridocc, sitesocc_id){
  occ_cov_df <- gridocc[sitesocc_id,] %>% 
    as.data.frame() %>% 
    select(bathymetry, dist_to_shore, slope)
  return(occ_cov_df)
}




# find which sites are sampled form the grid 
extract_sampled_cells <- function(do.pel = T, do.sam = T, do.pnm = T, do.mig = T, grid = grid){
  # find the grid-cells that intersect with sampling effort
  sites <- list()
  
  # pelmed
  if(do.pel){
    sites$pelmed <- find_intersecting_cells(pelmed_eff, grid, id_col = "id")
  }
  
  # pnmgdl
  if(do.pnm){
    sites$pnm <- find_intersecting_cells(transect, grid, id_col = "id.1")
  }
  
  # migralion
  if(do.mig){
    sites$migralion <- find_intersecting_cells(postnup2022_eff, grid, id_col = "id")
  }
  
  # SAM
  if(do.sam){
    # remove surveys from 2019 as no detection has been recorded
    effsamm1 <- effsamm %>% 
      mutate(year = lubridate::year(date)) %>% 
      st_transform(crs = st_crs(grid)) %>% 
      filter(year %in% c("2011","2012"))
    
    sites$samm <- find_intersecting_cells(effsamm1, grid, id_col = "id")
    sites$samm <- sites$samm[-c(1232:1235)]
  }
  
  IDs <- unique(unlist(sites))
  # gridocc stores the sites sampled at least once by a dataset
  gridocc <- grid[IDs,] %>% 
    mutate(ido = 1:length(IDs))
  
  # convert site IDs to have them between 1 and length(IDs)
  sitesocc <- convert_site_ID(sites, gridocc)
  
  # export the restricted grid
  sitesocc$grid <- gridocc
  
  return(sitesocc) 
  # sitesocc is the 'sites' list we needed 
}


find_intersecting_cells <- function(sampling_points, grid, id_col){
  
  intersection <- st_intersection(sampling_points, grid) %>%
    as.data.frame() %>%
    select(all_of(id_col)) %>%
    unique() 
  
  return(intersection[[1]])
}


convert_site_ID <- function(sites, gridocc){
  sitesocc <- list()
  
  id_df = gridocc %>% 
    as_tibble() %>% 
    select(id, ido) %>% 
    pivot_wider(names_from = id, values_from = ido)
  
  for (i in 1:length(sites)){
    ido_convert = id_df[as.character(sites[[i]])] %>% 
      pivot_longer(cols = everything(), names_to = "id", values_to = "ido") %>% 
      select(ido) 
    sitesocc[[i]] = ido_convert$ido
  }
  
  names(sitesocc) = names(sites)
  
  return(sitesocc)
}


