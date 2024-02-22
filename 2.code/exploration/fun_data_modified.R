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


