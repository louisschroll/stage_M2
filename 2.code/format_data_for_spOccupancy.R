

get_data_for_spOccupancy <- function(data_list, grid, species, add_coords=F){
  data_list <- filter_all_dataset(data_list, species)
  sitesocc <- extract_sampled_cells(data_list, grid = grid) # store which cells of `grid` are sampled, for each dataset
  sites_id <- get_sampled_cells_id(sitesocc)
  PA_and_detcov <- get_detection_data(data_list, sitesocc)
  occurence_covs <- get_occ_cov(sitesocc$grid)
  
  data.int <- list(y = PA_and_detcov$y, 
                   occ.covs = occurence_covs, 
                   det.covs = PA_and_detcov$det.covs, 
                   sites = sites_id)
  if (add_coords){
    data.int$coords = get_coords(sitesocc$grid)
  }
  return(data.int)
}


# ---- filter one species ----
filter_all_dataset <- function(dataset_list, species){
  filtered_dataset_list <- map(dataset_list, ~ filter_one_dataset(.x, species)) 

  removed_dataset = names(dataset_list)[unlist(map(filtered_dataset_list, is.null))]
  if (length(removed_dataset) > 0){
    print(paste0(removed_dataset, 
                   " contain no observation for this species (",
                  species, ")"))
  }
  # remove NULL values with compact (when no observation are made in a dataset)
  return(compact(filtered_dataset_list))
}


filter_one_dataset <- function(dataset, species){
  subset_obs <- dataset$obs %>% filter(nom_fr == species)
  if (nrow(subset_obs) == 0){
    return(NULL)
  }
  effective_session <- subset_obs %>% pull(session) %>% unique()
  sub_eff <- dataset$eff %>% filter(session %in% effective_session)
  return(list(obs = subset_obs, eff = sub_eff))
}


# ---- find which sites are sampled form the grid ----
extract_sampled_cells <- function(data_list, grid = grid){
  # find the grid-cells that intersect with sampling effort
  sites <- list()
  nb_of_datasets <- length(data_list)
  
  for (i in 1:nb_of_datasets){
    sites[[i]] <- find_intersecting_cells(data_list[[i]]$eff, grid)
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

find_intersecting_cells <- function(sampling_points, grid){
  
  intersection <- st_intersection(sampling_points, grid) %>%
    as.data.frame() %>%
    pull(id) %>%
    unique() 

  return(intersection)
}



convert_site_ID <- function(sites, gridocc){
  id_df = gridocc %>% 
    as_tibble() %>% 
    select(id, ido)
  
  sitesocc <- map(sites, ~ id_df %>%
                    filter(id %in% .x) %>%
                    pull(ido))  
  return(sitesocc)
}


get_sampled_cells_id <- function(sitesocc){
  # get sampled cells id for each dataset
  sites_id <- sitesocc[-length(sitesocc)] # or sitesocc[!names(sitesocc) %in% "grid"]
}


# ---- retrieve or compute PA data, covariates and coords ----
get_detection_data <- function(data_list, sitesocc){
  # Give the correct format to the data (as defined in spOccupancy)
  y_list <- list()
  det_cov_list <- list()
  
  for (i in 1:length(data_list)){
    PA_data <- get_y_and_detcov(obs_data = data_list[[i]]$obs, 
                                eff_data = data_list[[i]]$eff, 
                                gridocc = sitesocc$grid, 
                                sitesocc_id = sitesocc[[i]])
    y_list[[i]] <- PA_data$y
    det_cov_list[[i]] <- PA_data$det.covs
  }
  return(list(y = y_list, det.covs = det_cov_list))
}


get_y_and_detcov <- function(obs_data, eff_data, gridocc, sitesocc_id){
  
  session_list = unique(obs_data$session) %>% sort(decreasing = F)
  print(session_list)
  
  occurence_df <- gridocc[sitesocc_id,] %>% 
    mutate(id_data = 1:length(sitesocc_id))
  
  intersect_obs <- st_intersection(obs_data, occurence_df) 
  intersect_eff <- st_intersection(eff_data, occurence_df) 

  occurence_list = c()
  effort_list = c()
  store_session = c()
  effort_matrix = matrix(NA, ncol = length(session_list), nrow = nrow(occurence_df))
  #ii = 1
  for (S in session_list){
    intersect_eff_session <- intersect_eff %>% filter(session == S) #st_intersection(eff_data_y, occurence_df) 
    df_transect_length <- tibble(id = intersect_eff_session$id_data, len = st_length(intersect_eff_session)) %>% 
      group_by(id) %>% 
      summarise(total_length = sum(len))
    
    occurence_df$effort <- 0
    occurence_df$effort[df_transect_length$id] <- df_transect_length$total_length

    ### fill det/ non-det
    id_cells_with_obs <- intersect_obs %>% filter(session == S) %>% pull(id_data) %>% unique()
    occurence_df <- occurence_df %>% 
      mutate(y = ifelse(effort > 0, 0, NA)) %>% # fill cells with effort w/ 0
      mutate(y = ifelse(row_number() %in% id_cells_with_obs, 1, y)) %>% # fill cells with obs w/ 1
      mutate(y = ifelse(effort <= 0, NA, y)) # remove potential obs in cells absent from effort data this session
    
    effort_list = c(effort_list, occurence_df$effort)
    # effort_matrix[,ii] <- occurence_df$effort
    # ii = ii + 1
    occurence_list = c(occurence_list, occurence_df$y)
    store_session = c(store_session, rep(S, nrow(occurence_df)))
  }
  
  y <- matrix(occurence_list, ncol = length(session_list))
  
  detcov <- matrix(effort_list, ncol = length(session_list))
  detcov[detcov == 0] <- NA
  detcov <- scale(detcov)
  # check dimensions 
  if( !(length(which(is.na(y))) == length(which(is.na(detcov))) ) ){
    print("Warnings: inconsistent dimensions in data")
  }
  data <- list()
  data$y <- y
  data$det.covs <- list()
  data$det.covs$transect_length <- detcov
  data$det.covs$session <- matrix(store_session, ncol = length(session_list))
  return(data)
}


get_coords <- function(grid){
  coord = st_coordinates(st_centroid(grid))
  return(coord)
}


get_occ_cov <- function(grid){
  occ_cov_df <- grid %>% 
    as_tibble() %>% 
    select(-grid_c, -geometry, -id, -ido)
  return(occ_cov_df)
}


