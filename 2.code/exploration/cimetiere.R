get_y_and_detcov <- function(obs_data, eff_data, gridocc, sitesocc_id){
  
  session_list = unique(obs_data$session) %>% sort(decreasing = F)
  print(session_list)
  
  occurence_df <- gridocc[sitesocc_id,] %>% 
    mutate(id_data = 1:length(sitesocc_id))
  
  intersect_obs <- st_intersection(obs_data, occurence_df) 
  intersect_eff <- st_intersection(eff_data, occurence_df) 
  # intersect_eff$length <- st_length(intersect_eff)
  
  occurence_list = c()
  effort_list = c()
  transect_name_list = c()
  store_session = c()
  effort_matrix = matrix(NA, ncol = length(session_list), nrow = nrow(occurence_df))
  #ii = 1
  for (S in session_list){
    intersect_eff_session <- intersect_eff %>% filter(session == S) #st_intersection(eff_data_y, occurence_df) 
    df_transect_length <- tibble(id = intersect_eff_session$id_data, 
                                 len = st_length(intersect_eff_session),
                                 transect_name = intersect_eff_session$transect_name) %>% 
      group_by(id, transect_name) %>% 
      summarise(total_length = sum(len))
    
    occurence_df$effort <- 0
    occurence_df$effort[df_transect_length$id] <- df_transect_length$total_length
    occurence_df$transect_name <- NA
    occurence_df$transect_name[df_transect_length$id] <- df_transect_length$transect_name
    
    ### fill det/ non-det
    id_cells_with_obs <- intersect_obs %>% filter(session == S) %>% pull(id_data) %>% unique()
    occurence_df <- occurence_df %>% 
      mutate(y = ifelse(effort > 0, 0, NA)) %>% # fill cells with effort w/ 0
      mutate(y = ifelse(row_number() %in% id_cells_with_obs, 1, y)) %>% # fill cells with obs w/ 1
      mutate(y = ifelse(effort <= 0, NA, y)) # remove potential obs in cells absent from effort data this session
    
    effort_list = c(effort_list, occurence_df$effort)
    transect_name_list <- c(transect_name_list, occurence_df$transect_name)
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
  data$det.covs$transect_name <- matrix(transect_name_list, ncol = length(session_list))
  return(data)
}
