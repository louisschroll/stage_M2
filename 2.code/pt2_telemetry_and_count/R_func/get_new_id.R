#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  get_new_id
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-10
#'
#' Script description:
#' 
#' This function takes a data frame with a column 'id', assign new id ranging 
#' from 1 to the number of distinct id and return a vector with the same length
#' and order than df_w_id$id_column with the 'new' id.
#' 
#' This function is used in both format_rsf_data_for_nimble and prepare_data_Nmix.
#' 
#' @df_w_id: a data frame with a column containing an numeric id
#' @id_column: a character, the name of the column with the id
#' 
#' @output: a numeric vector with the re numbered id in the right order. 
#' -------------------------------------------------------------------------------


get_new_id <- function(df_w_id, id_column){
  # Pull the id and create a data frame with the current and the new id
  sampled_sites_id <- df_w_id %>% pull(id_column) %>% unique() %>% sort()
  id_correspondance <- tibble(old_id = sampled_sites_id, new_id = 1:length(sampled_sites_id))
  
  # Put the new id in the good row by correspondence with old id and pull them
  new_sites_id <- df_w_id %>% select(id_column) %>% 
    left_join(id_correspondance, by = join_by({{ id_column }} == "old_id")) %>%
    pull(new_id) 
  return(new_sites_id)
}