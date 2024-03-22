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




add_session_column <- function(data, min_time_btw_session = 30) {
  data <- data[order(data$date), ]
  
  # Initialize session counter and session labels
  session_counter <- 0
  session_labels <- c("A")
  
  # Initialize vector to store session labels
  session_column <- rep(NA, nrow(data))
  
  for (i in 1:nrow(data)) {
    # Check if it's time to start a new session
    if (i==1 || difftime(data$date[i], data$date[i - 1], units = "days") > min_time_btw_session) {
      session_counter <- session_counter + 1
      session_labels <- LETTERS[session_counter]
    }
    session_column[i] <- session_labels
  }
  
  # Add the column
  data$session <- as.factor(session_column)
  
  return(data)
}




addQuadraticEffect <- function(cov_combi_list, quadratic_cov){
  # Add the models with a quadratic effect on a selected covariates to the models
  # list
  new_combination <- map(cov_combi_list, function(x) x <- c(x, paste0("I(",quadratic_cov,")^2"))) %>% 
    purrr::keep(function(x) quadratic_cov %in% x)
  
  return(c(cov_combi_list, new_combination))
}




# Define UI for application
# ui <- fluidPage(
#   titlePanel("Choose Graph"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       selectInput("graph_type", "Select Graph Type:",
#                   choices = c("Traceplots occurence coeffs" = "TO",
#                                "Density occurence coeffs" = "DO",
#                                "Traceplots detection coeffs" = "TD",
#                                "Density detection coeffs" = "DD",
#                               "Chi-squared Group 1" = "CS1",
#                               "Chi-squared Group 2" = "CS2",
#                               "Freeman-Tukey Group 1" = "FT1",
#                               "Freeman-Tukey Group 2" = "FT2"))
#     ),
#     
#     mainPanel(
#       plotOutput("selected_plot")
#     )
#   )
# )
# 
# # Define server logic
# server <- function(input, output) {
#   output$selected_plot <- renderPlot({
#     if (input$graph_type == "CS1") {
#       res_fou$PPC_plots$pvalue_CS1
#     } else if (input$graph_type == "CS2") {
#       res_fou$PPC_plots$pvalue_CS2
#     } else if (input$graph_type == "FT1") {
#       res_fou$PPC_plots$pvalue_FT1
#     } else if (input$graph_type == "FT2") {
#       res_fou$PPC_plots$pvalue_FT2
#     } else if (input$graph_type == "TO"){
#       res_fou$coeff_plot$beta_traceplot
#     } else if (input$graph_type == "DO"){
#       res_fou$coeff_plot$beta_density
#     } else if (input$graph_type == "TD"){
#       res_fou$coeff_plot$alpha_traceplot
#     } else if (input$graph_type == "DD"){
#       res_fou$coeff_plot$alpha_density
#     }
#   })
# }
# 
# # Run the application 
# shinyApp(ui = ui, server = server)
