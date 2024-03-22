# HEADER ------------------------------------------------------------------------
#
# Script name:  model_selection_functions.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-16
#
# Script description:
# Functions to run a model with spOccupancy (intPGOcc()) without spatial 
# auto-correlation, and to retrieve the statistics associated with this model:
# Rhat, ESS, bayesian p-values, waic score, and cross validation score
#
# ------------------------------------------------------------------------------


test_all_models <- function(cov_combination, data.int){
  model_nb <- length(cov_combination)
  model_names <- map(cov_combination, function(x) paste(x, collapse = " + ")) %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(model = value)
  
  # Run models with spOccupancy
  model_result_list <- map(cov_combination, function(x){run_int_model(data.int, selected_cov = x)})
  
  # Retrieve statistics (Rhat, ESS, pvalues, WAIC, CV) for each model
  comparison_df <- model_names %>%
    bind_cols(map_dfr(model_result_list, ~{
                 get_model_stat(model_result = .x)}))
  
  # for (i in seq_along(cov_combination)){
  #   print(paste("Testing model", i, "/", model_nb))
  #   selected_cov = cov_combination[[i]]
  #   model_result <- run_int_model(data.int, selected_cov)
  #   comparison_df <- comparison_df %>% bind_rows(get_model_stat(model_result))
  #   beta_value_df <- beta_value_df %>% bind_rows(get_beta_values(model_result,selected_cov, model_nb=i))
  # }
  # 
  # comparison_df <- map(cov_combination, function(x) paste(x, collapse = " + ")) %>% 
  #   unlist() %>% 
  #   as_tibble() %>% 
  #   rename(model = value) %>% 
  #   bind_cols(comparison_df) 
  
  if (any(unlist(cov_combination) != "1")){
    # beta_df_completed <- complete(beta_value_df, covar, model, fill = list(beta = NA)) %>% 
    #   pivot_wider(names_from = covar, values_from = c(beta, sd_beta), names_sep = "_") %>% 
    #   select(-model)
  # Retrieve values of coefficients beta
  beta_value_df <- map_dfr(model_result_list, ~{
    model_result <- .x
    selected_cov <- model_result$beta.samples %>% as_tibble() %>% colnames()
    get_beta_values(model_result, selected_cov[-1], 
                    model_nb=paste(selected_cov[-1], 
                                   collapse ="+"))}) %>% 
    complete(covar, model, fill = list(beta = NA)) %>% 
    pivot_wider(names_from = covar, values_from = c(beta, sd_beta), names_sep = "_") %>% 
    select(-model)
  
    comparison_df <- bind_cols(comparison_df, beta_value_df)
  }
  return(list(comparison_df = comparison_df, model_result_list = model_result_list))
}


run_int_model <- function(data.int, selected_cov, add_spatial=FALSE, spatial_model="exponential"){
  # Wrapper for intPGOcc() function of spOccupancy
  # data.int: a list containing the data with the correct format for intPGOcc()
  # selected_cov: a character vector with the covariates to include in the model

  data_copie <- data.int
  occ.formula <- writeFormula(selected_cov)
  
  nb_datasets <- length(data_copie$y)
  total_sites_nb <- nrow(data.int$occ.covs)
  det.formula <- as.list(rep("~ scale(transect_length) + session", nb_datasets)) %>%  
    map(as.formula)
  
  if (!add_spatial){
    inits.list <- list(alpha = as.list(rep(0, nb_datasets)),
                       beta = 0, 
                       z = rep(1, total_sites_nb))
    
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                       alpha.normal = list(mean = as.list(rep(0, nb_datasets)), 
                                           var = as.list(rep(2.72, nb_datasets))))
    
    n.samples <- 9000
    n.burn <- 1500
    n.thin <- 1
    
    model_result <- intPGOcc(occ.formula = occ.formula,
                             det.formula = det.formula, 
                             data = data_copie,
                             inits = inits.list,
                             n.samples = n.samples, 
                             priors = prior.list, 
                             n.omp.threads = 1, 
                             verbose = FALSE, 
                             n.report = 2000, 
                             n.burn = n.burn, 
                             n.thin = n.thin, 
                             n.chains = 2,
                             k.fold = 6, 
                             k.fold.threads = 6,
                             k.fold.only = FALSE,
                             k.fold.seed = 42)
  }
  else {
    dist.int <- dist(data.int$coords)
    min.dist <- min(dist.int)
    max.dist <- max(dist.int)

        inits.list <- list(alpha = as.list(rep(0,nb_datasets)),
                       beta = 0, 
                       z = rep(1, total_sites_nb), 
                       sigma.sq = 2,
                       phi = 3 / mean(dist.int), 
                       w = rep(0, total_sites_nb))
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                       alpha.normal = list(mean = as.list(rep(0, nb_datasets)), 
                                           var = as.list(rep(2.72, nb_datasets))), 
                       sigma.sq.ig = c(2, 1), 
                       phi.unif = c(3 / max.dist, 3 / min.dist))
    
    batch.length <- 25
    n.batch <- 360
    n.burn <- 1500
    n.thin <- 3
    tuning <- list(phi = .2)
    model_result <- spIntPGOcc(occ.formula = occ.formula, 
                             det.formula = det.formula, 
                             data = data.int, 
                             inits = inits.list, 
                             priors = prior.list, 
                             tuning = tuning, 
                             cov.model = spatial_model, 
                             NNGP = T, 
                             n.neighbors = 5, 
                             n.batch = n.batch, 
                             n.burn = n.burn, 
                             n.chains = 3,
                             batch.length = batch.length, 
                             n.report = 200,
                             k.fold = 6, 
                             k.fold.threads = 6,
                             k.fold.only = FALSE,
                             k.fold.seed = 42) 
  }
  return(model_result)
}


run_non_int_model <- function(data, selected_cov){
  # Wrapper for intPGOcc() function of spOccupancy
  # data: a list containing the data with the correct format for intPGOcc()
  # selected_cov: a character vector with the covariates to include in the model
  
  occ.formula <- writeFormula(selected_cov)
  
  total_sites_nb <- nrow(data$occ.covs)
  det.formula <- ~ scale(transect_length) + session
  
  n.samples <- 9000
  n.burn <- 1500
  n.thin <- 3
    
  data <- list(y = data$y[[1]], 
                 occ.covs = data$occ.covs,
                 det.covs = data$det.covs[[1]])
    
  inits.list <- list(alpha = 0, 
                       beta = 0, 
                       z = apply(data$y, 1, max, na.rm = TRUE))
    
  prior.list <- list(alpha.normal = list(mean = 0, var = 2.72), 
                       beta.normal = list(mean = 0, var = 2.72))
    
    
  model_result <- PGOcc(occ.formula = occ.formula,
                             det.formula = det.formula, 
                             data = data,
                             inits = inits.list,
                             n.samples = n.samples, 
                             priors = prior.list, 
                             n.omp.threads = 1, 
                             verbose = FALSE, 
                             n.report = 2000, 
                             n.burn = n.burn, 
                             n.thin = n.thin, 
                             n.chains = 2,
                             k.fold = 6, 
                             k.fold.threads = 6,
                             k.fold.only = FALSE,
                             k.fold.seed = 42)

  return(model_result)
}


writeFormula <- function(covars){
  # covars : a character vector specifying the covariates/predictors
  as.formula(paste("~ ", paste(covars, collapse = " + ")))
}


get_model_stat <- function(model_result){
  # retrieve the statistics used to assess model performance
  # model_result: output of intPGOcc() function (spOccupancy)
  # Output: a tibble of dimension 20x1 storing max value of Rhat, min value
  # of ESS, bayesian p-values for each dataset (group 1 & 2), WAIC and 
  # k-fold deviance for each dataset and total
  
  convergence_df <- tibble(rhat_max = NA, ESS_min=NA)
  # Convergence: Rhat & ESS
  convergence_df$rhat_max <- model_result$rhat %>% unlist() %>% max(na.rm=T) %>% round(3)
  convergence_df$ESS_min <- model_result$ESS %>% unlist() %>% min(na.rm=T) %>% round(1)
  
  # Bayesian p-value
  ppc.out.1 <- ppcOcc(model_result, 'chi-squared', group = 1)
  ppc.out.2 <- ppcOcc(model_result, 'chi-squared', group = 2)
  
  len <- length(ppc.out.1$fit.y)
  nb_datasets <- ifelse(len<=4, len, 1)
  pvalue_df <- tibble(gr1 = compute_bayesian_pvalue(ppc.out.1), 
                      gr2 = compute_bayesian_pvalue(ppc.out.2)) %>% 
    pivot_longer(everything()) %>% 
    mutate(name = paste0("pvalue_d", rep(1:nb_datasets, each=2), "_", name)) %>% 
    pivot_wider(values_from = value, names_from = name)
    
  # WAIC
  waic_df <- waicOcc(model_result) %>% 
    bind_rows() %>% 
    select(WAIC) %>% 
    round(0) %>% 
    mutate(name = paste0("waic_d",1:nrow(.))) %>% 
    pivot_wider(names_from = name, values_from = WAIC) %>% 
    mutate(waic_tot = rowSums(.))
  
  # k-fold cross-validation
  kfold_df <- tibble(deviance = model_result$k.fold.deviance %>% round(0), 
                  name = paste0("CV_d", seq_along(model_result$k.fold.deviance))) %>% 
    pivot_wider(names_from = name, values_from = deviance) %>% 
    mutate(CV_tot = rowSums(.))
  
  return(cbind(convergence_df, pvalue_df, waic_df, kfold_df))
}


compute_bayesian_pvalue <- function(ppc.out){
  # Compute baysian p-values with the result of ppcOcc function (spOccupancy)
  # Output: a numeric vector containing one p-value for each data set
  if (length(ppc.out$fit.y)>1 & length(ppc.out$fit.y)<=4){
    obs_values <- ppc.out$fit.y
    sim_values <- ppc.out$fit.y.rep
    pvalues <- lapply(seq_along(obs_values), function(i) {
      obs <- obs_values[[i]]
      sim <- sim_values[[i]]
      sum(sim > obs) / length(obs)
    })
    return(round(unlist(pvalues), 2))
  }
  else {
    obs <- ppc.out$fit.y
    sim <- ppc.out$fit.y.rep
    pvalue <- sum(sim > obs) / length(obs)
    return(round(pvalue, 2))
  }
}


get_beta_values <- function(model_result, selected_cov, model_nb=1){
  if ("1" %in% selected_cov){
    return(tibble())
  }
  selected_cov[str_detect(selected_cov, pattern = "I")] <- str_replace(selected_cov[str_detect(selected_cov, pattern = "I")], "\\^2", "")
  
  mean_beta <- model_result$beta.samples %>% as_tibble() %>% 
    select(all_of(selected_cov)) %>% 
    summarise_all(mean) %>% 
    pivot_longer(all_of(selected_cov), names_to = "covar", values_to = "beta") %>% 
    mutate(model = model_nb,
           beta = round(beta, 2)) 
  
  sd_beta <- model_result$beta.samples %>% as_tibble() %>% 
    select(all_of(selected_cov)) %>% 
    summarise_all(sd) %>% 
    pivot_longer(all_of(selected_cov), names_to = "covar", values_to = "sd_beta") %>% 
    mutate(model = model_nb,
           sd_beta = round(sd_beta, 2))
  
  return(mean_beta %>% 
    full_join(sd_beta, by = join_by(covar, model)))
}


generateAllCombinations <- function(covars)  {
  # Function to generate all possible combinations of covariates.
  # covars: a character vector specifying the covariates/predictors
  if(!is.character(covars))
    stop("covars must be character vectors")
  covars <- unique(covars)
  
  ncovs <- length(covars)
  nformulas <- 2^ncovs
  tfmat <- matrix(FALSE, nformulas, ncovs) 
  for(i in 1:ncovs){
    tfmat[, i] <- rep(c(FALSE, TRUE), each=2^(i-1))
  }
  if(ncovs > 1)
    tfmat <- tfmat[order(rowSums(tfmat)), ]
  RHS <- apply(tfmat, 1, function(x) covars[x])
  RHS[1] <- "1"
  return(RHS)
}


addSelectionSheet <- function(workbook, sheet_name, df, datasets_nb){
  # Add the data in df to a new sheet of a workbook
  # workbook: a Workbook object, openxlsx package
  # sheet_name: a character giving the name of the sheet
  # df: a dataframe storing the results of model selection (formulas, 
  # rhat, ESS, WAIC, cross-validation)
  addWorksheet(workbook, sheet_name)
  writeData(workbook, sheet = sheet_name, x = df)
  style <- createStyle(valign = "center", halign = "center")
  addStyle(workbook, sheet = sheet_name, style = style, rows = 1:(nrow(df)+1), cols = 1:ncol(df), gridExpand = TRUE)
  for (cols in 1:(2*datasets_nb)+3) {
    conditionalFormatting(workbook, sheet = sheet_name, cols = cols, rows = 2:(nrow(df) + 1),
                          type = "between", rule = c(0.1, 0.9), 
                          style = createStyle(bgFill = "lightgreen", fontColour = "darkgreen"))
  }
  
  for (cols in 1:(2*(datasets_nb+1))+2*datasets_nb+3) {
    conditionalFormatting(workbook, sheet = sheet_name, cols = cols, rows = 2:(nrow(df) + 1),
                          type = "colorScale", style = c("red", "white"))
    min_value <- df %>% pull(all_of(cols)) %>% min()
    conditionalFormatting(workbook, sheet = sheet_name, cols = cols, rows = 2:(nrow(df) + 1),
                          type = "expression", 
                          rule = paste0("<", min_value+3),
                          style = createStyle(fontColour = "lightpink", textDecoration = "bold"))
  }
}


get_dataset_names <- function(data_list, species){
  # Get the names of datasets containing obsevations of the species
  # data_list: a list of the data sets, each subdivided in a sub-list of effort 
  # and observation
  # species: a character
  # output: a character vector indicating data sets names
  map(data_list, ~ filter_one_dataset(.x, species)) %>% 
    compact() %>% 
    names()
}


save_maps_as_pdf <- function(maps_list, path_to_save){
  pdf(path_to_save)
  for (i in seq_along(maps_list)){
    plot <- plot_grid(plotlist = maps_list[[i]], nrow = 2,
                      labels = names(maps_list)[i], #%>% str_to_upper(), 
                      label_y = 1,
                      scale = 0.92)
    print(plot)
  }
  dev.off() 
}
