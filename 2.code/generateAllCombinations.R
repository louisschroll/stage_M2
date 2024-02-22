

generateAllCombinations <- function(covars)  {
  # Function to generate all combinations of covariates.
  # covars : a character vector specifying the covariates/predictors
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





