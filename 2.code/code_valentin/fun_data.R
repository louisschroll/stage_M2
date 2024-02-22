# R functions to prepare data for spatial intergrated occupancy 

# find dates of sampling 
samp_dates <- function(do.pel = T, do.sam = F, do.pnm = T, do.mig = T, do.count = T){
  dates <- list()
  
  if(do.pel){
    dates$peldate <- range(pelmed_eff$date)}
  
  if(do.pnm){
    dates$pnmdate <- range(transect$date_tr)
  }
  
  if(do.sam){
    dates$samdate <- range(effsamm$date)
  }
  
  if(do.mig){
    dates$ migdate <- range(c(postnup2022_eff$day, prenup$Date_UTC))
  }
  
  if(do.count){
    dates$count <- unique(c(countcol$Year))
  }
  
  return(dates)
}


# find which sites are sampled form the grid 
sites.l <- function(do.pel = T, do.sam = T,
                    do.pnm = T, do.mig = T,
                    grid = grid){
  # find the grid-cells that intersect with sampling effort
  sites <- list()
  
  # pelmed
  if(do.pel == T){
    intpel <- st_intersection(pelmed_eff, grid) 
    sites$pelmed <- unique(intpel$id)
  }
  
  # pnmgdl
  if(do.pnm == T){
    intpnm <- st_intersection(transect, grid) 
    sites$pnm <- unique(intpnm$id.1)
  }
  
  # migralion
  if(do.mig == T){
    intmigr <- st_intersection(postnup2022_eff, grid) 
    sites$migralion <- unique(intmigr$id)
  }
  
  # migralion
  if(do.sam == T){
    # remove surveys from 2019 as no detection has been recorded
    effsamm1 <- effsamm %>% 
      mutate(year = lubridate::year(date)) %>% 
      st_transform(crs= st_crs(grid)) %>% 
      filter(year %in% c("2011","2012"))
    
    intsamm <- st_intersection(effsamm1, grid) 
    sites$samm <- unique(intsamm$id)
    sites$samm <- sites$samm[-c(1232:1235)]
  }
  
  
  IDs <- unique(unlist(sites))
  gridocc <- grid[IDs,] %>% 
    mutate(ido = 1:length(IDs))
  
  # gridocc stores the sites sampled at least once by a dataset
  sitesocc <- list()
  
  # pelmed
  if(do.pel == T){
    intpel <- st_intersection(pelmed_eff, gridocc) 
    sitesocc$pelmed <- unique(intpel$ido)
  }
  
  # samm
  if(do.sam == T){
    intsamm<- st_intersection(effsamm1, gridocc) 
    sitesocc$samm <- unique(intsamm$ido)
  }
  
  # migralion
  if(do.mig == T){
    intmigr<- st_intersection(postnup2022_eff, gridocc) 
    sitesocc$migralion <- unique(intmigr$ido)
  }
  
  # pnmgdl
  if(do.pnm == T){
    intpnm <- st_intersection(transect, gridocc) 
    sitesocc$pnm <- unique(intpnm$ido)
  }
  
  # export the restricted grid
  sitesocc$grid <- gridocc
  
  return(sitesocc) 
  # #SITESOCC is the 'sites' list we needed 
}


# data and covariates 
yndet <- function(do.pel = T, do.sam = T,
                  do.pnm = T, do.mig = T){
  
  yocc <- list()
  detcov <- list()
  yy <- list()
  pp <- list()
  
  
  ## ---- pelmed -----
  if(do.pel == T){
    pelmed <- pelmed %>% 
      mutate(year = lubridate::year(date)) 
    
    ypel <- gridocc[sitesocc$pelmed,] %>% 
      mutate(idpel = 1:length(sitesocc$pelmed),
             ytot = 0,
             y2017 = NA,
             y2018 = NA,
             y2019 = NA,
             y2020 = NA,
             efftot = 0,
             eff2017 = 0,
             eff2018 = 0,
             eff2019 = 0,
             eff2020 = 0)
    intypel <- st_intersection(pelmed, ypel) 
    inteffpel <- st_intersection(pelmed_eff, ypel) 
    inteffpel$length <- st_length(inteffpel)
    
    ### fill seff 
    for(i in 1:nrow(ypel)){
      id <- which(inteffpel$ido == ypel$ido[i])
      ypel$efftot[i] <- sum(inteffpel$length[id])
      ypel$eff2017[i] <- sum(inteffpel$length[id][year(inteffpel$date)[id] =="2017"])
      ypel$eff2018[i] <- sum(inteffpel$length[id][year(inteffpel$date)[id] =="2018"])
      ypel$eff2019[i] <- sum(inteffpel$length[id][year(inteffpel$date)[id] =="2019"])
      ypel$eff2020[i] <- sum(inteffpel$length[id][year(inteffpel$date)[id] =="2020"])
      
    }
    
    ### fill det/ non-det
    
    #### fill w/ 0
    ypel$y2017[ypel$eff2017 > 0] <- 0
    ypel$y2018[ypel$eff2018 > 0] <- 0
    ypel$y2019[ypel$eff2019 > 0] <- 0
    ypel$y2020[ypel$eff2020 > 0] <- 0
    
    #### fill w/ 1
    ypel$ytot[unique(intypel$idpel)] <- 1
    ypel$y2017[unique(intypel$idpel[intypel$year=="2017"])] <- 1
    ypel$y2018[unique(intypel$idpel[intypel$year=="2018"])] <- 1
    ypel$y2019[unique(intypel$idpel[intypel$year=="2019"])] <- 1
    ypel$y2020[unique(intypel$idpel[intypel$year=="2020"])] <- 1
    
    a <- matrix(c(ypel$y2017,
                  ypel$y2018,
                  ypel$y2019,
                  ypel$y2020), nrow = nrow(ypel), ncol = 4)
    
    pa <- matrix(c(ypel$eff2017,
                   ypel$eff2018,
                   ypel$eff2019,
                   ypel$eff2020),
                 nrow = nrow(ypel), ncol = 4)
    pa[pa ==0] <- NA
    pa <- scale(pa)
    
    # check dimensions 
    diff.na <-  length(which(is.na(a)==T)) - length(which(is.na(pa)==T))
    
    if( !(length(which(is.na(a))) == length(which(is.na(pa))) ) ){
      print("Warnings: inconsistent dimensions in pelmed data")
    }
    if( diff.na < 0){
      id <- setdiff(which(is.na(pa)) ,which(is.na(a)))
      pa[id] <-  0
    }
    if( diff.na > 0){
      id <- setdiff(which(is.na(a)) ,which(is.na(pa)))
      a[id] <-  0
    }
    if( length(which(is.na(a))) == length(which(is.na(pa)) ) ){
      print("Warnings: dimensions corrected")
    }
    
    # return pelmed matrix
    yy$a <- a
    pp$pa <- pa
    
  }
  ## ----  samm -----
  
  if(do.sam == T){
    
    samm <- samm %>% 
      mutate(year = lubridate::year(date)) %>% 
      st_transform(crs = st_crs(gridocc))
    
    # effsamm1
    effsamm1 <- effsamm %>% 
      mutate(year = lubridate::year(date)) %>% 
      st_transform(crs = st_crs(grid)) %>% 
      filter(year %in% c("2011","2012"))
    
    ysamm <- gridocc[sitesocc$samm,] %>% 
      mutate(idsamm = 1:length(sitesocc$samm),
             ytot = 0,
             y2011 = NA,
             y2012 = NA,
             efftot = 0,
             eff2011 = 0,
             eff2012 = 0)
    intysamm <- st_intersection(samm, ysamm) 
    inteffsamm <- st_intersection(effsamm1, ysamm) 
    inteffsamm$length <- st_length(inteffsamm)
    
    ### fill seff 
    for(i in 1:nrow(ysamm)){
      id <- which(inteffsamm$ido == ysamm$ido[i])
      ysamm$efftot[i] <-  sum(inteffsamm$length[id])
      ysamm$eff2011[i] <- sum(inteffsamm$length[id][year(inteffsamm$date)[id] =="2011"])
      ysamm$eff2012[i] <- sum(inteffsamm$length[id][year(inteffsamm$date)[id] =="2012"])
      
    }
    
    ###fill y
    #### fill w/ 0
    ysamm$y2011[ysamm$eff2011 > 0] <- 0
    ysamm$y2012[ysamm$eff2012 > 0] <- 0
    
    #### fill w/ 1
    ysamm$ytot[unique(intysamm$idsamm)] <- 1
    ysamm$y2011[unique(intysamm$idsamm[intysamm$year=="2011"])] <- 1
    ysamm$y2012[unique(intysamm$idsamm[intysamm$year=="2012"])] <- 1
    
    b <- matrix(c(ysamm$y2011,
                  ysamm$y2012),
                nrow = nrow(ysamm), ncol = 2)
    
    pb <- matrix(c(ysamm$eff2011,
                   ysamm$eff2012),
                 nrow = nrow(ysamm), ncol = 2)
    
    pb[pb==0] <- NA
    pb <- scale(pb)
    
    # check dimensions 
    diff.na <-  length(which(is.na(b)==T)) - length(which(is.na(pb)==T))
    
    if( !(length(which(is.na(b))) == length(which(is.na(pb))) ) ){
      print("Warnings: inconsistent dimensions in samm data")
    }
    if( diff.na < 0){
      id <- setdiff(which(is.na(pb)) ,which(is.na(b)))
      pb[id] <-  0
    }
    if( diff.na > 0){
      id <- setdiff(which(is.na(b)) ,which(is.na(pb)))
      b[id] <-  0
    }
    if( length(which(is.na(b))) == length(which(is.na(pb)) ) ){
      print("Warnings: dimensions corrected")
    }
    
    # return samm matrix
    yy$b <- b
    pp$pb <- pb
  }
  
  ## ---- migralion -----
  
  if(do.mig == T){
    migralion <- migralion %>% 
      mutate(year = lubridate::year(Date_UTC),
             month = lubridate::month(Date_UTC)) 
    
    ymigr <- gridocc[sitesocc$migralion,] %>% 
      mutate(idmigr = 1:length(sitesocc$migralion),
             ytot = 0,
             y2022pre = NA,
             y2022post = NA,
             efftot = 0,
             eff2022pre = 0,
             eff2022post = 0)
    intymigr <- st_intersection(migralion, ymigr) 
    inteffmigr <- st_intersection(postnup2022_eff, ymigr) 
    inteffmigr$length <- st_length(inteffmigr)
    
    ### fill seff 
    for(i in 1:nrow(ymigr)){
      id <- which(inteffmigr$ido == ymigr$ido[i])
      ymigr$efftot[i] <-  sum(inteffmigr$length[id])
      ymigr$eff2022post[i] <- sum(inteffmigr$length[id])
      
    }
    
    ### fill y 
    #### fill w/ 0
    ymigr$y2022pre[ymigr$eff2022post > 0] <- 0
    ymigr$y2022post[ymigr$eff2022post > 0] <- 0
    
    #### fill w/ 1
    ymigr$ytot[unique(intymigr$idmigr)] <- 1
    ymigr$y2022pre [unique(intymigr$idmigr[intymigr$month %in% c(3,4,5)])] <- 1
    ymigr$y2022post[unique(intymigr$idmigr[intymigr$month %in% c(9,10)])] <- 1
    
    c <- matrix(c(ymigr$y2022pre,
                  ymigr$y2022post),
                nrow = nrow(ymigr), ncol = 2)
    
    pc <- matrix(c(ymigr$eff2022pre,
                   ymigr$eff2022post),
                 nrow = nrow(ymigr), ncol = 2)
    
    pc[pc==0] <- NA
    pc <- scale(pc)
    pc[,1] <- pc[,2] # same sampling effort between prenup and postnup
    
    # check dimensions 
    diff.na <-  length(which(is.na(c)==T)) - length(which(is.na(pc)==T))
    
    if( !(length(which(is.na(c))) == length(which(is.na(pc))) ) ){
      print("Warnings: inconsistent dimensions in migralion data")
    }
    if( diff.na < 0){
      id <- setdiff(which(is.na(pc)) ,which(is.na(c)))
      pc[id] <-  0
    }
    if( diff.na > 0){
      id <- setdiff(which(is.na(c)) ,which(is.na(pc)))
      c[id] <-  0
    }
    if( length(which(is.na(c))) == length(which(is.na(pc)) ) ){
      print("Warnings: dimensions corrected")
    }
    
    # return migralion matrix
    yy$c <- c
    pp$pc <- pc
  }
  
  ## ---- pnm ----
  if(do.pnm == T){
    pnm <- pnm %>% 
      mutate(year = lubridate::year(dmy(date)),
             month = lubridate::month(dmy(date))) %>% 
      st_transform(crs = st_crs(gridocc))
    
    ypnm <- gridocc[sitesocc$pnm,] %>% 
      mutate(idpnm = 1:length(sitesocc$pnm),
             ytot = 0,
             y2019 = NA,
             y2020 = NA,
             y2021 = NA,
             efftot = 0,
             eff2019 = 0,
             eff2020 = 0,
             eff2021 = 0)
    intypnm <- st_intersection(pnm, ypnm) 
    inteffpnm <- st_intersection(transect, ypnm) 
    inteffpnm$length <- st_length(inteffpnm)
    
    ### fill seff 
    for(i in 1:nrow(ypnm)){
      id <- which(inteffpnm$ido == ypnm$ido[i])
      ypnm$efftot[i] <-  sum(inteffpnm$length[id])
      ypnm$eff2019[i] <- sum(inteffpnm$length[id][year(inteffpnm$date)[id] =="2019"])
      ypnm$eff2020[i] <- sum(inteffpnm$length[id][year(inteffpnm$date)[id] =="2020"])
      ypnm$eff2021[i] <- sum(inteffpnm$length[id][year(inteffpnm$date)[id] =="2021"])
    }
    
    ## fill y 
    #### fill w/ 0
    ypnm$y2019[ypnm$eff2019 > 0] <- 0
    ypnm$y2020[ypnm$eff2020 > 0] <- 0
    ypnm$y2021[ypnm$eff2021 > 0] <- 0
    
    #### fill w/ 1
    ypnm$ytot[unique(intypnm$idpnm)] <- 1
    ypnm$y2019[unique(intypnm$idpnm[intypnm$year=="2019"])] <- 1
    ypnm$y2020[unique(intypnm$idpnm[intypnm$year=="2020"])] <- 1
    ypnm$y2021[unique(intypnm$idpnm[intypnm$year=="2021"])] <- 1
    
    d <- matrix(c(ypnm$y2019,
                  ypnm$y2020,
                  ypnm$y2021),
                nrow = nrow(ypnm), ncol = 3)
    
    pd <- matrix(c(ypnm$eff2019,
                   ypnm$eff2020,
                   ypnm$eff2021),
                 nrow = nrow(ypnm), ncol = 3)
    
    pd[pd==0] <- NA
    pd <- scale(pd)
    
    # check dimensions 
    diff.na <-  length(which(is.na(d)==T)) - length(which(is.na(pd)==T))
    
    if( !(length(which(is.na(d))) == length(which(is.na(pd))) ) ){
      print("Warnings: inconsistent dimensions in samm data")
    }
    if( diff.na < 0){
      id <- setdiff(which(is.na(pd)) ,which(is.na(d)))
      pd[id] <-  0
    }
    if( diff.na > 0){
      id <- setdiff(which(is.na(d)) ,which(is.na(pd)))
      d[id] <-  0
    }
    if( length(which(is.na(d))) == length(which(is.na(pd)) ) ){
      print("Warnings: dimensions corrected")
    }
    
    # return pnm matrix
    yy$d <- d
    pp$pd <- pd
  }
  
  det.covv <- list()
  ndim <- length(yy)
  
  for(i in 1:ndim){
    det.covv[[i]] <- list(matrix(as.numeric(pp[[i]]), 
                                 nrow = nrow(pp[[i]]),
                                 ncol = ncol(pp[[i]])))
    
    names(det.covv[[i]]) <- paste0("det.cov.",i)
  }
  
  # return data 
  dat <- list(yy= yy,pp= det.covv)
  
  return(dat)
}


# sst <- sites.l(do.pel = T, do.sam = F, do.pnm = T, do.mig = T)
# gridocc <- sst$grid
# sitesocc <- sst
# sitesocc$grid <- NULL
# COV <- gridocc$depth.sc
# xoords <- st_coordinates(st_centroid(gridocc))
# ydt <- yndet(do.pel = T, do.sam = F,
#              do.pnm = T, do.mig = T)
# 