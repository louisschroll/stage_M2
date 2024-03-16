### V. Lauret

# Long piece of code to simulate tracks and fit RSF

##### packages ####
library(dagR)
library(tidyverse)
library(ggplot2)
library(spatialEco)
library(raster)
library(sf)
library(amt)
library(here)

library(localGibbs) # install from the Github of Theo Michelot
library(tidyverse)
library(stars)

# defined coolors palette  https://coolors.co/5f0f40-9a031e-fb8b24-e36414-0f4c5c
blue <-  "#3d405b"
orange <- "#f2cc8f"
red <- "#e07a5f"
purple <-  "#81b29a"

#### LANDSCAPE & TRANSECTS #####
# Create a landscape with 1 continuous cov and 1 discrete and spatial autocorrelation
# Landscape
# DISCRET
     covlist <- list()
     covlist[[1]] <- simRaster(rho=10, lim=c(-500,500,-500,500), res=1)
     covlist[[2]] <- simRaster(rho=30, lim=c(-500,500,-500,500), res=1)
     
     rr1 <- st_as_stars(covlist[[1]])
     rr2 <- st_as_stars(covlist[[2]])
     
     # discretize cov 1
     rr1$layer <- as_factor(ifelse(rr1$layer > 0.50, 1,0))
     # scale cov 2
     rr2.sc <- st_as_stars(scale(covlist[[2]]))


# plot the landscape      
if(FALSE){
  theme_set(theme_minimal())
  p_habd <- ggplot() +
    geom_stars(data= rr1) + 
    scale_fill_viridis_d(name = "Discrete", labels = c("A","B"))+
    theme(legend.position = 'top')
  
  p_habc <- ggplot() +
    geom_stars(data= rr2) +
    labs(fill = "Conitnuous")+
    theme(legend.position = 'top')
  library(patchwork)
  p_habd + p_habc
  }


# -- SIM TRACKS ----
    
# biased Correlated Random Walk towards colony

# Parameter settings
nindividual=1   # number of individuals
nstep=300       # number of steps by individual
p.sel=4         # preference of hab2 over hab1
s0=0.8          # param 1 for biased-CRW
k=0.4           # param 2 for biased-CRW

# TWO parameters
p.sel # parameter of preference for categorial veget 
# beta = 1.5 # RSF coefficient for continuous hab

ntmp = 10 # number of potential locations at each step
# Simulation of tracks with a biased CRW towards colony
data_tot=data.frame(matrix(NA,nrow=0,ncol=10))
colnames(data_tot)=c('id','time','x','y','dist2col','step','angle','direction','veget','cov')

# add individual heterogeneity
coef <- rnorm(nindividual, log(p.sel), 0.25) # individual heterogeneity 
# or no individual heterogeneity coef <- rep(log(p.sel), nindividual)
sd(coef)
# without heterogeneity sd(coef) = 0.17 over 1000 tracks. 

#### FUNCTION TO SIMULATE TRACKS 
for(j in 1:nindividual){ # j =1
  
  kcpt= coef[j] #rnorm(1, log(p.sel), 0.5)
  
  # Create dataset
  data = data.frame(matrix(NA,nrow=nstep,ncol=10))
  colnames(data) = c('id','time','x','y','dist2col','step','angle','direction','veget','cov')
  data$id = j
  data[1,]$time = 0
  data[1,]$x = 0
  data[1,]$y = 0
  data[1,]$dist2col = 0
  
  # First step from the colony
  
  # select next step
  #en faire une fonction
  tmp <- data.frame(matrix(NA,nrow=ntmp,ncol=7))
  colnames(tmp)=c('new.coord.x','new.coord.y','veg.tmp','cov.tmp', 'p.tmp','step','angle')
  
  for(i in 1:ntmp){
    tmp$step[i] = rgamma(n=1,shape=5,scale=2)
    tmp$angle[i] = runif(n=1,0,2*pi)
    tmp[i,1:2] = anglePoint(c(0,0), angl=tmp$angle[i]+pi, len=tmp$step[i])
    tmp$veg.tmp[i] = st_extract(rr1, cbind(tmp[i,1],tmp[i,2]))
    tmp$cov.tmp[i] = st_extract(rr2.sc, cbind(tmp[i,1],tmp[i,2]))
    # calculate a likelihood to select this location
    # discrete effect only
    tmp$p.tmp[i] <-  exp(kcpt * (as.numeric(tmp$veg.tmp[i][[1]])-1) ) # beta * tmp$cov.tmp[i][[1]] 
    
  }
  
  loc.f <- which(rmultinom(1,1, tmp$p.tmp) ==1)[1] # which(tmp$p.tmp == max(tmp$p.tmp))[1]
  
  # store step data
  data[2,]$time = 1
  data[2,]$x=tmp$new.coord.x[loc.f]
  data[2,]$y=tmp$new.coord.y[loc.f]
  data[2,]$dist2col=sqrt(tmp$new.coord.x[loc.f]^2+tmp$new.coord.y[loc.f]^2)
  data[2,]$step=tmp$step[loc.f]
  data[2,]$angle=tmp$angle[loc.f]
  data[2,]$direction=tmp$angle[loc.f]
  data[2,]$veget=tmp$veg.tmp[loc.f][[1]]
  data[2,]$cov=tmp$cov.tmp[loc.f][[1]]
  # if(data[2,]$veget==1){kcpt=kcpt-1}
  
  for(i in 3:nstep){
    # Dbefore=data[i-2,]$dist2col
    # D=data[i-1,]$dist2col
    # tmp <- loc.tmp(c(data$x[i-1], data$y[i-1]), ntmp)
    
    tmp <-data.frame(matrix(NA,nrow=ntmp,ncol=9))
    colnames(tmp)=c('new.coord.x','new.coord.y','veg.tmp','cov.tmp', 'p.tmp','step','angle','direction','dist2col')
    
    for( n in 1:ntmp){
      tmp$step[n]=rgamma(n=1,shape=5,scale=2)
      alpha=atan2(data[i-2,]$y-0,data[i-2,]$x-0) - atan2(data[i-2,]$y-data[i-1,]$y,data[i-2,]$x-data[i-1,]$x)
      alpha=2*((abs(alpha)-0)/(pi-0))-1
      angle=rnorm(n=1,mean=0,sd=s0*(1+k*alpha))
      angle=ifelse(angle>pi,angle-pi,ifelse(angle<(-pi),angle+pi,angle))
      tmp$angle[n]=angle
      tmp$direction[n]=data[i-1,]$direction+angle
      tmp$direction[n]=ifelse(tmp$direction[n]>2*pi,tmp$direction[n]-2*pi,ifelse(tmp$direction[n]<0,tmp$direction[n]+2*pi,tmp$direction[n]))
      tmp[n,1:2]=anglePoint(c(data[i-1,]$x,data[i-1,]$y),angl=tmp$direction[n]+pi,len=tmp$step[n])
      tmp$dist2col[n]=sqrt(tmp[n,1]^2+tmp[n,2]^2)
      tmp$veg.tmp[n] =st_extract(rr1,cbind(tmp[n,1],tmp[n,2]))
      tmp$cov.tmp[n] =st_extract(rr2.sc,cbind(tmp[n,1],tmp[n,2]))
      # calculate a likelihood to select this location
      # discrete effect onyl -- > 
      tmp$p.tmp[n] <- exp(kcpt * (as.numeric(tmp$veg.tmp[n][[1]])-1)) # beta * tmp$cov.tmp[i][[1]] +
      # continuous cov tmp$p.tmp[n] <- exp(beta * tmp$cov.tmp[n][[1]])
      
    }
    if(TRUE %in% is.na(tmp$p.tmp)){
      loc.f <-  which(is.na(tmp$p.tmp) ==F)[1] 
    }else{
    loc.f <-  which(rmultinom(1,1, na.omit(tmp$p.tmp)) ==1)[1] # which(tmp$p.tmp == max(tmp$p.tmp))[1]
    }
    
    data[i,]$step=tmp$step[loc.f]
    data[i,]$angle=tmp$angle[loc.f]
    data[i,]$direction=tmp$direction[loc.f]
    data[i,]$x=tmp[loc.f,1]
    data[i,]$y=tmp[loc.f,2]
    data[i,]$dist2col=tmp$dist2col[loc.f]
    data[i,]$time=i-1
    data[i,]$veget =tmp$veg.tmp[loc.f][[1]]
    data[i,]$cov =tmp$cov.tmp[loc.f][[1]]
    
  }
  
  
  data_tot=rbind(data_tot,data)
  
}

head(data_tot)
dim(data_tot) # nindividual * nstep
data_tot$id=factor(data_tot$id)

# plot the track
if(FALSE){
  p_trackd <- ggplot()+
    geom_stars(data= rr1, alpha = 0.6) +
    scale_fill_viridis_d(name ="Discrete cov")+
    geom_path(data = data_tot, aes(x, y, group = id, color = id), linewidth =2)+
    theme_classic()+
    theme(legend.position = 'top') +
    guides(color = "none") +
    labs(subtitle = paste0("Preference for hab = 1: p.sel = ", p.sel))
  
  # if continuous covariate used for simulation
  # p_trackc <- ggplot()+
  #   geom_stars(data= rr2, alpha = 0.9)+
  #   geom_path(data=data_tot,aes(x,y,group=id,color=id))+
  #   geom_point(aes(0,0),colour='black',show.legend = NA,inherit.aes=FALSE)+  # colony location
  #   theme_classic()+
  #   theme(legend.position = 'top') + 
  #   guides(color = "none") + 
  #   labs(subtitle = paste0("Preference for cov: beta = ", beta),
  #        fill = "Continuous cov")
  
  p_trackd #+ p_trackc
  
}

# 
p.est = table(data_tot$veget)[1] / table(data_tot$veget)[2]
-log(p.est)
log(p.sel)
# coefficient are not exactly the same --> Simulation stochasticity.

# --- SIM TRACKS 2 without attraction to colony ----
# nbObs <- 300 # number of locations to simulate
# beta <- c(4) # RSF coefficients of habitat covariates
# allr <- 10 # availability radius (affects perception and movement speed)
# 
# # simulations
# xy <- simLG(nbObs=nbObs, beta=beta, allr=allr, covlist=covlist[1])
# colnames(xy) <- c("x","y")
# 
# data_tot2 <- xy %>% as_tibble() %>% 
#   mutate(veget = st_extract(rr1,cbind(xy[,"x"],xy[,"y"]))) 
# 
# table(data_tot2$veget)

# ---- FIT RSF ----
head(data_tot)
nindividual

nind <- 1 # how many individuals
nptInd <- 300 # how many step per individual you want to keep
nptRand <- 3100 # how many random points



dfRSF <- tibble()
for(i in 1:nind){
  
  # subset true data for RSF 
  dat_g <- data_tot %>%
    mutate(id = as.numeric(id)) %>% 
    filter(id == i) %>% #sample(unique(data_tot$id),1)
    mutate(id = as_factor(i),
           dist = dist2col) 
  
  dat_samp <- dat_g[sample(nrow(dat_g),size = nptInd),] %>% 
    select(x,y,id,time,veget,dist) %>% 
    mutate(case = 1,
           veget = as_factor(veget - 1))
  
  # generat available data 
  meanDistanceColony <- 250
  nullCoord <- data.frame(
    dist = rexp(nptRand, 1/meanDistanceColony),
    angle = runif(nptRand, 0, 2*pi)
  ) %>% 
    mutate(
      x = cos(angle)*dist,
      y = sin(angle)*dist
    )
  
  nullCoord <- nullCoord %>% 
    as_tibble() %>% 
    arrange(dist) 
  
  nullCoord <- nullCoord[1:(0.95*nrow(nullCoord)),]
  
  rpts <- st_as_sf(nullCoord, coords= c("x","y")) %>% 
    st_crop(rr1)
  
  rpts2 <- rpts %>%
    mutate(veget = as_factor(st_extract(rr1, at = rpts)$layer),
           case = 0,
           id = as_factor(i),
           time = 0,
           x = st_coordinates(rpts)[,"X"],
           y = st_coordinates(rpts)[,"Y"]) %>% 
    select(veget,case,time,dist,x,y, id)
  
  dfRSF <- bind_rows(dfRSF, rpts2 %>% st_drop_geometry(), dat_samp)
  
}



head(dfRSF)

# description metrics
tab <- dfRSF %>% group_by(case) %>% count(veget)
tab

pest2 <- tab$n[tab$case==1 & tab$veget==0][1]/tab$n[tab$case==1 & tab$veget==1][1]

p.est = table(data_tot$veget)[1] / table(data_tot$veget)[2]
log(p.est)
log(4)
log(pest2)

# verify if points are outside study area
dfRSF$veget[is.na(dfRSF$veget)] <- 0
dfRSF$dist[is.na(dfRSF$dist)] <- 500
dfRSF$id %>% unique()
table(dfRSF %>% filter(case == 1) %>% pull(veget))

# --- RSF with AMT ----
library(amt)
w = dfRSF$case
w[w == 0] <- 1000
table(w)
#mrsf <- lme4::glmer(case ~ as_factor(veget) + (as_factor(veget)|id) ,data = dfRSF %>% filter(as.numeric(id) <=20), family = binomial(link = "logit"))
mod <- glmmTMB::glmmTMB(case ~ as_factor(veget), 
                        data = dfRSF,
                        family = binomial(link = "logit"), 
                        weights = w)
mod %>% summary()
#summary(mrsf)

## ----  RSF with NIMBLE ----
# Bayesian RSF w/ NIMBLE
library(nimble)

# rsf no random effect
rsf <- nimbleCode({
  
  ### PRIORS ###
  ## habitat cov
  beta_pop ~ dnorm(0,1)
  #betadist ~ dnorm(0,1)
  # intercept
  beta0  ~ dnorm(0,1)
  
  # no random effect  tauInd ~ dunif(0,1) # tau = precision = 1 / sd
  # no random effect  # individual random effects
  # no random effect  for( i in 1:nind){
  # no random effect    beta[i] ~ dnorm(beta_pop, tauInd)
  # no random effect  }
  
  # fit a logistic regression w/ fixed habitat param with random ind effect on slope
  for(t in 1:ntot){
    # logistic regression RSF
    # no random effect logit(omega[t]) <- beta0 + beta[idind[t]] * veget[t] + betadist * dist[t]
    logit(omega[t]) <- beta0 + beta_pop * veget[t] #+ betadist * dist[t]
    
    kase[t] ~ dbinom(omega[t],w[t])
  }
})
# rsf random effect 

head(dfRSF)

dfRSF$veget[is.na(dfRSF$veget)] <- 0

run.rsf <- function(dfRSF = dfRSF){
  
  w = dfRSF$case
  w[w == 0] <- 1000
# constants
constants.ni <-  list(ntot = nrow(dfRSF), 
                      nind = length(unique(dfRSF$id)),
                      veget = dfRSF$veget[],
                      idind = as.numeric(dfRSF$id),
                      w = w)#,
                     # dist = dfRSF$dist)

# data
data.ni <- list(kase = dfRSF$case)

# Inits
inits.ni <-  list(beta0 = 0, beta_pop = 0,
                  beta = rnorm(length(unique(dfRSF$id)),0,0.2),
                  tauInd = runif(1,0,1))# ,
                  #betadist = 0)

# --- bundle and run 
Rmodel2 <- nimbleModel(code= rsf, constants = constants.ni, data = data.ni, inits = inits.ni)
Rmodel2$initializeInfo()
print(Rmodel2$calculate())# - 10707
# configure model
conf2 <- configureMCMC(Rmodel2)
# conf2$printSamplers()
# conf2$removeSampler(target = c("beta0","betadist","beta_pop"))
# conf2$addSampler(target = c("beta0","betadist","beta_pop"), type= "RW_block")
#
## Build and compile MCMC
Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(Rmodel2)
Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)
#Cmcmc2$run(10)
# Run
samples <- runMCMC(Cmcmc2, niter = 50000, nburnin = 5000, nchains = 1, samplesAsCodaMCMC = TRUE) ## DT: use runMCMC
return(samples)
}

RSF <- run.rsf(dfRSF = dfRSF %>%
                   mutate(idd= as.numeric(id)) )

summary(RSF) # Nimble
mcmcplots::denplot(RSF)
coda::effectiveSize(RSF)


# plot
RSF %>% 
  as_tibble() %>% 
  mutate(model = "RSF") %>% 
  ggplot() + 
  geom_boxplot(aes(model, beta_pop), color = red) + 
  geom_hline(yintercept = - log(pest2), color = "darkgreen", lwd = 1, lty = 5)+
  labs(title = "Discrete habitat preference",
       subtitle = "green line is the true value") + 
  theme(plot.title = element_text(family = "Times", face = "bold"))


