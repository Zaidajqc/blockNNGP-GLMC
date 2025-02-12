############################################################################
## Date: 12.02.2025
##
## Description:
##
##    The code performs a Bayesian analysis of generalized linear  
##    model of coregionalization using block-NNGP and NNGP through 
##    Integrated Nested Laplace approximation (INLA). These models 
##    are applied to Birds sampling data. It is assumed a Poisson 
##    distribution for the response variables (bird abundances of  
##    American robin and Mourning dove). 
##    @author:  Zaida Quiroz (\email{zquiroz@pucp.edu.pe}).
##
#############################################################################

gc()
rm(list=ls())


#######################
## libraries required
####################### 
library(INLA)
library(fields)
library(lattice)
library(akima) 
library(Matrix)
library(slam)
library(igraph)
library(coda)
library(MBA)
library(mvtnorm)
library(ggforce)
library(Rcpp)
library(tidyverse)
library(raster)
library(tidyverse)
library(patchwork)
library(FNN)

#######################
## functions required
####################### 

setwd("/set_your_directory/")
# setwd("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/FINALCODES-Multivariate_blockNNGP/Github/App_Birds/")

source("blockNNGPfunctions.R")
source("blockNNGPrgeneric.R")
source("NNGPrgeneric.R")
source("utils.R")
source("summary_functions.R")

## Description:
##
## runblockNNGP() function
## Arguments:
## data.est: data for estimation, data.est=data.frame(loc=loc, y1=y1,y2=y2, X=X)
##    loc: locations i.e. (long, lat), y1: response variable 1, 
##    y2: response variable 2, X: matrix of covariates. 
## data.pred: data for prediction,  with the same format of  data.est.
## case:  'regular'  for blockNNGP models with regular blocks, 
##    and 'irregular' for blockNNGP models with irregular blocks (Voronoy poligons).  
## n.blocks: number of  blocks  (regular or irregular)
## num.nb: number of  neighbor blocks. 
## Value: List with the following objects
##      - res2:  returns an object of class "inla". 
##      - myres: Posterior summary  under reparameterization.
##      - LPML, 
##      - MSE.Y.res (mean square estimation error of Y),
##      - MSP.Y.res (mean square prediction error of Y), 
##      - data1 (reordered data for estimation), 
##      - datap (reordered data for prediction),
##      - predictions: saida.pred

runblockNNGP = function(data.est,data.pred,case,n.blocks,num.nb){
  #  case='regular' or 'irregular'
  
  # name.cov="exponential" by default
  name.cov="matern"
  
  datafill <- blockNNGP.struct(case, 
                                data.est,
                               n.blocks, num.nb, 
                               name.cov,par.cov = 0.5,data.pred)
  
  # the data are reordered...
  data1   <- datafill[[1]]
  names(data1)
  
  ## data for estimation of two fields
  nloc    <- dim(data1)[1]
  Y       <- matrix(NA, nloc + nloc,  2)
  Y[1:nloc, 1]         <- data1$y1
  Y[nloc + 1:nloc, 2]  <- data1$y2
  
  Intercept.y1 <- c(rep(1,nloc), rep(NA, nloc))
  Intercept.y2 <- c(rep(NA, nloc), rep(1,nloc))
  
  spatial.field.y1 <- c(1:nloc, rep(NA, nloc))
  spatial.field.y2 <- c(rep(NA, nloc), 1:nloc)
  
  long1 <- c(data1$loc1, rep(NA, nloc))
  long2 <- c(rep(NA, nloc), data1$loc1)
  
  long1.2 <- c(data1$loc1^2, rep(NA, nloc))
  long2.2 <- c(rep(NA, nloc), data1$loc1^2)
  
  lat1 <- c(data1$x.3, rep(NA, nloc))
  lat2 <- c(rep(NA, nloc), data1$x.3)
  
  lat1.2 <- c(data1$x.3^2, rep(NA, nloc))
  lat2.2 <- c(rep(NA, nloc), data1$x.3^2)
  
  tmax1 <- c(data1$x.4, rep(NA, nloc))
  tmax2 <- c(rep(NA, nloc), data1$x.4)
  
  tmax1.2 <- c(data1$x.4^2, rep(NA, nloc))
  tmax2.2 <- c(rep(NA, nloc), data1$x.4^2)
  
  vehic11 <-  exp(data1$x.5)-1
  v1.sd   <-  (vehic11-mean(vehic11))/sd(vehic11)
  
  vehic1 <- c(v1.sd, rep(NA, nloc))
  vehic2 <- c(rep(NA, nloc), v1.sd)
  
  base.copy.y2 <-  c(rep(NA, nloc), 1:nloc)
  
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  idx   <- 1:nloc
  
  data2 <- list(y1=data1$y1, y2=data1$y2, loc1=data1$loc1, loc2 = data1$loc2,idx=idx,
               y = Y, Intercept.y1=Intercept.y1, 
               lat1=lat1, tmax1=tmax1, vehic1=vehic1, 
               long1=long1, long2=long2,
               spatial.field.y1=spatial.field.y1,
               Intercept.y2= Intercept.y2, 
               lat2=lat2, tmax2=tmax2, vehic2=vehic2, 
               spatial.field.y2=spatial.field.y2, 
               base.copy.y2=base.copy.y2)
  
  
  ## pc.prior ####
  
  ## Calculate hyperparameters
  # range <- (8*nu)/2*phis
  # set by user
  # P(rho<rho_0) = alpha1
  prior.range  <- c(60,0.95)  # P(rho<0.5*0.15) = 0.05
  # P(sigma>sigma0) = alpha2
  prior.sigma  <- c(100 ,0.05) # P(sigma>5) = 0.05
  
  ## set the blockNNGP as latent effect

    blockNNGP.model1  <- inla.blockNNGP.model(case, name.cov="matern",
                                          name.prior="pc.prior",datafill,
                                          prior.set=c(prior.range,prior.sigma))

    blockNNGP.model2  <- inla.blockNNGP.model(case, name.cov="matern",
                                          name.prior="pc.prior",datafill,
                                          prior.set=c(prior.range,prior.sigma))
  
  
  #names(data2)
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  hyper   <- list(beta = list(prior = 'normal', param = c(0, 10)))
  
  f.blockNNGP <- y ~ -1 + Intercept.y1 + Intercept.y2 + 
    long1 + lat1 + lat1.2 + tmax1 +  lat2 + vehic2 + 
    f(spatial.field.y1, model = blockNNGP.model1) +
    f(base.copy.y2, copy = "spatial.field.y1", fixed = FALSE, hyper = hyper)+  
  f(spatial.field.y2, model = blockNNGP.model2) 

  
  # inla function to fit the model
  # The user can change the family and priors here.
  res2  <- inla(f.blockNNGP, data = (data2),
               family = rep("poisson", 2),
               control.predictor =list(link=c(rep(1, nloc), rep(2, nloc))),
               control.compute = list(waic = TRUE, cpo = TRUE,config=TRUE))

  myres <- summary.blockNNGP(name.prior="pc.prior", resf=res2, data1,n.blocks, num.nb)
  LPML  <- round(sum(log(res2$cpo$cpo),na.rm = T),3)
  
  # estimation of y
  est.y1 <- res2$summary.fitted.values[,1][1:nloc]
  est.y2 <- res2$summary.fitted.values[,1][(nloc+1):(nloc*2)]
  par(mfrow=c(1,2))
  plot(data1$y1, est.y1)
  abline(0,1)
  plot(data1$y2, est.y2)
  abline(0,1)
  
  MSE.Y1.res <- round(mean((est.y1-data1$y1)^2),3); MSE.Y1.res
  MSE.Y2.res <- round(mean((est.y2-data1$y2)^2),3); MSE.Y2.res
  MSE.Y.res  <- (MSE.Y1.res+MSE.Y2.res)/2
    
  par(mfrow=c(1,2))
  plot_birds(resf=res2, data1,variable.est=log(est.y1) ,n.blocks=n.blocks, num.nb=n.blocks)
  plot_birds(resf=res2, data1,variable.est=log(est.y2),n.blocks=n.blocks, num.nb=num.nb )
  
  # spatial random effects
  est.w1 <- res2$summary.random$spatial.field.y1$mean
  est.w2 <- res2$summary.random$spatial.field.y2$mean
  par(mfrow=c(1,2))
  plot_birds(resf=res2, data1,variable.est=est.w1 ,n.blocks=n.blocks, num.nb=num.nb)
  plot_birds(resf=res2, data1,variable.est=est.w2,n.blocks=n.blocks, num.nb=num.nb )
  
  datap   <- datafill[[11]]
  names(datap)
  
  saida.pred <- blockNNGP_pred(case, n.blocks, num.nb,
                               data.est= data1, pred.data=datap,
                               res=res2,n.sample=1000,family="poisson")
  
  MSP.Y1.res <- round(mean((saida.pred$mean.y1-saida.pred$Y.pred)^2),3); MSP.Y1.res
  MSP.Y2.res <- round(mean((saida.pred$mean.y2-saida.pred$Y2.pred)^2),3); MSP.Y2.res
  MSP.Y.res  <- (MSP.Y1.res+ MSP.Y2.res)/2
    
  par(mfrow=c(1,2))
  plot(saida.pred$mean.y1,saida.pred$Y.pred)
  abline(0,1,col="red")
  plot(saida.pred$mean.y2,saida.pred$Y2.pred)
  abline(0,1,col="red")
  
  return(list(res2, myres, LPML, MSE.Y.res, MSP.Y.res, data1,datap, saida.pred))
}

## Description:
##
## runblockNNGP() function
## Arguments:
## data.est: data for estimation, data.est=data.frame(loc=loc, y1=y1,y2=y2, X=X)
##    loc: locations i.e. (long, lat), y1: response variable 1, 
##    y2: response variable 2, X: matrix of covariates. 
## data.pred: data for prediction,  with the same format of  data.est.
## case:  'NNGP'  for NNGP models 
## num.nb: number of neighbors. 
## Value: List with the following objects
##      - res2:  returns an object of class "inla". 
##      - myres: Posterior summary  under reparameterization.
##      - LPML, 
##      - MSE.Y.res (mean square estimation error of Y),
##      - MSP.Y.res (mean square prediction error of Y), 
##      - data1 (reordered data for estimation), 
##      - datap (reordered data for prediction),
##      - predictions: saida.pred

runNNGP = function(data.est,data.pred,case,num.nb){
 
  # name.cov="exponential" by default
  name.cov <- "matern"
  
  datafill <- blockNNGP.struct(case= 'NNGP', 
                               data.est, 
                               n.blocks=1, num.nb, name.cov,par.cov = 0.5,data.pred)
  
  #the data are reordered...
  data1   <- datafill[[1]]
  names(data1)
  
  ## data for estimation of two fields
  nloc    <- dim(data1)[1]
  Y       <- matrix(NA, nloc + nloc,  2)
  Y[1:nloc, 1]         <- data1$y1
  Y[nloc + 1:nloc, 2]  <- data1$y2
  
  Intercept.y1 <- c(rep(1,nloc), rep(NA, nloc))
  Intercept.y2 <- c(rep(NA, nloc), rep(1,nloc))
  
  spatial.field.y1 <- c(1:nloc, rep(NA, nloc))
  spatial.field.y2 <- c(rep(NA, nloc), 1:nloc)
  
  long1 <- c(data1$loc1, rep(NA, nloc))
  long2 <- c(rep(NA, nloc), data1$loc1)
  
  long1.2 <- c(data1$loc1^2, rep(NA, nloc))
  long2.2 <- c(rep(NA, nloc), data1$loc1^2)
  
  lat1 <- c(data1$x.3, rep(NA, nloc))
  lat2 <- c(rep(NA, nloc), data1$x.3)
  
  lat1.2 <- c(data1$x.3^2, rep(NA, nloc))
  lat2.2 <- c(rep(NA, nloc), data1$x.3^2)
  
  tmax1 <- c(data1$x.4, rep(NA, nloc))
  tmax2 <- c(rep(NA, nloc), data1$x.4)
  
  tmax1.2 <- c(data1$x.4^2, rep(NA, nloc))
  tmax2.2 <- c(rep(NA, nloc), data1$x.4^2)
  
  vehic11 <-  exp(data1$x.5)-1
  v1.sd   <- (vehic11-mean(vehic11))/sd(vehic11)
  
  vehic1 <- c(v1.sd, rep(NA, nloc))
  vehic2 <- c(rep(NA, nloc), v1.sd)
  
  base.copy.y2 <- c(rep(NA, nloc), 1:nloc)
  
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  idx   <- 1:nloc
  
  data2 <- list(y1=data1$y1, y2=data1$y2, loc1=data1$loc1, loc2 = data1$loc2,idx=idx,
               y = Y, Intercept.y1=Intercept.y1, 
               lat1=lat1, tmax1=tmax1, vehic1=vehic1, 
               long1=long1, long2=long2,
               spatial.field.y1=spatial.field.y1,
               Intercept.y2= Intercept.y2, 
               lat2=lat2, tmax2=tmax2, vehic2=vehic2, 
               spatial.field.y2=spatial.field.y2, 
               base.copy.y2=base.copy.y2)
  
  
  ## pc.prior ####
  
  ## Calculate hyperparameters
  #range <- (8*nu)/2*phis
  # set by user
  # P(rho<rho_0) = alpha1
  prior.range  <- c(60,0.95)  # P(rho<0.5*0.15) = 0.05
  # P(sigma>sigma0) = alpha2
  prior.sigma  <- c(100 ,0.05) # P(sigma>5) = 0.05
  
  ## set the NNGP as latent effect
  NNGP.model1 <- inla.blockNNGP.model(case= 'NNGP', name.cov="matern",
                                     name.prior="pc.prior",datafill,
                                     prior.set=c(prior.range,prior.sigma))
  NNGP.model2 <- inla.blockNNGP.model(case= 'NNGP', name.cov="matern",
                                     name.prior="pc.prior",datafill,
                                     prior.set=c(prior.range,prior.sigma))
  
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  hyper   <- list(beta = list(prior = 'normal', param = c(0, 10)))
  
  f.NNGP <- y ~ -1 + Intercept.y1 + Intercept.y2 + 
    long1 + lat1 + lat1.2 + tmax1 +  lat2 + vehic2 + 
    f(spatial.field.y1, model = NNGP.model1) +
    f(base.copy.y2, copy = "spatial.field.y1", fixed = FALSE, hyper = hyper)+  
    f(spatial.field.y2, model = NNGP.model2) 
  
  
  # inla function to fit the model
  # The user can change the family and priors here.
  res2 <- inla(f.NNGP, data = (data2),
               family = rep("poisson", 2),
               control.predictor=list(link=c(rep(1, nloc), rep(2, nloc))),
               control.compute = list(waic = TRUE, cpo = TRUE,config=TRUE))
  
  myres <- summary.blockNNGP(name.prior="pc.prior", resf=res2, data1,n.blocks=1, num.nb)
  LPML  <- round(sum(log(res2$cpo$cpo),na.rm = T),3)
  
  # estimation of y
  est.y1 <- res2$summary.fitted.values[,1][1:nloc]
  est.y2 <- res2$summary.fitted.values[,1][(nloc+1):(nloc*2)]
  par(mfrow=c(1,2))
  plot(data1$y1, est.y1)
  abline(0,1)
  plot(data1$y2, est.y2)
  abline(0,1)
  
  MSE.Y1.res <- round(mean((est.y1-data1$y1)^2),3); MSE.Y1.res
  MSE.Y2.res <- round(mean((est.y2-data1$y2)^2),3); MSE.Y2.res
  MSE.Y.res  <- (MSE.Y1.res+MSE.Y2.res)/2
  
  par(mfrow=c(1,2))
  plot_birds(resf=res2, data1,variable.est=log(est.y1) ,n.blocks=n.blocks, num.nb=n.blocks)
  plot_birds(resf=res2, data1,variable.est=log(est.y2),n.blocks=n.blocks, num.nb=num.nb )
  
  # spatial random effects
  est.w1 <- res2$summary.random$spatial.field.y1$mean
  est.w2 <- res2$summary.random$spatial.field.y2$mean
  par(mfrow=c(1,2))
  plot_birds(resf=res2, data1,variable.est=est.w1 ,n.blocks=n.blocks, num.nb=num.nb)
  plot_birds(resf=res2, data1,variable.est=est.w2,n.blocks=n.blocks, num.nb=num.nb )
  
  datap <- datafill[[6]]
  names(datap)
  
  saida.pred <- NNGP_pred(case='NNGP',  num.nb,
                          data.est= data1, pred.data=datap,
                          res=res2,n.sample=1000,family="poisson")
  
  MSP.Y1.res <- round(mean((saida.pred$mean.y1-saida.pred$Y.pred)^2),3); MSP.Y1.res
  MSP.Y2.res <- round(mean((saida.pred$mean.y2-saida.pred$Y2.pred)^2),3); MSP.Y2.res
  MSP.Y.res  <- (MSP.Y1.res+ MSP.Y2.res)/2
  
  par(mfrow=c(1,2))
  plot(saida.pred$mean.y1,saida.pred$Y.pred)
  abline(0,1,col="red")
  plot(saida.pred$mean.y2,saida.pred$Y2.pred)
  abline(0,1,col="red")
  
  return(list(res2, myres, LPML, MSE.Y.res, MSP.Y.res, data1,datap, saida.pred))
}


###########################
### model M1
#y1 <- data_robin$robin
#y2 <- data_robin$mourning
###########################

## load birds abundance data
dir.save = getwd()
load( paste("data.est.Rdata",sep=""))
load(paste("data.pred.Rdata",sep=""))

## run blockNNGP-GLMC models

case='regular'
n.partition <- 7
n.blocks    <- n.partition^2
num.nb 	    <- 2
resR1       <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resR1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resR1       <- NULL


case='regular'
n.partition <- 8
n.blocks    <- n.partition^2
num.nb 	    <- 2
resR1       <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resR1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resR1       <- NULL

case='regular'
n.partition <- 8
n.blocks    <- n.partition^2
num.nb 	    <- 4
resR1       <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resR1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resR1       <- NULL

case='irregular'
n.blocks  <- 40
num.nb 	  <- 2
resI1     <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resI1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resI1     <- NULL


case='irregular'
n.blocks  <- 40
num.nb 	  <- 4
resI1     <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resI1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resI1     <- NULL

case='irregular'
n.blocks  <- 53
num.nb 	  <- 2
resI1     <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resI1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resI1     <- NULL


case='irregular'
n.blocks  <- 53
num.nb 	  <- 4
resI1     <- runblockNNGP(data.est,data.pred,case,n.blocks,num.nb)
save(resI1, file=paste("M1",case,n.blocks,"-",num.nb,"res.Rdata",sep=''))
resI1     <- NULL


## run NNGP-GLMC models

case='NNGP'
num.nb 	<- 30
resN1   <- runNNGP(data.est,data.pred,case,num.nb)
save(resN1, file=paste("M1",case,num.nb,"res.Rdata",sep=''))
resN1   <- NULL

