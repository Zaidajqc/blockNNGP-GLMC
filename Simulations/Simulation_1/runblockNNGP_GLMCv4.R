############################################################################
## Authors: Z. Quiroz and C. Gonzáles
## Date: 6.05.2026
##
## Description:
##
##    The code performs a full Bayesian analysis of NNGP-GLM and blockNNGP-GLMC
##    models using Integrated Nested Laplace approximation (INLA). The BlockNNGP
##    latent effect is implemented using the rgeneric funcion.

## Arguments:
## case:  'simNNGP' for NNGP models, 
## 'irregular' for blockNNGP models with voronoi polygons. 
## formula:  NNGP  runs  ‘inla’ formula like
##  ‘y ~ 1 + x + f(idx, model = NNGP.model)’ 
## formula: blockNNGP_reg and blockNNGP  run  ‘inla’ formula like
##  ‘y ~ 1 + x + f(idx, model = blockNNGP.model)’ 
## family: A string indicating the likelihood family. For a list of possible
##          alternatives and use ‘inla.doc’ for detailed docs for
##          individual families.
## loc:  locations
## y: observed data
## X: covariates
## n.blocks: number of  blocks  (regular or irregular)
## num.nb: number of neighbors or neighbor blocks. 
##  Value:
##      ‘inla’ returns an object of class ‘"inla"’. 
##
#############################################################################


rm(list=ls())

getwd()

setwd("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/")

#######################
## functions required
####################### 
source("blockNNGPfunctions.R")
source("newQGLMCrgeneric.R")
source("NNGPrgeneric.R")
source("utils.R")
source("summary_functions.R")

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


sim_poisson.data = function(m,n){
  
  dir.save 	<- getwd()

  set.seed(m)
  
  loc 	<- cbind(runif(n,0,1), runif(n,0,1))
  
  ##Matern 1##
  sigma.sq <- 0.5
  phi 	   <- 12
  nu 	     <- 0.5
  sd.1     <- 0.01
  range    <- sqrt(8 * nu) / phi

  D 	 <- rdist(loc)
  C1 	 <- sigma.sq * exp(-phi * D) 
  Q1   <- solve(C1)
  
  k 	<- 2
  lambda12 <- 0.7
  M <-  diag(1, k)
  M[upper.tri(M)] <-  lambda12

  sigma.sq2 <-  1
  phi2      <- 6
  nu2       <- 0.5
  sd.2      <- 0.01
  
  range2    <- sqrt(8 * nu2) / phi2
  
  C2 	 <- sigma.sq2 * exp(-phi2*D) 
  

  gamma1<- matrix(M[1,], 1,k)
  
  Cnew1 <- kronecker(C1,t(gamma1)%*%(gamma1))
  
  gamma2 <- matrix(M[2,], 1,k)
  
  Cnew2 <- kronecker(C2,t(gamma2)%*%(gamma2))
  
  C <- Cnew1+ Cnew2
  
  
  cholC 	 	<- chol(C)
  nloc 	 <- dim(loc)[1]
  set.seed(m)
  rnorm_n.obs 	<- rnorm(k*nloc)
  w 		<- t(matrix(rnorm_n.obs, ncol=(k*nloc))%*%cholC )
  set.seed(m)
  X 		<- as.matrix(cbind(1, rnorm(nloc,40,1))) ## X = intercept + covariate
  B         <- as.matrix(c(0.2, 0.05))
  B2 		<- as.matrix(c(-0.1,0.08))
  xbeta <- c( X%*%B , X%*%B2)
  
  
  sq <- seq(1,(2*nloc),2)
  ind <- rep(0,2*nloc)
  k1 <- 1
  for(i in(sq)){
    ind[i]= k1
    ind[i+1]=nloc + k1
    k1=k1+1
  }
  xbnew <-  xbeta[ind]
  
  set.seed(m)
  y 		<- rpois(k*nloc, exp(xbnew + w)  ) ## y= X beta + w(spatial) + nugget

  orig.param <- t(data.frame(beta00 = B[1,1], beta01 = B[2,1], 
                             sigma21 = sigma.sq, phi1 = phi, sd.1 = sd.1, 
                             beta02 = B2[1,1], beta22 = B2[2,1], 
                             sigma22 = sigma.sq2, phi2 = phi2, lambda12 = lambda12, sd.2 = sd.2))
  

  
  indv1 = seq(1,k*nloc,by=2)
  ynew1 = y[indv1]
  ynew2 = y[indv1+1]
  par(mfrow=c(1,2))
  hist(ynew1)
  hist(ynew2)
  
  data.original <- data.frame(loc = loc, y1 = ynew1, y2 = ynew2, 
                              X = X, w1 = w[indv1],w2 = w[indv1+1])
  names(data.original)
  
  data.est <- data.original[1:(0.8*nloc),]
  data.pred <- data.original[(0.8*nloc+1):nloc,]
  
  
  return(list(data.est, data.pred, orig.param))
}




run_poissonIrreg0 = function(all.data, n.blocks, num.nb){
  
  data.est <- all.data[[1]]
  data.pred <- all.data[[2]] 
  
  ##########################################
  #### Run irregular voronoi block-NNGP models 
  ##########################################
  
  case <- 'irregular'
  
  # name.cov <- "exponential" by default
  name.cov <- "matern"
  
  datafill <- blockNNGP.struct(case = 'irregular', 
                               data.est, 
                               n.blocks, num.nb, name.cov, par.cov = 0.5, data.pred)
  
  # the data is reordered...
  data1 <- datafill[[1]]
  names(data1)
  
  nloc <- dim(data1)[1]
  
  sq <- seq(1,2*nloc,2)
  ind <- rep(0, 2*nloc)
  k2 <- 1
  for(i in(sq)){
    ind[i] <- k2
    ind[i+1] <- nloc + k2
    k2 <- k2+1
  }
  
  k <- 2
  
  yf <- c(data1$y1, data1$y2)
  ynew <-  yf[ind]
  wnew <- c(data1$w1, data1$w2)
  wnew <-  wnew[ind]
  blocks <- c(data1$blocks, data1$blocks)
  blocks <-  blocks[ind]
  loc1 <- c(data1$loc.1, data1$loc.1)
  loc1 <- loc1[ind]
  loc2 <- c(data1$loc.2, data1$loc.2)
  loc2 <- loc2[ind]
  
  indv1 <- seq(1, k*nloc,by = 2)
  xnew1 <-  c(rep(NA, nloc),rep(NA, nloc))
  xnew1[indv1] <- data1$X.2
  xnew2 <-  c(rep(NA, nloc),rep(NA, nloc))
  xnew2[indv1+1] <- data1$X.2
  int1 <-  c(rep(NA, nloc),rep(NA, nloc))
  int1[indv1] <- 1
  int2 <-  c(rep(NA, nloc),rep(NA, nloc))
  int2[indv1+1] <- 1
  
  k <- 2
  spatial.field.y <- 1:(k*nloc)
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  indg <- rep(1:2, nloc) 
  index <- rep(1:nloc, each= k) 
  data2<- list(y = ynew, x1=xnew1, x2 = xnew2,int1 = int1, int2 = int2,
               spatial.field.y = spatial.field.y, w = wnew, 
               blocks = blocks, index = index, indg = indg,
               loc.1 = loc1, loc.2 = loc2)
  
  ## Calculate hyperparameters
  #range <- (8*nu)/2*phi
  # set by user
  # P(rho<rho_0) <- alpha1
  prior.range  <- c( 0.5*0.15,0.05)  # P(rho<0.5*0.15) = 0.05
  # P(sigma>sigma0) <- alpha2
  prior.sigma  <- c(5 ,0.05) # P(sigma>5) = 0.05
  
  
  ## set the block-NNGP as latent effect
  #  name.prior="pc.prior"
  prior.set <- c(prior.range,prior.sigma)
  k <- 2
  prior.range <- c(prior.set[1:2])
  prior.sigma <- c(prior.set[3:4])
  lam1 <- -log(prior.range[2])*prior.range[1]
  lam2 <- -log(prior.sigma[2])/prior.sigma[1]
  
  initial.range <- log(prior.range[1]) + 1
  initial.sigma <- log(prior.sigma[1]) - 1
  W1 = datafill[[2]]
  
  maskWk <- function(k,W1){
    M <- diag(1, k)
    M [upper.tri(M )] <- 0.5 # just some value for W mask
    invM <- solve(M)
    listQnew <- NULL
    Wn <- W1
    gamma1inv <- matrix(invM[,1], k,1)
    Qnew1 <- kronecker(W1,(gamma1inv)%*%t(gamma1inv))
    
    gamma2inv <- matrix(invM[,2], k,1)
    Qnew2 <- kronecker(W1,(gamma2inv)%*%t(gamma2inv))
    Q <- Qnew1 +  Qnew2
    W <- Q*0
    mask_matrix <- Q !=0
    W[mask_matrix] <-1
    W <- as(W, "sparseMatrix")
    return(W)
  }
  
  W <- maskWk(k,W1)
  
  blockNNGP.model <- inla.rgeneric.define(inla.rgeneric.LMCk2_QblockNNGP.model.pc.prior, 
                                          k = k,    
                                          W = W,
                                          lam1 = lam1,
                                          lam2 = lam2,
                                          n = datafill[[3]], n.blocks = datafill[[4]],nb = datafill[[5]],
                                          ind_obs1 = datafill[[6]],num1 = datafill[[7]],indb = datafill[[8]],
                                          coords.D = datafill[[9]], a = prior.set[1],b = prior.set[2],
                                          initial.range = initial.range,
                                          initial.sigma = initial.sigma)
  
  f.blockNNGP <- y ~ 0+ int1 + int2 + x1+ x2 +
    f(spatial.field.y, model = blockNNGP.model) 
 
  # inla function to fit the model
  res4 <- inla(f.blockNNGP, data = (data2),
               family = "poisson",
               control.predictor= list(compute = TRUE),
               control.compute = list(waic = TRUE, cpo = TRUE,config=TRUE),
               verbose=TRUE)
  

  res <- summary.blockNNGP_LMC(name.prior="pc.prior", resf=res4, 
                               data1,n.blocks, num.nb, family="poisson")
  
  #  the data is reordered...
  datap <- datafill[[11]]
  names(datap)
  
  saida.pred <- blockNNGP_predLMC(case='irregular', n.blocks, num.nb,
                                  data.est= data2, pred.data=datap,
                                  res=res4,n.sample=1000,
                                  family="poisson")
  
  orig.param <- all.data[[3]]
  MSE.Y1.res <- round(mean((res4$summary.fitted.values[indv1,1]-data2$y[indv1])^2),3)
  MSE.Y2.res <- round(mean((res4$summary.fitted.values[(indv1+1),1]-data2$y[(indv1+1)])^2),3)
  MSP.Y1.res <- round(mean((saida.pred$mean.y1-saida.pred$Y.pred)^2),3)
  MSP.Y2.res <- round(mean((saida.pred$mean.y2-saida.pred$Y2.pred)^2),3)
  
  return(list(data.est=data2, data.pred=datap, res4, res, orig.param, saida.pred, MSP.Y1.res, MSP.Y2.res, MSE.Y1.res, MSE.Y2.res))  
  
  
}


run_poisson_NNGP0 = function(all.data,  num.nb){
  
  data.est <- all.data[[1]] 
  data.pred <- all.data[[2]] 
  
  ##########################################
  #### Run irregular voronoi block-NNGP models 
  ##########################################
  
  ############
  ## pc.prior
  ############
  
  
  case <- 'NNGP'
  

  # name.cov="exponential" by default
  name.cov <- "matern"
  
  datafill <- blockNNGP.struct(case = 'NNGP', data.est,
                               n.blocks = 1, num.nb, 
                               name.cov,par.cov = 0.5,data.pred)

  # the data is reordered...
  data1 <- datafill[[1]]
  names(datafill[[1]])
  
  nloc <- dim(data1)[1]
  
  
  sq <- seq(1,(2*nloc),2)
  ind <- rep(0,2*nloc)
  k2 <- 1
  for(i in(sq)){
    ind[i] <- k2
    ind[i+1] <- nloc + k2
    k2 <- k2+1
  }
  
  k <- 2
  
  yf <- c(data1$y1,data1$y2)
  ynew <-  yf[ind]
  wnew <- c(data1$w1, data1$w2)
  wnew <-  wnew[ind]
  blocks <- c(data1$blocks, data1$blocks)
  blocks <-  blocks[ind]
  loc1 <- c(data1$loc.1,data1$loc.1)
  loc1 <- loc1[ind]
  loc2 <- c(data1$loc.2,data1$loc.2)
  loc2 <- loc2[ind]
  
  indv1 <-seq(1,k*nloc,by=2)
  xnew1 <-  c(rep(NA, nloc),rep(NA, nloc))
  xnew1[indv1] <- data1$X.2
  xnew2 <-  c(rep(NA, nloc),rep(NA, nloc))
  xnew2[indv1+1] <- data1$X.2
  int1 <-  c(rep(NA, nloc),rep(NA, nloc))
  int1[indv1] <- 1
  int2 <-  c(rep(NA, nloc),rep(NA, nloc))
  int2[indv1+1] <- 1
  
  k <- 2
  spatial.field.y=1:(k*nloc)
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  
  indg <- rep(1:2, nloc) 
  index <- rep(1:nloc, each= k) 
  data2<- list(y = ynew, x1 = xnew1, x2 = xnew2, int1 = int1, int2 = int2,
               spatial.field.y = spatial.field.y, w = wnew, 
               blocks = blocks, index = index, indg = indg,
               loc.1 = loc1, loc.2 = loc2)
  
  ## Calculate hyperparameters
  # P(rho<rho_0) = alpha1
  prior.range  = c( 0.5*0.15,0.05)  # P(rho<0.5*0.15) = 0.05
  # P(sigma>sigma0) = alpha2
  prior.sigma  = c(5 ,0.05) # P(sigma>5) = 0.05
  
  
  ## set the block-NNGP as latent effect
  #  name.prior="pc.prior"
  prior.set <- c(prior.range,prior.sigma)
  k <- 2
  prior.range <- c(prior.set[1:2])
  prior.sigma <- c(prior.set[3:4])
  lam1 <- -log(prior.range[2])*prior.range[1]
  lam2 <- -log(prior.sigma[2])/prior.sigma[1]
  
  initial.range <- log(prior.range[1]) + 1
  initial.sigma <- log(prior.sigma[1]) - 1
  W1 <- datafill[[2]]
  
 
  maskWk <- function(k,W1){
    M <- diag(1, k)
    M [upper.tri(M )] <- 0.5 # just some value for W mask
    invM <- solve(M)
    listQnew <- NULL
    Wn <- W1
    gamma1inv <- matrix(invM[,1], k,1)
    Qnew1 <- kronecker(W1,(gamma1inv)%*%t(gamma1inv))
    
    gamma2inv <- matrix(invM[,2], k,1)
    Qnew2 <- kronecker(W1,(gamma2inv)%*%t(gamma2inv))
    Q <- Qnew1 +  Qnew2
    W <- Q*0
    mask_matrix <- Q !=0
    W[mask_matrix] <-1
    W <- as(W, "sparseMatrix")
    return(W)
  }
  
  W <- maskWk(k,W1)
  

  NNGP.model <- inla.rgeneric.define(inla.rgeneric.NNGP.model.pc.prior2,
                                     k = 2, W = W,
                                     lam1 = lam1, lam2 = lam2,
                                     n = nloc,
                                     coords.D = as.matrix(datafill[[3]], ncol,ncol), 
                                     AdjMatrix = as.matrix(datafill[[4]], ncol, ncol), 
                                     initial.range = initial.range,
                                     initial.sigma = initial.sigma)
  
  
  f.NNGP <- y ~ 0+ int1 + int2 + x1+ x2 + 
    f(spatial.field.y, model = NNGP.model)
  
  # inla function to fit the model
  
  res6 <- inla(f.NNGP, data = data2,
               family = "poisson",
                control.predictor = list(compute = TRUE),
               control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
               verbose = TRUE)
  
  res <- summary.blockNNGP_LMC(name.prior = "pc.prior", resf = res6, data2, 
                               n.blocks=1, num.nb, family = "poisson")
  
  # the data is reordered...
  datap <- datafill[[5]]
  names(datap)
  
  saida.pred <- NNGP_predLMC(case ='NNGP', num.nb,
                             data.est = data2, pred.data = datap,
                             res = res6, n.sample = 1000, family = "poisson")
  
  
  orig.param <- all.data[[3]]
  MSE.Y1.res <- round(mean((res6$summary.fitted.values[indv1,1]-data2$y[indv1])^2), 3)
  MSE.Y2.res <- round(mean((res6$summary.fitted.values[(indv1+1),1]-data2$y[(indv1+1)])^2), 3)
  MSP.Y1.res <- round(mean((saida.pred$mean.y1-saida.pred$Y.pred)^2), 3)
  MSP.Y2.res <- round(mean((saida.pred$mean.y2-saida.pred$Y2.pred)^2), 3)
  
  return(list(data.est=data2, data.pred=datap, res6, res, orig.param, saida.pred, MSP.Y1.res, MSP.Y2.res, MSE.Y1.res, MSE.Y2.res))  
  
}



# posterior results

post_results = function(resf){
  resinla <- resf[[3]]
  waic <- resinla$waic$waic
  time <- sum(resinla$cpu.used)
  LPML <- sum(log(resinla$cpo$cpo))
  
  
  MSE.Y1.res <- resf[[9]]
  MSE.Y2.res <- resf[[10]]
  MSP.Y1 <- resf[[7]]
  MSP.Y2 <- resf[[8]]
  
  orig.param <- resf[[5]]
  orig.param<- orig.param[c(1:4,6:10)] #for poisson
  
  hyperpar <- resf[[4]]$hyperpar
  fixed <- resf[[4]]$fixed
  lambda <- resf[[4]]$hyperparLMC
  
  est.beta00 <-fixed[1,1]
  est.beta01 <-fixed[3,1]
  est.beta02 <-fixed[2,1]
  est.beta22 <-fixed[4,1]
  est.sigma21 <-hyperpar[1,1]
  est.phi1 <-hyperpar[2,1]
  est.sigma22 <-hyperpar[3,1]
  est.phi2 <-hyperpar[4,1]
  est.lambda <- lambda[,1]
  
  covbeta00=0; if(orig.param[1]>= fixed[1,3] & orig.param[1]<= fixed[1,5]) covbeta00 =1
  covbeta01=0; if(orig.param[2]>= fixed[3,3] & orig.param[2]<= fixed[3,5]) covbeta01 =1
  covsigma21=0; if(orig.param[3]>= hyperpar[1,3] & orig.param[3]<= hyperpar[1,5]) covsigma21 =1
  covphi1=0; if(orig.param[4]>=  hyperpar[2,3] & orig.param[4]<=  hyperpar[2,5]) covphi1 =1
  
  covbeta02=0; if(orig.param[5]>= fixed[2,3] & orig.param[5]<= fixed[2,5]) covbeta02 =1
  covbeta22=0; if(orig.param[6]>= fixed[4,3] & orig.param[6]<= fixed[4,5]) covbeta22 =1
  covsigma22=0; if(orig.param[7]>= hyperpar[3,3] & orig.param[7]<= hyperpar[3,5]) covsigma22 =1
  covphi2=0; if(orig.param[8]>=  hyperpar[4,3] & orig.param[8]<=  hyperpar[4,5]) covphi2 =1
  
  covlambda=0; if(orig.param[9]>=  lambda[,3] & orig.param[9]<=  lambda[,5]) covlambda =1
  
  names.mean.estparam <- c("est.beta00", "est.beta01", "est.beta02","est.beta22", "est.sigma21", "est.phi1", "est.sigma22", "est.phi2", "est.lambda")
  mean.estparam <- round(c(est.beta00, est.beta01, est.beta02, est.beta22,est.sigma21, est.phi1, est.sigma22, est.phi2, est.lambda),3)
  names.cov.estparam <- c("covbeta00", "covbeta01", "covbeta02","covbeta22","covsigma21", "covphi1", "covsigma22", "covphi2", "covlambda")
  cov.param <- c(covbeta00, covbeta01, covbeta02, covbeta22,covsigma21, covphi1, covsigma22, covphi2,covlambda)
  
  
  name.crit <- c("waic", "LPML","time(sec)","MSP.Y1","MSP.Y2","MSE.Y1.res", "MSE.Y2.res") #, "RMSPy")
  crit.model<- round(c(  waic, LPML,  time, MSP.Y1, MSP.Y2,MSE.Y1.res,MSE.Y2.res),3) #, RMSP
  name.full.results <- c(t(names.mean.estparam), t(names.cov.estparam), t(name.crit))
  full.results <- c(t(mean.estparam), t(cov.param), t(crit.model) )
  final.res <- data.frame(res=name.full.results, val=full.results)
  #save(res=final.res, file="multilevel_siblockNNGP.Rdata")
  
  return(final.res)
  
}



dir.save = getwd()

########################################
#### Run regular block-NNGP models 
#######################################

m = 2
n = 2500
all.data <- sim_poisson.data(m,n)
hist(all.data[[1]]$y1)
hist(all.data[[1]]$y2)


##########################################
#### Run irregular voronoi block-NNGP models 
##########################################

## Irregular 64-2 ####
res5 <- run_poissonIrreg0(all.data, n.blocks=30, num.nb=2)
full.results5 <- post_results(res5)
save(res5=res5, file="simblockNNGP_res5_30_2.Rdata")
## Irregular 64-2 ####
res6 <- run_poissonIrreg0(all.data, n.blocks=30, num.nb=4)
full.results6 <- post_results(res6)
save(res6=res6, file="simblockNNGP_res6_30_4.Rdata")
## Irregular 64-2 ####
res7 <- run_poissonIrreg0(all.data, n.blocks=40, num.nb=2)
full.results7 <- post_results(res7)
save(res7=res7, file="simblockNNGP_res7_40_2.Rdata")
## Irregular 64-2 ####
res8 <- run_poissonIrreg0(all.data, n.blocks=50, num.nb=2)
full.results8 <- post_results(res8)
save(res8=res8, file="simblockNNGP_res8_50_2.Rdata")

res8 <- run_poissonIrreg0(all.data, n.blocks=40, num.nb=4)
full.results8 <- post_results(res8)
save(res8=res8, file="simblockNNGP_res8_40_4.Rdata")

res8 <- run_poissonIrreg0(all.data, n.blocks=50, num.nb=4)
full.results8 <- post_results(res8)
save(res8=res8, file="simblockNNGP_res8_50_4.Rdata")

##########################################
#### Run NNGP models 
##########################################

res9 <- run_poisson_NNGP0(all.data,  num.nb=10)
full.results9 <- post_results(res9)
save(res9=res9, file="simNNGP_res9_10.Rdata")

res10 <- run_poisson_NNGP0(all.data,  num.nb=20)
full.results10 <- post_results(res10)
save(res10=res10, file="simNNGP_res10_20.Rdata")

res11 <- run_poisson_NNGP0(all.data,  num.nb=30)
full.results11 <- post_results(res11)
save(res11=res11, file="simNNGP_res11_30.Rdata")


##################################
######## resultados #############
#################################

#load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res1_64_2.Rdata")
#load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res2_64_4.Rdata")
#load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res4_100_6.Rdata")

load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res5_30_2.Rdata")
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res6_30_4.Rdata")
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res8_40_4.Rdata")

load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simNNGP_res9_10.Rdata")
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simNNGP_res10_20.Rdata")
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simNNGP_res11_30.Rdata")

setwd("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/")

returnMarg <- function(res){
  marg.sigmasq1 <- inla.tmarginal(function(x) { (exp(x))^2},
                                  res$marginals.hyperpar[[1]])
  marg.phi1 <- inla.tmarginal(function(x) { 2/exp(x)},
                              res$marginals.hyperpar[[2]])
  marg.sigmasq2 <- inla.tmarginal(function(x) { (exp(x))^2},
                                  res$marginals.hyperpar[[3]])
  marg.phi2 <- inla.tmarginal(function(x) { 2/exp(x)},
                              res$marginals.hyperpar[[4]])
  marg.lambda <- res$marginals.hyperpar[[5]]
  
  all.marg <- list(marg.phi1, marg.sigmasq1,
                   marg.phi2, marg.sigmasq2,marg.lambda)
return(all.marg)  
}


res30.2I <- NULL
res30.2I$all.marg <- returnMarg(res5[[3]])
res30.4I <- NULL
res30.4I$all.marg <- returnMarg(res6[[3]])
res40.4I <- NULL
res40.4I$all.marg <- returnMarg(res8[[3]])

resNNGP_10 <- NULL
resNNGP_10$all.marg <- returnMarg(res9[[3]])
resNNGP_20 <- NULL
resNNGP_20$all.marg <- returnMarg(res10[[3]])
resNNGP_30 <- NULL
resNNGP_30$all.marg <- returnMarg(res11[[3]])


library(gridExtra)


## marginal for irregular

df <-data.frame(phi1.orig = 12,  phi2.orig = 6,
                phi1= c(res30.2I$all.marg[[1]][,1],res30.4I$all.marg[[1]][,1],res40.4I$all.marg[[1]][,1]),
                marg.phi1= c(res30.2I$all.marg[[1]][,2],res30.4I$all.marg[[1]][,2],res40.4I$all.marg[[1]][,2]),
                phi2= c(res30.2I$all.marg[[3]][,1],res30.4I$all.marg[[3]][,1],res40.4I$all.marg[[3]][,1]),
                marg.phi2= c(res30.2I$all.marg[[3]][,2],res30.4I$all.marg[[3]][,2],res40.4I$all.marg[[3]][,2]),
                model =c(rep("30-2", dim(res30.2I$all.marg[[1]])[1]), rep("30-4", dim(res30.4I$all.marg[[1]])[1]),  rep("40-4", dim(res40.4I$all.marg[[1]])[1])))


plot21 <- ggplot(df) + 
  geom_line(aes(x=phi1, y = marg.phi1 , color = model))+
  geom_vline(data = df,
             aes(xintercept = phi1.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(phi[1]))

plot22 <- ggplot(df) + 
  geom_line(aes(x=phi2, y = marg.phi2 , color = model))+
  geom_vline(data = df,
             aes(xintercept = phi2.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(phi[2]))


df1 <- data.frame(sigma1.orig = 0.5,  sigma2.orig = 1,
                  sigma1 = c(res30.2I$all.marg[[2]][,1],res30.4I$all.marg[[2]][,1], res40.4I$all.marg[[2]][,1]),
                  marg.sigma1 = c(res30.2I$all.marg[[2]][,2],res30.4I$all.marg[[2]][,2],res40.4I$all.marg[[2]][,2]),
                  sigma2 = c(res30.2I$all.marg[[4]][,1],res30.4I$all.marg[[4]][,1],res40.4I$all.marg[[4]][,1]),
                  marg.sigma2 = c(res30.2I$all.marg[[4]][,2],res30.4I$all.marg[[4]][,2],res40.4I$all.marg[[4]][,2]),
                  model = c(rep("30-2", dim(res30.2I$all.marg[[1]])[1]), rep("30-4", dim(res30.4I$all.marg[[1]])[1]), rep("40-4", dim(res40.4I$all.marg[[1]])[1])))

plot23 <- ggplot(df1) + 
  geom_line(aes(x=sigma1, y = marg.sigma1 , color = model))+
  geom_vline(data = df1,
             aes(xintercept = sigma1.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(sigma[1]^2))

plot24 <- ggplot(df1) + 
  geom_line(aes(x=sigma2, y = marg.sigma2 , color = model))+
  geom_vline(data = df1,
             aes(xintercept = sigma2.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y = "marginal", x = bquote(sigma[2]^2))



df2 <- data.frame(lambda1.orig = 0.7,  
                  lambda1= c(res30.2I$all.marg[[5]][,1],res30.4I$all.marg[[5]][,1], res40.4I$all.marg[[5]][,1]),
                  marg.lambda1= c(res30.2I$all.marg[[5]][,2],res30.4I$all.marg[[5]][,2],res40.4I$all.marg[[5]][,2]),
                  model =c(rep("30-2", dim(res30.2I$all.marg[[5]])[1]), rep("30-4", dim(res30.4I$all.marg[[5]])[1]), rep("40-4", dim(res40.4I$all.marg[[5]])[1])))


plot25 <- ggplot(df2) + 
  geom_line(aes(x=lambda1, y = marg.lambda1 , color = model))+
  geom_vline(data = df2,
             aes(xintercept = lambda1.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(lambda[12]))


## marginal for nngp

df <-  data.frame(phi1.orig = 12,  phi2.orig = 6,
                  phi1 = c(resNNGP_10$all.marg[[1]][,1],resNNGP_20$all.marg[[1]][,1],resNNGP_30$all.marg[[1]][,1]),
                  marg.phi1 = c(resNNGP_10$all.marg[[1]][,2],resNNGP_20$all.marg[[1]][,2],resNNGP_30$all.marg[[1]][,2]),
                  phi2 = c(resNNGP_10$all.marg[[3]][,1],resNNGP_20$all.marg[[3]][,1],resNNGP_30$all.marg[[3]][,1]),
                  marg.phi2 = c(resNNGP_10$all.marg[[3]][,2],resNNGP_20$all.marg[[3]][,2],resNNGP_30$all.marg[[3]][,2]),
                  model = c(rep("10", dim(resNNGP_10$all.marg[[1]])[1]), rep("20", dim(resNNGP_20$all.marg[[1]])[1]), rep("30", dim(resNNGP_30$all.marg[[1]])[1])) )



plot31 <- ggplot(df) + 
  geom_line(aes(x=phi1, y = marg.phi1 , color = model))+
  geom_vline(data = df,
             aes(xintercept = phi1.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(phi[1]))

plot32 <- ggplot(df) + 
  geom_line(aes(x=phi2, y = marg.phi2 , color = model))+
  geom_vline(data = df,
             aes(xintercept = phi2.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(phi[2]))



df1 <- data.frame(sigma1.orig = 0.5,  sigma2.orig = 1,
                  sigma1 = c(resNNGP_10$all.marg[[2]][,1],resNNGP_20$all.marg[[2]][,1],resNNGP_30$all.marg[[2]][,1]),
                  marg.sigma1 = c(resNNGP_10$all.marg[[2]][,2],resNNGP_20$all.marg[[2]][,2],resNNGP_30$all.marg[[2]][,2]),
                  sigma2 = c(resNNGP_10$all.marg[[4]][,1],resNNGP_20$all.marg[[4]][,1],resNNGP_30$all.marg[[4]][,1]),
                  marg.sigma2 = c(resNNGP_10$all.marg[[4]][,2],resNNGP_20$all.marg[[4]][,2],resNNGP_30$all.marg[[4]][,2]),
                  model = c(rep("10", dim(resNNGP_10$all.marg[[1]])[1]), rep("20", dim(resNNGP_20$all.marg[[1]])[1]), rep("30", dim(resNNGP_30$all.marg[[1]])[1])))


plot33 <- ggplot(df1) + 
  geom_line(aes(x = sigma1, y = marg.sigma1 , color = model))+
  geom_vline(data = df1,
             aes(xintercept = sigma1.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(sigma[1]^2))

plot34 <- ggplot(df1) + 
  geom_line(aes(x = sigma2, y = marg.sigma2 , color = model))+
  geom_vline(data = df1,
             aes(xintercept = sigma2.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(sigma[2]^2))



df2 <- data.frame(lambda1.orig = 0.7,  
                  lambda1 = c(resNNGP_10$all.marg[[5]][,1],resNNGP_20$all.marg[[5]][,1],resNNGP_30$all.marg[[5]][,1]),
                  marg.lambda1 = c(resNNGP_10$all.marg[[5]][,2],resNNGP_20$all.marg[[5]][,2],resNNGP_30$all.marg[[5]][,2]),
                  model =c(rep("10", dim(resNNGP_10$all.marg[[5]])[1]), rep("20", dim(resNNGP_20$all.marg[[5]])[1]), rep("30", dim(resNNGP_30$all.marg[[5]])[1])))


plot35 <- ggplot(df2) + 
  geom_line(aes(x = lambda1, y = marg.lambda1 , color = model))+
  geom_vline(data = df2,
             aes(xintercept = lambda1.orig), linewidth = 1, color="gray50",linetype = "dashed")+
  labs(y= "marginal", x = bquote(lambda[12]))



plot1 <- grid.arrange(plot31, plot21, plot32, plot22,
                      plot33, plot23, plot34, plot24,
                      plot35, plot25, nrow = 5)
plot1
library(ggplot2)
ggsave("plot1.eps", plot = plot1, device = "eps") 


#############################
## plot for best models..
############################
library(gridExtra)
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim1/simblockNNGP_res6_30_4.Rdata")

res1 <- res6

#call data.est
data.est <- res1[[1]]
#call data.pred
data.pred <- res1[[2]]
# call res from inla
res <- res1[[3]]

k     <- 2
nloc  <- length(data.est$y)/2
indv1 <- seq(1,k*nloc,by=2)

M <- 30
num.nb <- 4
v1  <- data.est$w[indv1]
v2  <-  data.est$w[indv1] * res$summary.hyperpar[5,1] + data.est$w[indv1+1]
rz1 <-res$summary.random$spatial.field.y$mean[indv1]
rz2 <-res$summary.random$spatial.field.y$mean[indv1+1] 
library(patchwork)
plot3 <- data.frame(v1,rz1,v2,rz2)
v1.plot <- ggplot(data = plot3,mapping = aes(x=v1,y=rz1))+
  geom_point()+ geom_abline(colour = "firebrick")+
  labs(x=expression(v[1]),y=expression(E(v[1]*"|"*y)))+
  ggtitle(paste("M =", M, ", nb =", num.nb))
v2.plot <- ggplot(data = plot3,mapping = aes(x=v2,y=rz2))+
  geom_point()+ geom_abline(colour = "firebrick")+
  labs(x=expression(v[2]),y=expression(E(v[2]*"|"*y)))+
  ggtitle(paste("M =", M, ", nb =", num.nb))

grid.arrange(v1.plot, v2.plot, nrow = 2)


#plot 4
y1 <- data.est$y[indv1]
y2 <- data.est$y[indv1+1]
esty1 <- res$summary.fitted.values$mean[indv1]
esty2 <- res$summary.fitted.values$mean[indv1+1]

plot4 <- data.frame( y1= y1, y2= y2,
                     y1est=esty1,y2est=esty2)
# buscar para poner los dos gráficos en uno
y1.plot <- ggplot(data = plot4,mapping = aes(x=y1,y=y1est))+
  geom_point()+ geom_abline(colour = "firebrick")+
  labs(x=expression(y[1]),y=expression(E(y[1]*"|"*y)))+
  ggtitle(paste("M =", M, ", nb =", num.nb))+
  theme_bw(base_size=16)
y2.plot <- ggplot(data = plot4,mapping = aes(x=y2,y=y2est))+
  geom_point()+ geom_abline(colour = "firebrick")+
  labs(x=expression(y[2]),y=expression(E(y[2]*"|"*y)))+
  ggtitle(paste("M =", M, ", nb =", num.nb))+
  theme_bw(base_size=16)


grid.arrange(y1.plot, y2.plot, nrow = 2)

#plot 5
saida.pred <- res1[[6]]
y1p <- saida.pred$Y.pred
y2p <- saida.pred$Y2.pred

predy1 <- saida.pred$mean.y1
predy2 <- saida.pred$mean.y2

cor(y1p, predy1)
cor(y2p, predy2)

plot5 <- data.frame( y1p = y1p, y2p = y2p,
                     predy1 = predy1, predy2 = predy2)

y1.plotP <- ggplot(data = plot5,mapping = aes(x=y1p, y=predy1))+
  geom_point()+ geom_abline(colour = "firebrick")+
  labs(x=expression(y[1]),y=expression(E(y[10]*"|"*y)))+
  ggtitle(paste("M =", M, ", nb =", num.nb)) +
  theme_bw(base_size=16)
y2.plotP <- ggplot(data = plot5,mapping = aes(x=y2p, y=predy2))+
  geom_point()+ geom_abline(colour = "firebrick")+
  labs(x=expression(y[2]),y=expression(E(y[20]*"|"*y)))+
  ggtitle(paste("M =", M, ", nb =", num.nb))+ 
  theme_bw(base_size=16)

   
library(ggpubr)
plot2 <- ggarrange(y1.plot, y2.plot,
          y1.plotP, y2.plotP,nrow = 2, ncol=2,
          labels = c("a)", "b)", "c)", "d)"),
          font.label = list(size = 16, color = "black", face = "bold", family = NULL)) 
plot2          
          
ggsave("plot2.eps", plot = plot2, device = "eps") 


# ## posterior mean of spatial random  effects v
loc <- cbind(data.est$loc.1, data.est$loc.2)[indv1,]
w1  <- data.est$w[indv1]
w2  <- data.est$w[indv1+1]
lambda <- res$summary.hyperpar[5,1]
v1  <- w1
v2  <- w1*lambda + w2
rz1 <- res$summary.random$spatial.field.y$mean[indv1]
rz2 <- res$summary.random$spatial.field.y$mean[indv1+1] 


int.elev1 <- mba.surf(cbind(loc,v1), 100, 100, extend=TRUE)$xyz.est
int.elev2 <- mba.surf(cbind(loc,v2), 100, 100, extend=TRUE)$xyz.est
int.elev12 <- mba.surf(cbind(loc,rz1), 100, 100, extend=TRUE)$xyz.est
int.elev22 <- mba.surf(cbind(loc,rz2), 100, 100, extend=TRUE)$xyz.est

grid <- expand.grid(x = int.elev1$x, y = int.elev1$y)

pred_dfv1 <- data.frame(x1=grid[,1],x2=grid[,2],v1= c(int.elev1$z))
pred_dfv2 <- data.frame(x1=grid[,1],x2=grid[,2],v2=c(int.elev2$z))
pred_dfv1est <- data.frame(x1=grid[,1],x2=grid[,2],v1=c(int.elev12$z))
pred_dfv2est <- data.frame(x1=grid[,1],x2=grid[,2],v2=c(int.elev22$z))
#pred_dfw3est <- data.frame(x1=grid[,1],x2=grid[,2],w3=c(int.elev23$z))

mapv1 <- ggplot(pred_dfv1, aes(x = x1, y = x2, fill = v1)) +
geom_raster() +
  scale_fill_gradientn(colours=tim.colors(100)) + 
labs(x="Longitude",y="Latitude",fill=expression(v[1]))+
  theme_bw(base_size=16)
  
mapv2 <- ggplot(pred_dfv2, aes(x = x1, y = x2, fill = v2)) +
  geom_raster() +
  scale_fill_gradientn(colours=tim.colors(100))+ 
  labs(x="Longitude",y="Latitude",fill=expression(v[2]))+
  theme_bw(base_size=16)

mapv1est <- ggplot(pred_dfv1est, aes(x = x1, y = x2, fill = v1)) +
  geom_raster() +
  scale_fill_gradientn(colours=tim.colors(100))+ 
  labs(x="Longitude",y="Latitude",fill=expression(E(v[1]*"|"*y)))+
  theme_bw(base_size=16)

mapv2est <- ggplot(pred_dfv2est, aes(x = x1, y = x2, fill = v2)) +
  geom_raster() +
  scale_fill_gradientn(colours=tim.colors(100))+ 
  labs(x="Longitude",y="Latitude",fill=expression(E(v[2]*"|"*y)))+
  theme_bw(base_size=16)

plot3 <- ggarrange(mapv1, mapv2,
                    mapv1est, mapv2est,
                   nrow = 2, ncol=2,
                   labels = c("a)", "b)", "c)", "d)"),
                   font.label = list(size = 16, color = "black", face = "bold", family = NULL)) 
plot3

ggsave("plot3.eps", plot = plot3, device = "eps") 

