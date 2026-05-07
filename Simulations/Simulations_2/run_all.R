############################################################################
## Authors: Z. Quiroz, C. Pizango
## Date: 6.05.2026
##
## Description:
##
##    The code performs a full Bayesian analysis of NNGP-GLMC and blockNNGP-GLMC
##    models using Integrated Nested Laplace approximation (INLA). The BlockNNGP
##   latent effect is implemented using the rgeneric funcion.

## Arguments:
## case:  'simNNGP' for NNGP models, 'regular' for blockNNGP models with regular blocks, 
## and 'irregular' for blockNNGP models with irregular blocks. 
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

setwd("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/")

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
library(meshed)

##blockNNGP


sim_poisson.data2 = function(m, n){
  
  dir.save=getwd()

  set.seed(m)
  
  loc 	<- cbind(runif(n, 0, 1), runif(n, 0, 1))
  
  ##Matern 1##
  sigma.sq <- 0.5
  phi 	   <- 12
  nu 	     <- 0.5
  sd.1     <- 0.01
  
  range    <- sqrt(8 * nu) / phi
  D 	     <- rdist(loc)
  C1 	     <- sigma.sq * exp(-phi * D) 

  k <- 2
  lambda12 <- 0.7
  M <-  diag(1, k)
  M[upper.tri(M)] <-  lambda12

  sigma.sq2 <-  1
  phi2      <- 6
  nu2       <- 0.5
  sd.2      <- 0.01
  
  range2    <- sqrt(8 * nu2) / phi2
  
  C2 	 <- sigma.sq2 * exp(-phi2*D) 
  
  
  # C original
  gamma1 <- matrix(M[1,], 1,k)
  
  Cnew1  <- kronecker(C1,t(gamma1)%*%(gamma1))
  
  gamma2 <- matrix(M[2,], 1,k)
  
  Cnew2  <- kronecker(C2,t(gamma2)%*%(gamma2))
  
  C <- Cnew1+ Cnew2
 
  
  cholC 	 	<- chol(C)
  nloc 	 <- dim(loc)[1]
  set.seed(m)
  rnorm_n.obs 	<- rnorm(k*nloc)
  w 		<- t(matrix(rnorm_n.obs, ncol = (k*nloc))%*%cholC )
  set.seed(m)
  X 		<- as.matrix(cbind(1, rnorm(nloc,40,1))) 
  B     <- as.matrix(c(0.2, 0.05))
  B2 		<- as.matrix(c(-0.1,0.08))
  xbeta <- c( X%*%B , X%*%B2)
  
  
  sq  <- seq(1, 2*nloc, 2)
  ind <- rep(0, 2*nloc)
  k1  <- 1
  for(i in(sq)){
    ind[i]   <- k1
    ind[i+1] <-nloc + k1
    k1       <- k1+1
  }
  xbnew      <-  xbeta[ind]
  
  set.seed(m)
  y 		<- rpois(k*nloc, exp(xbnew + w)  ) 

  orig.param <- t(data.frame(beta00=B[1,1], beta01= B[2,1], 
                             sigma21=sigma.sq, phi1=phi, sd.1=sd.1, 
                             beta02=B2[1,1], beta22=B2[2,1], 
                             sigma22=sigma.sq2,phi2=phi2,lambda12=lambda12, sd.2=sd.2))
  
  indv1 <- seq(1, k*nloc, by=2)
  ynew1 <- y[indv1]
  ynew2 <- y[indv1+1]
  par(mfrow = c(1,2))
  hist(ynew1)
  hist(ynew2)
  
  data.original <- data.frame(loc = loc, y1 = ynew1, y2 = ynew2, 
                              X = X, w1 = w[indv1], w2 = w[indv1+1])
  names(data.original)
  
  data.est  <- data.original[1:(0.8*nloc),]
  data.pred <- data.original[(0.8*nloc+1):nloc,]
  
  
  return(list(data.est, data.pred, orig.param))
}



run_poissonIrreg0 = function(all.data, n.blocks, num.nb){

  data.est <- all.data[[1]]
  data.pred <- all.data[[2]] 
  
  ##########################################
  #### Run irregular voronoi block-NNGP models 
  ##########################################
  
  case='irregular'
  
  # name.cov="exponential" by default
  name.cov="matern"
  
  datafill <- blockNNGP.struct(case= 'irregular', 
                               data.est, 
                               n.blocks, num.nb, name.cov,par.cov = 0.5,data.pred)
  
  #  the data is reordered...
  data1 <- datafill[[1]]
  names(data1)
  
  nloc <- dim(data1)[1]
  
  sq  <- seq(1,(2*nloc),2)
  ind <- rep(0,2*nloc)
  k2  <- 1
  for(i in(sq)){
    ind[i]   <- k2
    ind[i+1] <- nloc + k2
    k2       <- k2+1
  }
  
  k <-2
  
  yf   <- c(data1$y1, data1$y2)
  ynew <-  yf[ind]
  wnew <- c(data1$w1, data1$w2)
  wnew <-  wnew[ind]
  blocks <- c(data1$blocks, data1$blocks)
  blocks <-  blocks[ind]
  loc1 <- c(data1$loc.1, data1$loc.1)
  loc1 <- loc1[ind]
  loc2 <- c(data1$loc.2, data1$loc.2)
  loc2 <- loc2[ind]
  
  indv1 <- seq(1, k*nloc, by=2)
  xnew1 <-  c(rep(NA, nloc), rep(NA, nloc))
  xnew1[indv1] <- data1$X.2
  xnew2 <-  c(rep(NA, nloc), rep(NA, nloc))
  xnew2[indv1+1] <- data1$X.2
  int1 <-  c(rep(NA, nloc), rep(NA, nloc))
  int1[indv1] <- 1
  int2 <-  c(rep(NA, nloc), rep(NA, nloc))
  int2[indv1+1] <- 1
  
  k <- 2
  spatial.field.y = 1:(k*nloc)
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  indg  <- rep(1:2, nloc) 
  index <- rep(1:nloc, each = k) 
  data2 <- list(y = ynew, x1 = xnew1, x2 = xnew2, int1 = int1, int2 = int2,
               spatial.field.y = spatial.field.y, w = wnew, 
               blocks = blocks, index = index, indg = indg,
               loc.1 = loc1, loc.2 = loc2)
  
  # set by user
  # P(rho<rho_0) = alpha1
  prior.range  = c( 0.5*0.15,0.05)  # P(rho<0.5*0.15) = 0.05
  # P(sigma>sigma0) = alpha2
  prior.sigma  = c(5 ,0.05) # P(sigma>5) = 0.05
  
  
  ## set the block-NNGP as latent effect
  prior.set=c(prior.range,prior.sigma)
  k=2
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
  
  blockNNGP.model <- inla.rgeneric.define(inla.rgeneric.LMCk2_QblockNNGP.model.pc.prior, 
                                          k = k,    
                                          W = W,
                                          lam1 = lam1,
                                          lam2 = lam2,
                                          n = datafill[[3]], n.blocks = datafill[[4]], nb = datafill[[5]],
                                          ind_obs1 = datafill[[6]], num1 = datafill[[7]], indb = datafill[[8]],
                                          coords.D = datafill[[9]], a = prior.set[1],b = prior.set[2],
                                          initial.range = initial.range,
                                          initial.sigma = initial.sigma)
  
  f.blockNNGP <- y ~ 0+ int1 + int2 + x1+ x2 +
    f(spatial.field.y, model = blockNNGP.model) 

  # inla function to fit the model
  res4 <- inla(f.blockNNGP, data = (data2),
               family = "poisson",
               control.predictor= list(compute = TRUE),
               control.compute = list(waic = TRUE, cpo = TRUE, config=TRUE),
               verbose=TRUE)
  

  res <- summary.blockNNGP_LMC(name.prior = "pc.prior", resf = res4, 
                               data1, n.blocks, num.nb, family = "poisson")
  
  #  the data is reordered...
  datap <- datafill[[11]]
  names(datap)

  saida.pred <- blockNNGP_predLMC(case = 'irregular', n.blocks, num.nb,
                                  data.est = data2, pred.data = datap,
                                  res = res4,n.sample = 1000,
                                  family = "poisson")
  
  orig.param <- all.data[[3]]
  MSE.Y1.res <- round(mean((res4$summary.fitted.values[indv1,1]-data2$y[indv1])^2),3)
  MSE.Y2.res <- round(mean((res4$summary.fitted.values[(indv1+1),1]-data2$y[(indv1+1)])^2),3)
  MSP.Y1.res <- round(mean((saida.pred$mean.y1-saida.pred$Y.pred)^2),3)
  MSP.Y2.res <- round(mean((saida.pred$mean.y2-saida.pred$Y2.pred)^2),3)
  

  return(list(data.est = data2, data.pred = datap, res4, res, orig.param, 
              saida.pred, MSP.Y1.res, MSP.Y2.res, MSE.Y1.res, MSE.Y2.res))  

}




run_poisson_NNGP0 = function(all.data,  num.nb){
  
  data.est <- all.data[[1]] 
  data.pred <- all.data[[2]] 
  
  ##########################################
  #### Run irregular voronoi block-NNGP models 
  ##########################################

  case='NNGP'
  
  name.cov="matern"
  
  datafill <- blockNNGP.struct(case= 'NNGP', data.est,
                               n.blocks=1, num.nb, 
                               name.cov,par.cov = 0.5,data.pred )
  
  #  the data is reordered...
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
  
  indv1 <- seq(1,k*nloc,by=2)
  xnew1 <-  c(rep(NA, nloc),rep(NA, nloc))
  xnew1[indv1] <- data1$X.2
  xnew2 <-  c(rep(NA, nloc),rep(NA, nloc))
  xnew2[indv1+1] <- data1$X.2
  int1 <-  c(rep(NA, nloc),rep(NA, nloc))
  int1[indv1] <- 1
  int2 <-  c(rep(NA, nloc),rep(NA, nloc))
  int2[indv1+1] <- 1
  
  k <- 2
  spatial.field.y = 1:(k*nloc)
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  
  indg  <- rep(1:2, nloc) 
  index <- rep(1:nloc, each= k) 
  data2 <- list(y = ynew, x1 = xnew1, x2 = xnew2, int1 = int1, int2 = int2,
               spatial.field.y = spatial.field.y, w = wnew, 
               blocks = blocks, index = index, indg = indg,
               loc.1 = loc1, loc.2 = loc2)
  
  # set by user
  # P(rho<rho_0) = alpha1
  prior.range   <-  c( 0.5*0.15,0.05)  # P(rho<0.5*0.15) = 0.05
  # P(sigma>sigma0) = alpha2
  prior.sigma  <-  c(5 ,0.05) # P(sigma>5) = 0.05
  
  
  ## set the block-NNGP as latent effect
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
  
 
  NNGP.model <- inla.rgeneric.define(inla.rgeneric.NNGP.model.pc.prior2,
                                     k = 2, W = W,
                                     lam1 = lam1, lam2 = lam2,
                                     n = nloc,
                                     coords.D = as.matrix(datafill[[3]], ncol, ncol), 
                                     AdjMatrix = as.matrix(datafill[[4]], ncol, ncol), 
                                     initial.range = initial.range,
                                     initial.sigma = initial.sigma)
  
  
  f.NNGP <- y ~ 0+ int1 + int2 + x1+ x2 + 
    f(spatial.field.y, model = NNGP.model) 
  
  # inla function to fit the model
  
  res6 <- inla(f.NNGP, data = (data2),
               family = "poisson",
                control.predictor= list(compute = TRUE),
               control.compute = list(waic = TRUE, cpo = TRUE,config=TRUE),
               verbose=TRUE)
  
  res <- summary.blockNNGP_LMC(name.prior="pc.prior", resf=res6, data2,n.blocks=1, num.nb,
                               family="poisson")
  
  #  the data is reordered...
  datap <- datafill[[5]]
  names(datap)
  
  saida.pred <- NNGP_predLMC(case = 'NNGP', num.nb,
                             data.est = data2, pred.data = datap,
                             res = res6,n.sample = 1000, family = "poisson")
  
  
  orig.param <- all.data[[3]]
  MSE.Y1.res <- round(mean((res6$summary.fitted.values[indv1,1]-data2$y[indv1])^2),3)
  MSE.Y2.res <- round(mean((res6$summary.fitted.values[(indv1+1),1]-data2$y[(indv1+1)])^2),3)
  MSP.Y1.res <- round(mean((saida.pred$mean.y1-saida.pred$Y.pred)^2),3)
  MSP.Y2.res <- round(mean((saida.pred$mean.y2-saida.pred$Y2.pred)^2),3)
  
  return(list(data.est = data2, data.pred = datap, res6, res, orig.param, 
              saida.pred, MSP.Y1.res, MSP.Y2.res, MSE.Y1.res, MSE.Y2.res))  
  
}


run_meshpoisson = function(all.data,d){

  data.est <- all.data[[1]]
  data.pred <- all.data[[2]] 
  orig.param <- all.data[[3]]

  
  n <- dim(data.est)[1] # number of locations
  q <- 2 # number of outcomes
  k <- 2 # true number of spatial factors used to make the outcomes
  p <- 2 # number of covariates
  coords <-  data.est[,1:2]
  colnames(coords) <- c("Var1", "Var2")
  XX <- matrix(rnorm(n*p), ncol=p)
  XX[,1] <- data.est[,5]
  XX[,2] <- data.est[,6]
  
  
  YY_full <- matrix(rnorm(n*p), ncol=p)
  YY_full[,1] <- data.est[,3]
  YY_full[,2] <- data.est[,4]
  
  YY <- YY_full
  
  simdata <- coords %>%
    cbind(data.frame(Outcome_full=YY_full, 
                     Outcome_obs = YY, 
                     w1 = data.est[,7], w2 = data.est[,8])) 
  
  data.est <- NULL

  mcmc_keep <- 30000 
  mcmc_burn <- 10000
  mcmc_thin <- 1
  
  library(meshed)
  mesh_total_time <- system.time({
    meshout <- spmeshed(y = YY, x = XX, coords = coords, k = 2,
                        family = c("poisson","poisson"), 
                        axis_partition = c(d,d), # 8*8=64 number of blocks
                        n_samples = mcmc_keep, 
                        n_burn = mcmc_burn, 
                        n_thin = mcmc_thin, 
                        settings = list(forced_grid=FALSE), 
                        prior = list(phi=c(1, 30)),
                        verbose=0
    )})
  

  mesh_total_time 
  
  par(mfrow=c(2,2))
  plot(meshout$beta_mcmc[1,1,],col="white", ylab=expression(beta[11]))
  lines(meshout$beta_mcmc[1,1,])
  plot(meshout$beta_mcmc[1,2,],col="white", ylab=expression(beta[12]))
  lines(meshout$beta_mcmc[1,2,])
  plot(meshout$beta_mcmc[2,1,],col="white", ylab=expression(beta[21]))
  lines(meshout$beta_mcmc[2,1,])
  plot(meshout$beta_mcmc[2,2,],col="white", ylab=expression(beta[22]))
  lines(meshout$beta_mcmc[2,2,])
  
  par(mfrow=c(2,2))
  plot(meshout$lambda_mcmc[1,1,],col="white", ylab=expression(lambda[11]))
  lines(meshout$lambda_mcmc[1,1,])
  plot(meshout$lambda_mcmc[2,1,],col="white", ylab=expression(lambda[12]))
  lines(meshout$lambda_mcmc[2,1,])
  plot(meshout$lambda_mcmc[2,2,],col="white", ylab=expression(lambda[22]))
  lines(meshout$lambda_mcmc[2,2,])
  
  par(mfrow=c(1,2))
  plot(meshout$theta_mcmc[1,1,],col="white", ylab=expression(theta[1]))
  lines(meshout$theta_mcmc[1,1,])
  plot(meshout$theta_mcmc[1,2,],col="white", ylab=expression(theta[2]))
  lines(meshout$theta_mcmc[1,2,])
  
  beta00 <- c(summary(meshout$beta_mcmc[1,1,])[c(4)], 
              sd(meshout$beta_mcmc[1,1,]), summary(meshout$beta_mcmc[1,1,])[c(2,3,5)])
  beta01 <- c(summary(meshout$beta_mcmc[2,1,])[c(4)], 
              sd(meshout$beta_mcmc[2,1,]), summary(meshout$beta_mcmc[2,1,])[c(2,3,5)])
  beta02 <- c(summary(meshout$beta_mcmc[1,2,])[c(4)], 
              sd(meshout$beta_mcmc[1,2,]), summary(meshout$beta_mcmc[1,2,])[c(2,3,5)])
  beta22 <- c(summary(meshout$beta_mcmc[2,2,])[c(4)], 
              sd(meshout$beta_mcmc[2,2,]), summary(meshout$beta_mcmc[2,2,])[c(2,3,5)])
  beta <- rbind(beta00, beta01, beta02,beta22)
  colnames(beta)[2] <- "sd"
  
  lambda11 <- c(summary(meshout$lambda_mcmc[1,1,])[c(4)], 
                sd(meshout$lambda_mcmc[1,1,]), summary(meshout$lambda_mcmc[1,1,])[c(2,3,5)])
  lambda21 <- c(summary(meshout$lambda_mcmc[2,1,])[c(4)], 
                sd(meshout$lambda_mcmc[2,1,]), summary(meshout$lambda_mcmc[2,1,])[c(2,3,5)])
  lambda22 <- c(summary(meshout$lambda_mcmc[2,2,])[c(4)], 
                sd(meshout$lambda_mcmc[2,2,]), summary(meshout$lambda_mcmc[2,2,])[c(2,3,5)])
  lambda <- rbind(lambda11,  lambda21,lambda22)
  colnames(lambda)[2] <- "sd"
  
  
  theta11 <- c(summary(meshout$theta_mcmc[1,1,])[c(4)], 
               sd(meshout$theta_mcmc[1,1,]), summary(meshout$theta_mcmc[1,1,])[c(2,3,5)])
  theta12 <- c(summary(meshout$theta_mcmc[2,1,])[c(4)], 
               sd(meshout$theta_mcmc[2,1,]), summary(meshout$theta_mcmc[2,1,])[c(2,3,5)])
  theta21 <- c(summary(meshout$theta_mcmc[1,2,])[c(4)], 
               sd(meshout$theta_mcmc[1,2,]), summary(meshout$theta_mcmc[1,2,])[c(2,3,5)])
  theta22 <- c(summary(meshout$theta_mcmc[2,2,])[c(4)], 
               sd(meshout$theta_mcmc[2,2,]), summary(meshout$theta_mcmc[2,2,])[c(2,3,5)])
  theta <- rbind(theta11, theta12, theta21, theta22)
  colnames(theta)[2] <- "sd"
  
  ret <- rbind(beta, lambda, theta)
  
  my_abind = function(arg.list){
    m <- 0
    for (i in seq(length.out=length(arg.list))) {
      print(i)
      m1 <- arg.list[[i]]
      m <- m + m1
    }
    mf <- m/length(arg.list)  
    return(mf)
  }
  
  
  arg.list = meshout$yhat_mcmc
  y_post_sample <- my_abind(arg.list)
  
  perf1 <- meshout$coordsdata %>% 
    cbind(y_pm_1 = y_post_sample[,1]) %>%
    left_join(simdata, by = c("Var1", "Var2"))
  
  y1= perf1$Outcome_full.1 
  esty1=perf1$y_pm_1 
  
  MSE.Y1.res <- round(mean((y1-esty1)^2),3)
  
  perf2 <- meshout$coordsdata %>% 
    cbind(y_pm_2 = y_post_sample[,2]) %>%
    left_join(simdata, by = c("Var1", "Var2"))
  
  y2 <- perf2$Outcome_full.2 
  
  esty2 <- perf2$y_pm_2 
  
  MSE.Y2.res <- round(mean((y2-esty2)^2),3)
  
  par(mfrow = c(1,2))
  plot(y1,esty1)
  plot(y2,esty2)
  
  MSE.Y.res <- (MSE.Y1.res + MSE.Y2.res)/2
  
  
  
 
  res.plot <- data.frame(y1,  esty1, y2, esty2 )
  
  p1 <- ggplot(res.plot, aes(y1, esty1)) + 
    geom_point()+ 
    labs(x=expression(y[1]),y=expression(E(y[1]*"|"*y)))
  p2 <- ggplot(res.plot, aes(y2, esty2)) + 
    geom_point()+
    labs(x=expression(y[2]),y=expression(E(y[2]*"|"*y)))
 
  library(ggpubr)
  ggarrange( p1, p2, ncol = 2, nrow = 1)
  
  MSE.w1.res <- NA
  MSE.w2.res <- NA
  
  return(list(data.est, data.pred, meshout, res=ret, orig.param, saida.pred=NULL, MSP.Y1.res=NULL, MSP.Y2.res=NULL, MSE.Y1.res, MSE.Y2.res, MSE.w1.res, MSE.w2.res))  
}

# posterior results for blockNNGP and NNGP 

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
  
  # orig.param<- orig.param[c(1:4, 6:10)] #for poisson
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
  
  
  name.crit <- c("waic", "LPML","time(sec)","MSP.Y1","MSP.Y2","MSE.Y1.res", "MSE.Y2.res") 
  crit.model<- round(c(  waic, LPML,  time, MSP.Y1, MSP.Y2,MSE.Y1.res,MSE.Y2.res),3)
  name.full.results <- c(t(names.mean.estparam), t(names.cov.estparam), t(name.crit))
  full.results <- c(t(mean.estparam), t(cov.param), t(crit.model) )
  final.res <- data.frame(res=name.full.results, val=full.results)
   
  return(final.res)
  
}

# posterior results for QMGP
post_resultsmeshed = function(resf){
  resmesh <- resf[[3]]
   time <- resmesh$mcmc_time
   
  
  MSE.Y1.res <- resf[[9]]
  MSE.Y2.res <- resf[[10]]
  MSP.Y1 <- resf[[7]]
  MSP.Y2 <- resf[[8]]
  
   
  orig.param <- resf[[5]]
  hyperpar <- resf[[4]][5:11,]
  fixed <- resf[[4]][1:4,]
  lambda <- resf[[4]][6,]
  
  est.beta00 <-fixed[1,1]
  est.beta01 <-fixed[2,1]
  est.beta02 <-fixed[3,1]
  est.beta22 <-fixed[4,1]
  est.sigma21 <-hyperpar[1,1]
  est.phi1 <-hyperpar[4,1]
  est.sigma22 <-hyperpar[3,1]
  est.phi2 <-hyperpar[6,1]
  est.lambda <- lambda[1]
   
  covbeta00  <- 0; if(orig.param[1]>= fixed[1,3] & orig.param[1]<= fixed[1,5]) covbeta00 <- 1
  covbeta01  <- 0; if(orig.param[2]>= fixed[2,3] & orig.param[2]<= fixed[2,5]) covbeta01 <- 1
  covsigma21 <- 0; if(orig.param[3]>= hyperpar[1,3] & orig.param[3]<= hyperpar[1,5]) covsigma21 <- 1
  covphi1    <- 0; if(orig.param[4]>=  hyperpar[4,3] & orig.param[4]<=  hyperpar[4,5]) covphi1 <- 1
  
  covbeta02  <- 0; if(orig.param[5]>= fixed[3,3] & orig.param[5]<= fixed[3,5]) covbeta02 <- 1
  covbeta22  <- 0; if(orig.param[6]>= fixed[4,3] & orig.param[6]<= fixed[4,5]) covbeta22 <- 1
  covsigma22 <- 0; if(orig.param[7]>= hyperpar[3,3] & orig.param[7]<= hyperpar[3,5]) covsigma22 <- 1
  covphi2    <- 0; if(orig.param[8]>=  hyperpar[6,3] & orig.param[8]<=  hyperpar[6,5]) covphi2 <- 1
  
  covlambda <- 0; if(orig.param[9]>=  lambda[3] & orig.param[9]<=  lambda[5]) covlambda <- 1
  names.mean.estparam <- c("est.beta00", "est.beta01", "est.beta02","est.beta22", "est.sigma21", "est.phi1", "est.sigma22", "est.phi2", "est.lambda")
  mean.estparam <- round(c(est.beta00, est.beta01, est.beta02, est.beta22,est.sigma21, est.phi1, est.sigma22, est.phi2, est.lambda),3)
  names.cov.estparam <- c("covbeta00", "covbeta01", "covbeta02","covbeta22","covsigma21", "covphi1", "covsigma22", "covphi2", "covlambda")
  cov.param <- c(covbeta00, covbeta01, covbeta02, covbeta22,covsigma21, covphi1, covsigma22, covphi2,covlambda)
  
  
  name.crit <- c("waic", "LPML","time(sec)","MSP.Y1","MSP.Y2","MSE.Y1.res", "MSE.Y2.res","MSE.w1.res", "MSE.w2.res") #, "RMSPy")
  crit.model<- round(c(  0, 0,  time, 0, 0,MSE.Y1.res,MSE.Y2.res,0, 0),3) #, RMSP
  name.full.results <- c(t(names.mean.estparam), t(names.cov.estparam), t(name.crit))
  full.results <- c(t(mean.estparam), t(cov.param), t(crit.model) )
  final.res <- data.frame(res=name.full.results, val=full.results)
 
  return(final.res)
  
}


# blockNNGP
runblockNNGP_M_n = function(all.data, n.blocks,num.nb){
  res1 <- run_poissonIrreg0(all.data, n.blocks, num.nb)
  final.res <- post_results(res1)
  name1 <-  paste(n,"_", n.blocks, "_", num.nb,sep ="")
  save(res1 = res1, file=paste("simblockNNGP_",name1,"_pois.Rdata",sep =""))
  save(final.res=final.res, file=paste("ALLblockNNGP_sim_",name1,"_pois.Rdata", sep =""))
  res1 <- NULL
  return(final.res)
}

#NNGP
runNNGP_M_n = function(all.data, num.nb){
res1 <- run_poisson_NNGP0(all.data,  num.nb)
final.res <- post_results(res1)
name1 <-  paste(n,"_", num.nb, sep = "")
save(res1 = res1, file = paste("simNNGP_", name1, "_pois.Rdata", sep = ""))
save(final.res = final.res, file = paste("ALLNNGP_sim_",name1,"_pois.Rdata", sep = ""))
res1 <- NULL
return(final.res)
}

# QMGP
runQMGP_d = function(all.data, d){
  res1 <- run_meshpoisson(all.data, d)
  final.res <- post_resultsmeshed(res1)
  name1 <-  paste(n, "_", d, sep = "")
  save(res1 = res1, file = paste("simQMGP_",name1,"_pois.Rdata", sep = ""))
  save(final.res = final.res, file = paste("ALLQMGP_sim_", name1, "_pois.Rdata", sep = ""))
  res1 <- NULL
  return(final.res)
}



##############################################
## run ##
##############################################
set.seed <- 2

# sample size n 
n         <- 2500
all.data  <- sim_poisson.data2(m=set.seed, n)
block2500 <- runblockNNGP_M_n(all.data,n.blocks = 64,num.nb = 2)
NNGP2500 <- runNNGP_M_n(all.data, num.nb = 10)
meshed2500 <- runQMGP_d(all.data, d = 8)

n         <- 5000
all.data  <- sim_poisson.data2(m=set.seed, n)
block5000 <- runblockNNGP_M_n(all.data,n.blocks = 144,num.nb = 2)
NNGP5000 <- runNNGP_M_n(all.data, num.nb = 20)
meshed5000 <- runQMGP_d(all.data, d = 12)


n         <- 10000
all.data  <- sim_poisson.data2(m = set.seed, n)
block10000 <- runblockNNGP_M_n(all.data,n.blocks = 225,num.nb = 2)
NNGP10000 <- runNNGP_M_n(all.data, num.nb = 30)
meshed10000 <- runQMGP_d(all.data, d = 15)

n         <- 15000
all.data  <- sim_poisson.data2(m = set.seed, n)
block12000 <- runblockNNGP_M_n(all.data,n.blocks = 400,num.nb = 2)
NNGP12000 <- runNNGP_M_n(all.data, num.nb = 30)
meshed12000 <- runQMGP_d(all.data, d = 20)


# plots

load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLblockNNGP_sim_2500_64_2_pois.Rdata")
block2500 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLblockNNGP_sim_5000_144_2_pois.Rdata")
block5000 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLblockNNGP_sim_10000_225_2_pois.Rdata")
block10000 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLblockNNGP_sim_15000_400_2_pois.Rdata")
block15000 <- final.res
#block15000 <- final.res
final.res <- NULL
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLNNGP_sim_2500_10_pois.Rdata")
NNGP2500 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLNNGP_sim_5000_20_pois.Rdata")
NNGP5000 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLNNGP_sim_10000_30_pois.Rdata")
NNGP10000 <- final.res
final.res <- NULL
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLQMGP_sim_2500_8_pois.Rdata")
meshed2500 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLQMGP_sim_5000_12_pois.Rdata")
meshed5000 <- final.res
load("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/sim3/ALLQMGP_sim_10000_15_pois.Rdata")
meshed10000 <- final.res

final.res <- NULL

time_block2500 <- block2500[21,2]
time_block5000 <- block5000[21,2]
time_block10000 <- block10000[21,2]
time_block15000 <- block15000[21,2]


MSE_block2500 <- block2500[24,2]
MSE_block5000 <- block5000[24,2]
MSE_block10000 <- block10000[24,2]
MSE_block15000 <- block15000[24,2]
#MSE_block15000 <- block15000[24,2]

MSE2_block2500 <- block2500[25,2]
MSE2_block5000 <- block5000[25,2]
MSE2_block10000 <- block10000[25,2]
MSE2_block15000 <- block15000[25,2]


MSP_block2500 <- block2500[22,2]
MSP_block5000 <- block5000[22,2]
MSP_block10000 <- block10000[22,2]
MSP_block15000 <- block15000[22,2]

MSP2_block2500 <- block2500[23,2]
MSP2_block5000 <- block5000[23,2]
MSP2_block10000 <- block10000[23,2]
MSP2_block15000 <- block15000[23,2]

MSEw1_block2500 <- block2500[26,2]
MSEw1_block5000 <- block5000[26,2]
MSEw1_block10000 <- block10000[26,2]
MSEw1_block15000 <- block15000[26,2]
#MSEw1_block15000 <- block15000[26,2]

MSEw2_block2500 <- block2500[27,2]
MSEw2_block5000 <- block5000[27,2]
MSEw2_block10000 <- block10000[27,2]
MSEw2_block15000 <- block15000[27,2]


time_NNGP2500 <- NNGP2500[21,2]
time_NNGP5000 <- NNGP5000[21,2]
time_NNGP10000 <- NNGP10000[21,2]


MSE_NNGP2500 <- NNGP2500[24,2]
MSE_NNGP5000 <- NNGP5000[24,2]
MSE_NNGP10000 <- NNGP10000[24,2]


MSE2_NNGP2500 <- NNGP2500[25,2]
MSE2_NNGP5000 <- NNGP5000[25,2]
MSE2_NNGP10000 <- NNGP10000[25,2]


MSP_NNGP2500 <- NNGP2500[22,2]
MSP_NNGP5000 <- NNGP5000[22,2]
MSP_NNGP10000 <- NNGP10000[22,2]


MSP2_NNGP2500 <- NNGP2500[23,2]
MSP2_NNGP5000 <- NNGP5000[23,2]
MSP2_NNGP10000 <- NNGP10000[23,2]


MSEw1_NNGP2500 <- NNGP2500[26,2]
MSEw1_NNGP5000 <- NNGP5000[26,2]
MSEw1_NNGP10000 <- NNGP10000[26,2]


MSEw2_NNGP2500 <- NNGP2500[27,2]
MSEw2_NNGP5000 <- NNGP5000[27,2]
MSEw2_NNGP10000 <- NNGP10000[27,2]


time_meshed2500 <- meshed2500[21,2]
time_meshed5000 <- meshed5000[21,2]
time_meshed10000 <- meshed10000[21,2]


MSE_meshed2500 <- meshed2500[24,2]
MSE_meshed5000 <- meshed5000[24,2]
MSE_meshed10000 <- meshed10000[24,2]

MSE2_meshed2500 <- meshed2500[25,2]
MSE2_meshed5000 <- meshed5000[25,2]
MSE2_meshed10000 <- meshed10000[25,2]

MSP_meshed2500 <- meshed2500[22,2]
MSP_meshed5000 <- meshed5000[22,2]
MSP_meshed10000 <- meshed10000[22,2]

MSP2_meshed2500 <- meshed2500[23,2]
MSP2_meshed5000 <- meshed5000[23,2]
MSP2_meshed10000 <- meshed10000[23,2]

MSEw1_meshed2500 <- meshed2500[26,2]
MSEw1_meshed5000 <- meshed5000[26,2]
MSEw1_meshed10000 <- meshed10000[26,2]

MSEw2_meshed2500 <- meshed2500[27,2]
MSEw2_meshed5000 <- meshed5000[27,2]
MSEw2_meshed10000 <- meshed10000[27,2]

# join times
all.time <- c(time_block2500, time_block5000, time_block10000, #time_block12000, time_block15000,
              time_NNGP2500, time_NNGP5000, time_NNGP10000, #time_NNGP12000,
              time_meshed2500, time_meshed5000, time_meshed10000)

all.MSE <- c(MSE_block2500, MSE_block5000, MSE_block10000, #MSE_block12000, MSE_block15000,
             MSE_NNGP2500, MSE_NNGP5000, MSE_NNGP10000, #MSE_NNGP12000,
             MSE_meshed2500, MSE_meshed5000, MSE_meshed10000)

all.MSE2 <- c(MSE2_block2500, MSE2_block5000, MSE2_block10000, #MSE2_block12000, MSE2_block15000,
              MSE2_NNGP2500, MSE2_NNGP5000, MSE2_NNGP10000, #MSE2_NNGP12000,
              MSE2_meshed2500, MSE2_meshed5000, MSE2_meshed10000)

all.MSP <- c(MSP_block2500, MSP_block5000, MSP_block10000, #MSP_block12000, MSP_block15000,
             MSP_NNGP2500, MSP_NNGP5000, MSP_NNGP10000, #MSP_NNGP12000,
             MSP_meshed2500, MSP_meshed5000, MSP_meshed10000)

all.MSP2 <- c(MSP2_block2500, MSP2_block5000, MSP2_block10000, #MSP2_block12000, MSP2_block15000,
              MSP2_NNGP2500, MSP2_NNGP5000, MSP2_NNGP10000, #MSP2_NNGP12000,
              MSP2_meshed2500, MSP2_meshed5000, MSP2_meshed10000)

all.MSEw1 <- c(MSEw1_block2500, MSEw1_block5000, MSEw1_block10000,# MSEw1_block12000, MSEw1_block15000,
               MSEw1_NNGP2500, MSEw1_NNGP5000, MSEw1_NNGP10000, #MSEw1_NNGP12000,
               MSEw1_meshed2500, MSEw1_meshed5000, MSEw1_meshed10000)

all.MSEw2 <- c(MSEw2_block2500, MSEw2_block5000, MSEw2_block10000, #MSEw2_block12000, MSEw2_block15000,
               MSEw2_NNGP2500, MSEw2_NNGP5000, MSEw2_NNGP10000, #MSEw2_NNGP12000,
               MSEw2_meshed2500, MSEw2_meshed5000, MSEw2_meshed10000)

model <- c("blockNNGP-GLMC","blockNNGP-GLMC","blockNNGP-GLMC",#"GLMC-blockNNGP","GLMC-blockNNGP",
           "NNGP-GLMC","NNGP-GLMC","NNGP-GLMC",#"GLMC-NNGP",
           "QMGP", "QMGP", "QMGP")

samplesize <- c(4000,8000,16000,  
                4000,8000,16000, 
                4000,8000,16000)

plot.time <- data.frame(model=model, time=all.time)

plot.MSE <- data.frame(model=model, MSE1=all.MSE)

plot.MSE2 <- data.frame(model=model, MSE2=all.MSE2)

plot.MSP <- data.frame(model=model, MSE1=all.MSP)

plot.MSP2 <- data.frame(model=model, MSE2=all.MSP2)

all.plot <- data.frame(model=model, time=all.time, 
                       MSEy1=all.MSE, MSEy2=all.MSE2,
                       MSPy1=all.MSP, MSPy2=all.MSP2,
                       N=samplesize)

library(ggplot2)

p1 <- ggplot(all.plot, aes(x=N, y = time, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)
#p1


p2 <- ggplot(all.plot, aes(x=N, y = MSEy1, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)
#p4

p3 <- ggplot(all.plot, aes(x=N, y = MSEy2, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)
#p5

p4 <- ggplot(all.plot, aes(x=N, y = MSPy1, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)
#p4

p5 <- ggplot(all.plot, aes(x=N, y = MSPy2, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)

library(ggpubr)

ggarrange( p2, p3, ncol = 2, nrow = 1,
           common.legend=T,
           legend="right")

ggarrange( p4, p5, ncol = 2, nrow = 1,
           common.legend=T,
           legend="right")

finalplot0 <- ggarrange( p1, p2, p3, ncol = 3, nrow = 1,
           common.legend=T,
           legend="right",
           labels = c("a)", "b)", "c)"),
           font.label = list(size = 16, color = "black", face = "bold", family = NULL)) 

ggsave("finalplot0.eps", plot = finalplot0, device = "eps") 


all.plot2<- all.plot[-(4:6),]


p1N <- ggplot(all.plot2, aes(x=N, y = time, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)+
  theme(legend.position = "bottom")
#p1

p2N <- ggplot(all.plot2, aes(x=N, y = MSEy1, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)+
  theme(legend.position = "bottom")
#p4

p3N <- ggplot(all.plot2, aes(x=N, y = MSEy2, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)+
  theme(legend.position = "bottom")
#p5
p4N <- ggplot(all.plot2, aes(x=N, y = MSPy1, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)+
  theme(legend.position = "bottom")
#p4

p5N <- ggplot(all.plot2, aes(x=N, y = MSPy2, color=model)) +
  geom_point(aes(shape=model), size = 3) +
  geom_line(aes(linetype=model), linewidth = 1)+
  theme_bw(base_size=12)+
  theme(legend.position = "bottom")
#p5

library(ggpubr)

finalplot <- ggarrange(p1N, p2N, p3N, ncol = 3, nrow = 1, 
           common.legend=T,
           legend="right")
ggsave("finalplot.eps", plot = finalplot, device = "eps") 
