############################################################################
## Date: 06.05.2026
##
## Description:
##
##    The code performs a Bayesian analysis generalized block-NNGP 
##    models  through INLA. These models are applied to Birds 
##    sampling data. It is assumed a Poisson 
##    distribution for the response variables (bird abundances of  
##    robin, Sparrow, starling). 
##    @author:  Zaida Quiroz (\email{zquiroz@pucp.edu.pe}).
##
rm(list=ls())

getwd()

setwd("~/my_directory")

#######################
## functions required
####################### 
source("blockNNGPfunctions.R")
source("newQGLMCrgeneric.R")
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
library(plotly)



run_poissonIrregK = function(datafill, n.blocks, num.nb, k, model){
  
  # Ordering data
  data1 <- datafill[[1]]
  names(data1)[c(1:6, 11)] <- c("loc.1", "loc.2","y1", "y2","y3",
                   "X.1","blocks")
  nloc <- dim(data1)[1]
  
  sq  <- seq(1, k*nloc, k)
  ind <- rep(0, k*nloc)
  k2  <- 1
  for(i in(sq)){
    ind[i:(i+(k-1))] <- seq(k2, k*nloc, nloc)
    k2 <- k2+1
  }
  
  L1 <- 3:(2+k)
  L2 <- (5+k):(4+2*k)
  yf <- NULL
  for(j in(1:k)){
    yf1 <- data1[,L1[j]]
    yf  <- c(yf, yf1)
  }
  
  ynew   <-  yf[ind]
  blocks <- rep(data1$blocks, k)
  blocks <-  blocks[ind]
  loc1 <- rep(data1$loc.1, k)
  loc1 <- loc1[ind]
  loc2 <- rep(data1$loc.2, k)
  loc2 <- loc2[ind]
  
  indv1 <-seq(1, k*nloc, by=k)
  xnew1 <-  rep(rep(NA, nloc), k)
  xnew1[indv1] <- data1$loc.2  
  xnew2 <-  rep(rep(NA, nloc), k)
  xnew2[indv1+1] <- (data1$loc.2) 
  
  int1 <-  rep(rep(NA, nloc), k)
  int1[indv1] <- 1
  int2 <-  rep(rep(NA, nloc), k)
  int2[indv1+1] <- 1
  
  xnew3 <- 0
  int3  <- 0

    xnew3 <-  rep(rep(NA, nloc), k)
    xnew3[indv1+2] <- data1$loc.2  
    int3 <-  rep(rep(NA, nloc), k)
    int3[indv1+2] <- 1

  spatial.field.y <- 1:(k*nloc)
  
  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  indg  <- rep(1:k, nloc) 
  index <- rep(1:nloc, each= k) 
  data2 <- list(y = ynew, 
               x1 = xnew1, x2 = xnew2, x3 = xnew3, 
               x1.2 = xnew1^2, x2.2 = xnew2^2, x3.2 = xnew3^2, 
               int1 = int1, int2 = int2, int3 = int3, 
               spatial.field.y = spatial.field.y,
               blocks = blocks, index = index, indg = indg,
               loc.1 = loc1, loc.2 = loc2)
  
  # set by user
  # P(rho<rho_0) = alpha1
  prior.range  = c( 80,0.95)  # P(rho<80) = 0.95
  # P(sigma>sigma0) = alpha2
  prior.sigma  = c(1000 ,0.05) # P(sigma>1000) = 0.05
  
  
  ## set the block-NNGP as latent effect
  #  name.prior="pc.prior"
  prior.set   <- c(prior.range,prior.sigma)
  prior.range <- c(prior.set[1:2])
  prior.sigma <- c(prior.set[3:4])
  
  lam1 <- -log(prior.range[2])*prior.range[1]
  lam2 <- -log(prior.sigma[2])/prior.sigma[1]
  
  initial.range <- log(prior.range[1]) + 1
  initial.sigma <- log(prior.sigma[1]) - 1
  
  W1  <- datafill[[2]]
  

  maskWk = function(k,W1){
    M <- diag(1, k)
    M [upper.tri(M )] <- rep(0.5,k*(k-1)/2) # just some value for W mask
    invM <- solve(M)
    listQnew <- NULL
    Wn <- W1
    Q  <- 0
    for (j in(1:k)){
      gammajinv <- matrix(invM[,j], k,1)
      Qnewj <- kronecker(W1, gammajinv%*%t(gammajinv))
      Q <- Q + Qnewj
    }
    W <- Q*0
    mask_matrix <- Q!=0
    W[mask_matrix] <-1
    W <- as(W, "sparseMatrix")
    return(W)
  }
  
  W <- maskWk(k,W1)
  rm(W1)
  rm(all.data)
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
  
  rm(W)
  

    f.blockNNGP <- y ~ 0 + int1 + int2 + int3 + 
      x1 + x2 +   x3 + x1.2 + x2.2 + x3.2 +
      f(spatial.field.y, model = blockNNGP.model) 
 
  
  # inla function to fit the model  with cpo  and waic 
  res4 <- inla(f.blockNNGP, data = (data2),
               family = "poisson",
               control.compute = list(waic = TRUE, cpo = TRUE),
               verbose = TRUE)
  summary(res4)

  # Ordering data
  datap <- datafill[[11]]
  names(datap)[1:5] <- c("loc.1", "loc.2", "y.1", "y.2", "y.3")
  
  res6 <- list(data.est=data2, data.pred=datap, res4,datafill)
  
  return( res6)  
  
}

runforpred =function(datafill,data2,k,model){
  
  # Ordering data
  data1 <- datafill[[1]]
  names(data1)
  names(data1)[c(1:6, 11)] <- c("loc.1", "loc.2","y1", "y2","y3",
                   "X.1","blocks")
  nloc <- dim(data1)[1]
  
  sq  <- seq(1, k*nloc, k)
  ind <- rep(0, k*nloc)
  k2  <- 1
  for(i in(sq)){
    ind[i:(i+(k-1))] <- seq(k2,k*nloc, nloc)
    k2 <-k2+1
  }
  
  # set by user
  # P(rho<rho_0) = alpha1
  prior.range  <- c( 80,0.95)  # P(rho<80) = 0.95
  # P(sigma>sigma0) = alpha2
  prior.sigma  <- c(1000 ,0.05) # P(sigma>1000) = 0.05
  
  
  ## set the block-NNGP as latent effect
  prior.set   <- c(prior.range, prior.sigma)
  prior.range <- c(prior.set[1:2])
  prior.sigma <- c(prior.set[3:4])
  
  lam1 <- -log(prior.range[2])*prior.range[1]
  lam2 <- -log(prior.sigma[2])/prior.sigma[1]
  
  initial.range <- log(prior.range[1]) + 1
  initial.sigma <- log(prior.sigma[1]) - 1
  
  W1 <- datafill[[2]]
  

  maskWk = function(k,W1){
    M <- diag(1, k)
    M [upper.tri(M )] <- rep(0.5,k*(k-1)/2) # just some value for W mask
    invM <- solve(M)
    listQnew <- NULL
    Wn <- W1
    Q <- 0
    for (j in(1:k)){
      gammajinv <- matrix(invM[,j], k, 1)
      Qnewj <- kronecker(W1, gammajinv%*%t(gammajinv))
      Q <- Q + Qnewj
    }
    W <- Q*0
    mask_matrix <- Q!=0
    W[mask_matrix] <-1
    W <- as(W, "sparseMatrix")
    return(W)
  }
  
  W <- maskWk(k,W1)
  rm(W1)
  
  blockNNGP.model <- inla.rgeneric.define(inla.rgeneric.LMCk2_QblockNNGP.model.pc.prior, 
                                          k = k,    
                                          W = W,
                                          lam1 = lam1,
                                          lam2 = lam2,
                                          n = datafill[[3]], n.blocks = datafill[[4]],nb = datafill[[5]],
                                          ind_obs1 = datafill[[6]], num1 = datafill[[7]],indb = datafill[[8]],
                                          coords.D = datafill[[9]], a = prior.set[1], b = prior.set[2],
                                          initial.range = initial.range,
                                          initial.sigma = initial.sigma)
  
  rm(W)
  
    f.blockNNGP <- y ~ 0 + int1 + int2 + int3 + 
      x1 + x2 +   x3 + x1.2 + x2.2 + x3.2 +
      f(spatial.field.y, model = blockNNGP.model) 
  
    # for predictions calling config=TRUE
  res5 <- inla(f.blockNNGP, data = (data2),
               family = "poisson",
               control.compute = list(config=TRUE),
               verbose = TRUE)

  n.sample <- 1000

  param.pos <-   inla.posterior.sample(n.sample, res5)
  
  save(param.pos, file = paste(model,"param.pos", case, n.blocks,"-", num.nb,".Rdata",sep = ''))
  return(res5)
}

# summary for estimations
summary_block1 = function(res4, res5,n.blocks, num.nb,data2){
  k   <- 3
  res <- summary.blockNNGP_LMC(name.prior = "pc.prior", 
                               resf = res4, 
                               n.blocks, num.nb, 
                               family = "poisson")
  nloc <- length(data2$y)/3
  indv1 <- seq(1, k*nloc, by = k)
  
  MSE.Y <- NULL
  saida.est <- NULL
  for(j in (1:k)){
    MSE.Yj.res <- round(mean((res4$summary.fitted.values[indv1+j-1,1]-data2$y[indv1+j-1])^2), 3)
    MSE.Y <- c(MSE.Y, MSE.Yj.res)
    saida.estj <- cbind(data2$y[indv1+j-1],res4$summary.fitted.values[indv1+j-1,1])
    saida.est<- cbind(saida.est, saida.estj)
  }
  
  par(mfrow = c(1,3))
  plot(saida.est[,1],saida.est[,2])
  plot(saida.est[,3],saida.est[,4])
  plot(saida.est[,5],saida.est[,6])
  
  
  res6 <- list( res,  MSE.Y)
  return(res6)
}

# summary for predictions
summary_block2 = function(n.blocks, num.nb, param.pos,res4){
  saida.pred <- blockNNGP_predLMC(case = 'irregular', n.blocks, num.nb,
                                  data.est = res4[[1]], pred.data = res4[[2]],
                                  n.sample = 1000,
                                  family = "poisson", k = 3, param.pos,
                                  AdjMatrix = res4[[4]]$AdjMatrix)

  k     <- 3
  MSP.Y <- NULL
  for(j in (1:k)){
    MSP.Yj.res <- round(mean((saida.pred[,j+k]-saida.pred[,j])^2, na.rm=T), 3)
    MSP.Y <- c(MSP.Y, MSP.Yj.res)
  }
  
  par(mfrow=c(1,3))
  j <- 1
  plot(saida.pred[,j+k], saida.pred[,j])
  j <- 2
  plot(saida.pred[,j+k], saida.pred[,j])
  j <- 3
  plot(saida.pred[,j+k], saida.pred[,j])
  
  res6 <- list(saida.pred, MSP.Y)
  return(res6)
}


##########################################
#M1: robin, Sparrow, starling
#M2: robin, starling, Sparrow
#M3: Sparrow, robin , starling 
#M4: Sparrow, starling, robin 
#M5: starling, robin, Sparrow
#M6: starling, Sparrow, robin
###########################################

setwd("~/Documents/ZAIDA_MAC/PAPERS/multivariateblockNNGP/newCODES/app4-1/M1/")
dir.save = getwd()
load( paste("data.original0.Rdata",sep=""))

#  data for estimation and validation/prediction
nloc <- dim(data.original)[1]
indsample <- sample( 1:nloc, round(0.8*nloc))
data.est <- data.original[indsample,]
data.pred <- data.original[-indsample,]

all.data1 <- list(data.est, data.pred)
save(all.data1, file = paste("all.data1.Rdata", sep = ""))
rm(data.est, data.pred, data.original, nloc, indsample)

# creating blocks
case <- 'irregular'
n.blocks    	<- 100
num.nb 	<- 2
k <- 3

dataf <-  data_proc(all.data = all.data1, n.blocks, num.nb, k)
rm(all.data1)

# M1: robin, Sparrow, starling
resI1 <- run_poissonIrregK(dataf, n.blocks, num.nb, k, model = "M1")
save(resI1, file = paste("M1",case, n.blocks,"-", num.nb,"resI1.Rdata",sep = ''))
rm(resI1)


# M2: robin(3),starling(5),Sparrow(4)
dataf2 <- dataf
dataf2[[1]] <- dataf[[1]][,c(1,2, 3,5,4,6:11)]
resI2 <- run_poissonIrregK(dataf2, n.blocks, num.nb,k,model = "M2")
save(resI2, file=paste("M2", case, n.blocks,"-", num.nb,"resI2.Rdata", sep = ''))
rm(resI2, dataf2)


#M3: Sparrow(4),robin(3),starling(5)
dataf3 <- dataf
dataf3[[1]] <- dataf[[1]][,c(1,2,4,3,5,6:11)]
resI3 <- run_poissonIrregK(dataf3, n.blocks, num.nb, k, model = "M3")
save(resI3, file = paste("M3", case, n.blocks, "-", num.nb, "resI3.Rdata", sep = ''))
rm(resI3, dataf3)

#M4: Sparrow(4),starling(5),robin(3)
dataf4 <- dataf
dataf4[[1]] <- dataf[[1]][,c(1,2,4,5,3,6:11)]
resI4 <- run_poissonIrregK(dataf4, n.blocks, num.nb, k,model = "M4")
save(resI4, file = paste("M4", case, n.blocks, "-", num.nb, "resI4.Rdata", sep = ''))
rm(resI4, dataf4)


#M5: starling(5),robin(3),Sparrow(4)
dataf5 <- dataf
dataf5[[1]] <- dataf[[1]][,c(1,2,5,3,4,6:11)]
resI5 <- run_poissonIrregK(dataf5, n.blocks, num.nb, k, model = "M5")
save(resI5, file=paste("M5", case, n.blocks,"-", num.nb,"resI5.Rdata",sep = ''))
rm(resI5,dataf5)


#M6: starling(5),Sparrow(4),robin(3)
dataf6 <- dataf
dataf6[[1]] <- dataf[[1]][,c(1,2,5,4,3,6:11)]
resI6 <- run_poissonIrregK(dataf6, n.blocks, num.nb, k, model = "M6")
save(resI6, file=paste("M6", case, n.blocks, "-", num.nb,"resI6.Rdata", sep = ''))
rm(resI6, dataf6)


# predictions
load(paste(getwd(),"/M1irregular60-2resI1.Rdata", sep = ""))
predI1 <- runforpred(datafill = resI1[[4]], data2 = resI1[[1]], k = 3, model = 'M1')
save(predI1, file = paste("M1", case, n.blocks, "-", num.nb, "predI1.Rdata",sep = ''))
rm(resI1, predI1)

load(paste(getwd(),"/M2irregular60-2resI2.Rdata", sep = ""))
predI2 <- runforpred(datafill = resI2[[4]], data2 = resI2[[1]], k = 3, model = 'M2')
save(predI2, file = paste("M2", case, n.blocks, "-", num.nb, "predI2.Rdata", sep=''))
rm(resI2, predI2)

load(paste(getwd(),"/M3irregular60-2resI3.Rdata", sep = ""))
predI3 <- runforpred(datafill = resI3[[4]], data2 = resI3[[1]], k = 3, model = 'M3')
save(predI3, file = paste("M3", case, n.blocks,"-", num.nb,"predI3.Rdata",sep = ''))
rm(resI3, predI3)

load(paste(getwd(),"/M4irregular60-2resI4.Rdata", sep = ""))
predI4 <- runforpred(datafill = resI4[[4]], data2 = resI4[[1]], k = 3, model = 'M4')
save(predI4, file = paste("M4", case, n.blocks,"-", num.nb, "predI4.Rdata", sep = ''))
rm(resI4, predI4)

load(paste(getwd(),"/M5irregular60-2resI5.Rdata", sep = ""))
predI5 <- runforpred(datafill = resI5[[4]], data2 = resI5[[1]], k = 3, model = 'M5')
save(predI5, file = paste("M5", case, n.blocks, "-", num.nb, "predI5.Rdata",sep = ''))
rm(resI5, predI5)

load(paste(getwd(),"/M6irregular60-2resI6.Rdata", sep = ""))
predI6 <- runforpred(datafill = resI6[[4]], data2 = resI6[[1]], k = 3, model = 'M6')
save(predI6, file = paste("M6", case, n.blocks, "-", num.nb, "predI6.Rdata", sep = ''))
rm(resI6, predI6)

# summary
load(paste(getwd(),"/M1irregular60-2resI1.Rdata", sep = ""))
resF11 <- summary_block1(res4 = resI1[[3]],  num.nb, data2 = resI1[[1]])
save(resF11, file = paste("M1", case, n.blocks,"-", num.nb, "resF11.Rdata",sep = ''))

load(paste(getwd(),"/M2irregular60-2resI2.Rdata", sep = ""))
resF21 <- summary_block1(res4 = resI2[[3]],  num.nb, data2 = resI2[[1]])
save(resF21, file = paste("M2", case, n.blocks, "-", num.nb, "resF21.Rdata", sep = ''))

load(paste(getwd(),"/M3irregular60-2resI3.Rdata", sep = ""))
resF31 <- summary_block1(res4 = resI3[[3]],  num.nb, data2 = resI3[[1]])
save(resF31, file = paste("M3", case, n.blocks, "-", num.nb, "resF31.Rdata", sep = ''))

load(paste( getwd(),"/M4irregular60-2resI4.Rdata", sep = ""))
resF41 <- summary_block1(res4 = resI4[[3]],  num.nb, data2 = resI4[[1]])
save(resF41, file=paste("M4", case, n.blocks, "-", num.nb, "resF41.Rdata", sep=''))

load(paste(getwd(),"/M5irregular60-2resI5.Rdata", sep = ""))
resF51 <- summary_block1(res4 = resI5[[3]],  num.nb, data2 = resI5[[1]] )
save(resF51, file=paste("M5", case, n.blocks, "-", num.nb, "resF51.Rdata", sep = ''))

load(paste(getwd(),"/M6irregular60-2resI6.Rdata", sep = ""))
resF61 <- summary_block1(res4 = resI6[[3]],  num.nb, data2 = resI6[[1]])
save(resF61, file = paste("M6", case, n.blocks,"-", num.nb, "resF61.Rdata", sep = ''))

rm(resI1, resF11,resI2, resF21,resI3, resF31,
   resI4, resF41,resI5, resF51,resI6, resF61)

#predictions
load(paste(getwd(),"/M1param.posirregular60-2.Rdata", sep = ""))
load(paste(getwd(),"/M1irregular60-2resI1.Rdata", sep = ""))
resF12 <- summary_block2(n.blocks, num.nb, param.pos, res4 = resI1)
save(resF12, file = paste("M1", case, n.blocks,"-", num.nb,"resF12.Rdata",sep = ''))
rm(resI1, resF12, param.pos)

load(paste(getwd(),"/M2param.posirregular60-2.Rdata", sep = ""))
load(paste(getwd(),"/M2irregular60-2resI2.Rdata", sep = ""))
resF22 <- summary_block2( n.blocks,  num.nb, param.pos, res4 = resI2)
save(resF22, file = paste("M2", case, n.blocks, "-", num.nb, "resF22.Rdata", sep = ''))
rm(resI2, resF22, param.pos)

load(paste(getwd(),"/M3param.posirregular60-2.Rdata", sep = ""))
load(paste(getwd(),"/M3irregular60-2resI3.Rdata", sep = ""))
resF23 <- summary_block2( n.blocks,  num.nb, param.pos, res4 = resI3)
save(resF23, file = paste("M3", case, n.blocks, "-", num.nb, "resF32.Rdata", sep = ''))
rm(resI3, resF23, param.pos)


load(paste(getwd(),"/M4param.posirregular60-2.Rdata", sep = ""))
load(paste(getwd(),"/M4irregular60-2resI4.Rdata", sep = ""))
resF24 <- summary_block2( n.blocks,  num.nb, param.pos, res4 = resI4)
save(resF24, file = paste("M4", case, n.blocks, "-", num.nb, "resF42.Rdata",sep = ''))
rm(resI4, resF24, param.pos)

load(paste(getwd(),"/M5param.posirregular60-2.Rdata", sep = ""))
load(paste(getwd(),"/M5irregular60-2resI5.Rdata", sep = ""))
resF25 <- summary_block2( n.blocks,  num.nb, param.pos, res4 = resI5)
save(resF25, file = paste("M5", case, n.blocks, "-", num.nb, "resF52.Rdata",sep = ''))
rm(resI5, resF25, param.pos)

load(paste(getwd(), "/M6param.posirregular60-2.Rdata", sep = ""))
load(paste(getwd(), "/M6irregular60-2resI6.Rdata", sep = ""))
resF26 <- summary_block2( n.blocks,  num.nb, param.pos, res4 = resI6)
save(resF26, file = paste("M6", case, n.blocks, "-", num.nb, "resF62.Rdata",sep = ''))
rm(resI6, resF26, param.pos)

## open results

load(paste(getwd(),"/M1irregular60-2resI1.Rdata", sep = ""))
load(paste(getwd(),"/M1irregular60-2resF11.Rdata", sep = ""))
load(paste(getwd(),"/M1irregular60-2resF12.Rdata", sep = ""))

waic1 <- resI1[[3]]$waic$waic
lpml1 <- sum(log(resI1[[3]]$cpo$cpo))
time1 <- resI1[[3]]$cpu.used[4]
rmse1 <- sum(resF11[[2]])/3
rmsp1 <- sum(resF12[[2]])/3
res1  <- c(waic1, lpml1, rmse1,rmsp1, time1)

load(paste(getwd(),"/M2irregular60-2resI2.Rdata", sep = ""))
load(paste(getwd(),"/M2irregular60-2resF21.Rdata", sep = ""))
load(paste(getwd(),"/M2irregular60-2resF22.Rdata", sep = ""))

waic2 <- resI2[[3]]$waic$waic
lpml2 <- sum(log(resI2[[3]]$cpo$cpo))
time2 <- resI2[[3]]$cpu.used[4]
rmse2 <- sum(resF21[[2]])/3
rmsp2 <- sum(resF22[[2]])/3
res2  <- c(waic2, lpml2, rmse2, rmsp2, time2)

load(paste(getwd(),"/M3irregular60-2resI3.Rdata", sep = ""))
load(paste(getwd(),"/M3irregular60-2resF31.Rdata", sep = ""))
load(paste(getwd(),"/M3irregular60-2resF32.Rdata", sep = ""))
waic3 <- resI3[[3]]$waic$waic
lpml3 <- sum(log(resI3[[3]]$cpo$cpo))
time3 <- resI3[[3]]$cpu.used[4]
rmse3 <- sum(resF31[[2]])/3
rmsp3 <- sum(resF23[[2]])/3
res3  <- c(waic3, lpml3, rmse3,rmsp3, time3)

load(paste(getwd(),"/M4irregular60-2resI4.Rdata", sep = ""))
load(paste(getwd(),"/M4irregular60-2resF41.Rdata", sep = ""))
load(paste(getwd(),"/M4irregular60-2resF42.Rdata", sep = ""))
waic4 <- resI4[[3]]$waic$waic
lpml4 <- sum(log(resI4[[3]]$cpo$cpo))
time4 <- resI4[[3]]$cpu.used[4]
rmse4 <- sum(resF41[[2]])/3
rmsp4 <- sum(resF24[[2]])/3
res4  <- c(waic4, lpml4, rmse4,rmsp4, time4)

load(paste(getwd(),"/M5irregular60-2resI5.Rdata", sep = ""))
load(paste(getwd(),"/M5irregular60-2resF51.Rdata", sep = ""))
load(paste(getwd(),"/M5irregular60-2resF52.Rdata", sep = ""))
waic5 <- resI5[[3]]$waic$waic
lpml5 <- sum(log(resI5[[3]]$cpo$cpo))
time5 <- resI5[[3]]$cpu.used[4]
rmse5 <- sum(resF51[[2]])/3
rmsp5 <- sum(resF25[[2]])/3
res5  <- c(waic5, lpml5, rmse5,rmsp5, time5)

load(paste(getwd(),"/M6irregular60-2resI6.Rdata", sep = ""))
load(paste(getwd(),"/M6irregular60-2resF61.Rdata", sep = ""))
load(paste(getwd(),"/M6irregular60-2resF62.Rdata", sep = ""))
waic6 <- resI6[[3]]$waic$waic
lpml6 <- sum(log(resI6[[3]]$cpo$cpo))
time6 <- resI6[[3]]$cpu.used[4]
rmse6 <- sum(resF61[[2]])/3
rmsp6 <- sum(resF26[[2]])/3
res6  <- c(waic6, lpml6, rmse6, rmsp6, time6)

resfinal <- data.frame(rbind(res1, res2, res3, res4, res5, res6))
names(resfinal) <- c("WAIC", "LPML", "MSE", "MSP", "TIME")

print(resfinal)

# if best model 
resF61[[1]]
# show predictions 
saida.pred <- resF26[[1]]


