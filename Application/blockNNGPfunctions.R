## Functions to run blockNNGP-GLMC 
## @author:  Zaida Quiroz (\email{zquiroz@pucp.edu.pe}) and Carlos Gonzáles.

#' @title blockNNGP_IRREG()
#' @description function  to order the data, create the mask for Q, 
#' and create blocks and indexed required for  blockNNGP-GLMC models
#' through  Voronoi polygons.

blockNNGP_IRREG = function(case = "irregular", data.est , n.blocks, 
                           num.nb, name.cov, par.cov = NULL, data.pred){
  

  loc <-data.est[,1:2]
  if(length(which(duplicated(loc))==T)>0)  {
    print("duplicated locations...") 
    inddup   <- which(duplicated(loc)==TRUE)
    data.est <- data.est[-inddup,]
  }
    locp <- data.pred[,1:2]
    if(length(which(duplicated(locp))==T)>0)  {
      print("duplicated locations...") 
      inddup    <- which(duplicated(locp)==TRUE)
      data.pred <- data.pred[-inddup,]
    }
    
    data.est.orig  <- data.est 
    data.pred.orig <- data.pred
    all.data <- rbind(data.est,data.pred)
    data.est <- all.data
    loc <-data.est[,1:2]
  
 
  ###%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## The next codes implement the Adjacency matrix of blockNNGP and objects that the
  ## 'inla.rgeneric.blockNNGP.model' function needs, for specific number of
  ## blocks ( n.blocks) and num.nb neighbor blocks.
  ## The user does not need to change it.
  
  nloc <- dim(loc)[1]
  
  
  ds <-  data.frame(x = loc[ ,1], y = loc[ ,2])  
  fit <- kmeans(ds, centers = n.blocks)
  ds <- cbind(ds, cluster = as.factor(fit$cluster))
  
 
  # Call function ---------------------------------------------------------------
  VoronoiPlotly(fit, ds, n.sd.x = 3, n.sd.y = 3, print.ggplot = F)

  
  center.b <- fit$centers
  
  par(mfrow=c(1,2))
  plot(center.b)
  text(center.b[,1]+0.05, center.b[,2], rownames(center.b), col='red')
  
  ind.block <- sort.int(center.b[,2], index.return=TRUE)$ix
  
  newvindex <- cbind( ind.block, 1:n.blocks)
  
  ordered_polygons <- center.b[ind.block, ]
  
  plot(ordered_polygons)
  text(ordered_polygons[,1]+0.05, ordered_polygons[,2], 1:n.blocks, col='red')
  rownames(ordered_polygons) <- 1:n.blocks
  
  ds1 <- ds
  for( l in (1:n.blocks)){
    indv  <- which(ds$cluster==l)
    indv1 <- newvindex[which(newvindex[,1]==l), 2]
    ds1$cluster[indv] <- indv1
  }
  
  blocks <-  ds1$cluster
  loc.blocks <- ordered_polygons
  

  dist.mat <- rdist(loc.blocks)
  
  AdjMatrix <-  matrix(0, n.blocks, n.blocks)
  AdjMatrix[1,1] <-0
  
  for (j in 2:n.blocks){
    if (j <= num.nb+1) AdjMatrix[1:(j-1),j] = 1
    if (j > num.nb+1){
      ind1 <-  (sort.int(dist.mat[,j], index.return=TRUE))$ix
      ind <- (ind1[which(ind1<j)])[1:num.nb]
      AdjMatrix[ind,j] <-  1
    }
  }
  
  #AdjMatrix
  
  
  nloc1 <- dim(data.est.orig)[1]
  
  blocks.pred <- blocks[(nloc1+1):dim(all.data)[1]]  
  data.pred$blocks <- blocks.pred
  
  blocks <- blocks[1:nloc1]
  data.est <- data.est.orig
  loc1 <- loc[1:nloc1,]
  
  ## to make things simple reorder loc respect to the blocks
  ind1 <- sort.int(blocks, index.return=TRUE)
  newdata <- data.est[(ind1$ix),]
  loc2 <- loc1[(ind1$ix),]
  blocks <- blocks[(ind1$ix)]
  newdata$blocks <- blocks

  meanloc=mean(table(blocks))
  print(paste("approx points per block is" , meanloc))
  
  
  ### needed indexes to built the precision matrix of block_NNGP
  newindex     <- NULL
  nb           <- matrix(NA,n.blocks, 1)
  nb[1]        <- length(which(blocks==1))
  for (j in 1:n.blocks){
    ind_obs      <- which(blocks==j)
    newindex     <- c(newindex, ind_obs)
    if(j>1){
      nbj=length(ind_obs)
      nb[j] <- nb[j-1]+nbj
    }
  }
  nloc2         <- dim(loc2)[1]
  ind_obs1      <- which(blocks==1)
  num1          <- seq(1:length(ind_obs1))
  
  
  indb <- NULL
  for (k in 1:(n.blocks-1)){
    indb[[k]] <- util.index(k+1,blocks, AdjMatrix,newindex)
  }
  
  
  ## mask for precision-blockNNGP
  coords.D 	<- rdist(loc2)
  if(name.cov==""){name.cov="exponential"}
  if(name.cov=="exponential"){
    C1 <-  exp(-0.04*coords.D)
    diag(C1) <- 1
    }
  if(name.cov=="matern"){
    nu <- par.cov
    C1 <-  (coords.D*0.04)^nu/(2^(nu-1)*gamma(nu))*besselK(x=coords.D*0.04, nu=nu)
    diag(C1) <- 1
  }
  invC   <-  PrecblockNNGP(nloc2, n.blocks, C1, nb, ind_obs1, num1, indb)
  invCsp <- as.matrix(invC)
  invCsp[which(invC>0)]  <- 1
  invCsp[which(invC<0)]  <- 1
  invCsp[which(invC==0)] <- 0
  
  W <- invCsp
  W <- as(W, "sparseMatrix")
  
  ###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  data1 <- newdata
 
  datablock <-list(data1=data1, W = W, n= nloc2, n.blocks= n.blocks, nb =nb,
                   ind_obs1=ind_obs1, num1=num1, indb=indb, coords.D=coords.D, 
                   nu=par.cov, data.pred, AdjMatrix=AdjMatrix)
  
  
  return(datablock)
}

#' @title blockNNGP.struct()
#' @description This function calls the blockNNGP_IRREG() function.

blockNNGP.struct <- function(case,data.est, n.blocks, num.nb,name.cov, par.cov, data.pred){
  if(case== 'regular'){
    build.struct <- blockNNGP_reg(case= 'regular', data.est, n.blocks, num.nb,name.cov,par.cov,data.pred)
  }else  if (case== 'irregular'){
      build.struct <- blockNNGP_IRREG(case="irregular", data.est, n.blocks, num.nb,name.cov, par.cov,data.pred)
  }else if(case=="NNGP"){
    build.struct <- NNGP (case="NNGP", data.est, n.blocks=1, num.nb, name.cov, par.cov, data.pred )
  } else{
    print("Please define in case: regular, irregular or NNGP")
  }
return(build.struct)  
}
  
#' @title data_proc()
#' @description  function to process (order) the data and
#' create block and indexes for  blockNNGP-GLMC rgeneric
#' it calls blockNNGP.struct() function.

data_proc = function(all.data, n.blocks, num.nb, k){
  data.est  <- all.data[[1]]
  data.pred <- all.data[[2]] 
  
  ##########################################
  #### Run irregular voronoi block-NNGP models 
  ##########################################
  
  case='irregular'
  
  name.cov="matern"
  
  datafill <- blockNNGP.struct(case= 'irregular', 
                               data.est, 
                               n.blocks, num.nb, name.cov, par.cov = 0.5, data.pred)
  return(datafill)
}


#' @title choose.model()
#' @description This function is called by inla.blockNNGP.model.

choose.model <- function(case,name.cov,name.prior){
  if(case== 'regular'&& name.cov=="exponential" && name.prior=="prior1"){
    model.rgeneric = inla.rgeneric.blockNNGP.model
  } else   if(case== 'regular'&& name.cov=="matern" && name.prior=="prior1"){
    model.rgeneric = inla.rgeneric.blockNNGP.model.matern
  } else   if(case== 'irregular'&& name.cov=="exponential" && name.prior=="prior1"){
    model.rgeneric = inla.rgeneric.blockNNGP.model
  } else   if(case== 'irregular'&& name.cov=="matern" && name.prior=="prior1"){
    model.rgeneric = inla.rgeneric.blockNNGP.model
  } else   if(case== 'NNGP'&& name.cov=="exponential" && name.prior=="prior1"){
    model.rgeneric = inla.rgeneric.NNGP.model
  }else  if(case== 'NNGP'&& name.cov=="matern" && name.prior=="prior1"){
    model.rgeneric = inla.rgeneric.NNGP.model.matern
  }else  if(case== 'regular'&& name.cov=="exponential" && name.prior=="pc.prior"){
    model.rgeneric = inla.rgeneric.blockNNGP.model.pc.prior
  }else   if(case== 'regular'&& name.cov=="matern" && name.prior=="pc.prior"){
    model.rgeneric = inla.rgeneric.blockNNGP.model.pc.prior
  }else  if(case== 'irregular'&& name.cov=="exponential" && name.prior=="pc.prior"){
    model.rgeneric = inla.rgeneric.blockNNGP.model.pc.prior
  }else  if(case== 'irregular'&& name.cov=="matern" && name.prior=="pc.prior"){
    model.rgeneric = inla.rgeneric.blockNNGP.model.pc.prior
  }else  if(case== 'NNGP'&& name.cov=="exponential" && name.prior=="pc.prior"){
    model.rgeneric = inla.rgeneric.NNGP.model.pc.prior
  }else  if(case== 'NNGP'&& name.cov=="matern" && name.prior=="pc.prior"){
    model.rgeneric = inla.rgeneric.NNGP.model.pc.prior.matern
  }else{
    print("set case, name.cov and name.prior")
  }
  return(model.rgeneric)
  
}


#' @title inla.blockNNGP.model()
#' @description runs inla.rgeneric function to fit blockNNGP and NNGP models using INLA

inla.blockNNGP.model <- function(case, name.cov, name.prior, datafill, prior.set) {
  inla.rgeneric.blockNNGP = choose.model(case,
                                         name.cov,
                                         name.prior)
  if(name.prior=="prior1"&& case!= 'NNGP'&& name.cov=="exponential"){
  ## set the blockNNGP as latent effect
  INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                             n = datafill[[3]], n.blocks = datafill[[4]], nb = datafill[[5]],
                             ind_obs1 = datafill[[6]],num1 = datafill[[7]], indb = datafill[[8]],
                             coords.D = datafill[[9]], a = prior.set[1], b = prior.set[2])
  }else if(name.prior=="prior1"&& case!= 'NNGP'&& name.cov=="matern"){
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               n = datafill[[3]], n.blocks = datafill[[4]], nb = datafill[[5]],
                               ind_obs1 = datafill[[6]], num1 = datafill[[7]], indb = datafill[[8]],
                               coords.D = datafill[[9]], a = prior.set[1], b = prior.set[2], nu = datafill[[10]])
    }else if(name.prior=="pc.prior" && case!= 'NNGP'&& name.cov=="exponential"){
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1

    ## set the blockNNGP as latent effect
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               n = datafill[[3]], n.blocks = datafill[[4]], nb = datafill[[5]],
                               ind_obs1 = datafill[[6]], num1 = datafill[[7]], indb = datafill[[8]],
                               coords.D = datafill[[9]], 
                               lam1 = lam11, lam2 = lam22,
                               initial.range = initial.range,
                               initial.sigma = initial.sigma)
  } else if(name.prior=="pc.prior" && case!= 'NNGP'&& name.cov=="matern"){
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1
    
    ## set the blockNNGP as latent effect
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               n = datafill[[3]], n.blocks = datafill[[4]], nb = datafill[[5]],
                               ind_obs1 = datafill[[6]], num1 = datafill[[7]], indb = datafill[[8]],
                               coords.D = datafill[[9]], 
                               lam1 = lam11, lam2 = lam22,
                               initial.range = initial.range,
                               initial.sigma = initial.sigma, nu = datafill[[10]])
  } else if(name.prior=="prior1"&& case== 'NNGP'&& name.cov=="exponential"){
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D = datafill[[3]], AdjMatrix = datafill[[4]],
                               a = prior.set[1], b = prior.set[2])
  }else if(name.prior=="prior1"&& case== 'NNGP'&& name.cov=="matern"){
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D = datafill[[3]], AdjMatrix = datafill[[4]],
                               a=prior.set[1],b=prior.set[2])
  }else if(name.prior=="pc.prior"&& case== 'NNGP'&& name.cov=="exponential"){ 
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1

    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D = datafill[[3]], AdjMatrix = datafill[[4]], 
                               lam1 = lam11, lam2=lam22,
                               initial.range=initial.range,
                               initial.sigma=initial.sigma)
  }else { #pc.prior, NNGP, matern
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1
    
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D = datafill[[3]], AdjMatrix = datafill[[4]], 
                               lam1 = lam11, lam2 = lam22,
                               initial.range = initial.range,
                               initial.sigma = initial.sigma)
  }
}





###############################
## 2. function for predictions
##############################


#' @title blockNNGP_pred()
#' @description This function computes predictions for blockNNGP models.

blockNNGP_predLMC = function(case, n.blocks, num.nb, data.est, pred.data, n.sample,
                             family = "gaussian", k = 2 , param.pos, AdjMatrix){
 
  Y.pred <- pred.data$y.1
  n.pred <- length(Y.pred)
  X.pred <- cbind(1, pred.data$loc.2, pred.data$loc.2^2)
  Y2.pred <- pred.data$y.2
  X2.pred <- cbind(1, pred.data$loc.2, pred.data$loc.2^2)
  pred.coords <- cbind(pred.data$loc.1, pred.data$loc.2) 
  if(k==3){
    Y3.pred <- pred.data$y.3
    X3.pred <- cbind(1, pred.data$loc.2, pred.data$loc.2^2)
  }
  blocks.pred <- pred.data$blocks
  
  ## to make things simple reorder loc respect to the blocks
  ind1 <- sort.int(blocks.pred, index.return=TRUE)
  
  Y.pred      <- Y.pred[(ind1$ix)]
  blocks.pred <- blocks.pred[(ind1$ix)]
  X.pred      <- X.pred[(ind1$ix),]
  pred.coords <- pred.coords[(ind1$ix),]
  Y2.pred <- Y2.pred[(ind1$ix)]
  X2.pred <- X2.pred[(ind1$ix),]
  if(k==3){
    Y3.pred <- Y3.pred[(ind1$ix)]
    X3.pred <- X3.pred[(ind1$ix),]
  }
  set.seed(1)

  ## computing ypred
  obs.coords <- cbind(data.est$loc.1, data.est$loc.2)
  n <- length(data.est$y)/k
  indx <- seq(1, k*n, k)
  obs.coords <- obs.coords[indx,]
  blocks <- data.est$blocks[indx]
  mean.y1 <- NULL
  mean.y2 <- NULL
  mean.y3 <- NULL
  for(i in 1:n.pred){
    print(i)
     blocks.predi <- blocks.pred[i]
     blocki0 <- which(blocks==blocks.predi)
     if(as.numeric(blocks.predi)==1) blocki <- blocki0
    if(as.numeric(blocks.predi)>1){
         neigblocksi <- which(AdjMatrix[,blocks.predi]==1)
         if(length(neigblocksi)==1) {
           blocki1 <- which(blocks==neigblocksi)
           blocki  <- sort(c(blocki0, blocki1))
         }else{ #now for 2 neighbor blocks
           blocki1 <- which(blocks==neigblocksi[1])
           blocki2 <- which(blocks==neigblocksi[2])
           blocki  <-  sort(c(blocki0, blocki1, blocki2))
         }
    }
    if(length(blocki)>1){
      obs.D <- rdist(obs.coords[blocki,])
    }else{
      obs.D <- 0
    }
    obs.pred.D <- rdist( rbind(pred.coords[i,], obs.coords[blocki,]))[1,]
    obs.pred.D <- obs.pred.D[-1]
    y0 <- matrix(NA, n.sample, 1)
    y2 <- matrix(NA, n.sample, 1)
    y3 <- matrix(NA, n.sample, 1)
    if(k==3)  y3 <- matrix(NA, n.sample, 1)
    for(s in 1:n.sample){
      if(k==2){
        n1 <- seq((2*n+1),(4*n), 2)
        n2 <- seq((2*n+2),(4*n), 2)
        ncoef <- (4*n+1):(4*n+4) # for poisson         
      }
      if(k==3){
        n1 <- seq((k*n+1),(k*2*n), k)
        n2 <- seq((k*n+2),(k*2*n), k)
        n3 <- seq((k*n+3),(k*2*n), k)
         # if  covariates
        ncoef <- (k*2*n+1):(k*2*n+9)
        }
      if(family=="gaussian"){
       postng.tausq     <- sqrt(1/param.pos[[s]]$hyperpar[1])
        postng.tausq2   <- sqrt(1/param.pos[[s]]$hyperpar[3])
        postng.sigmasq  <- (exp(param.pos[[s]]$hyperpar[4]))^2
        postng.phi      <- 2/exp(param.pos[[s]]$hyperpar[2])
        postng.w1       <- param.pos[[s]]$latent[n1,]
        postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[5]))^2
        postng.phi2     <- 2/exp(param.pos[[s]]$hyperpar[6])
        postng.w2       <- param.pos[[s]]$latent[n2,]
        postng.lambda   <- param.pos[[s]]$hyperpar[7]
      }else if(family=="poisson"){
        if(k==2){
        postng.sigmasq  <- (exp(param.pos[[s]]$hyperpar[1]))^2
        postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[2]))^2
        postng.phi      <- 2/exp(param.pos[[s]]$hyperpar[3])
        postng.phi2     <- 2/exp(param.pos[[s]]$hyperpar[4])
        postng.lambda   <- param.pos[[s]]$hyperpar[5]
        postng.v1       <- param.pos[[s]]$latent[n1,]
        postng.v2       <- param.pos[[s]]$latent[n2,]
        postng.beta     <- c(param.pos[[s]]$latent[ncoef[1]], param.pos[[s]]$latent[ncoef[3]])
        postng.beta     <-as.matrix(postng.beta, 2, 1)
        postng.beta2    <- c(param.pos[[s]]$latent[ncoef[2]], param.pos[[s]]$latent[ncoef[4]])
        postng.beta2    <-as.matrix(postng.beta2, 2, 1)
        }else if(k==3){
          postng.sigmasq  <- (exp(param.pos[[s]]$hyperpar[1]))^2
          postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[2]))^2
          postng.sigmasq3 <- (exp(param.pos[[s]]$hyperpar[3]))^2
          postng.phi      <- 2/exp(param.pos[[s]]$hyperpar[4])
          postng.phi2     <- 2/exp(param.pos[[s]]$hyperpar[5])
          postng.phi3     <- 2/exp(param.pos[[s]]$hyperpar[6])
          postng.lambda   <- param.pos[[s]]$hyperpar[7:9]
          postng.v1       <- param.pos[[s]]$latent[n1,]
          postng.v2       <- param.pos[[s]]$latent[n2,]
          postng.v3       <- param.pos[[s]]$latent[n3,]
          # if two covariates
          postng.beta  <- c(param.pos[[s]]$latent[ncoef[1]], param.pos[[s]]$latent[ncoef[4]], param.pos[[s]]$latent[ncoef[7]])
          postng.beta  <- as.matrix(postng.beta, 3, 1)
          postng.beta2 <- c(param.pos[[s]]$latent[ncoef[2]], param.pos[[s]]$latent[ncoef[5]], param.pos[[s]]$latent[ncoef[8]])
          postng.beta2 <- as.matrix(postng.beta2, 3, 1)
          postng.beta3 <- c(param.pos[[s]]$latent[ncoef[3]], param.pos[[s]]$latent[ncoef[6]], param.pos[[s]]$latent[ncoef[9]])
          postng.beta3 <- as.matrix(postng.beta3, 3, 1)
        }
      }
          
      # w0 for k=1
      postng.w2 <- postng.v2  - postng.lambda[1]*postng.v1 
      C1        <- postng.sigmasq*exp(-postng.phi*obs.pred.D)
      Cn0       <- postng.sigmasq*exp(-postng.phi*obs.D)
      invCn0    <- solve(Cn0)
      m  <- t(C1)%*%(invCn0 %*%postng.v1[blocki])
      v  <- postng.sigmasq - t(C1)%*%(invCn0 %*%C1)
      w0 <- m + sqrt(v)*rnorm(1)
      # if one covariate
      XX1 <- t(as.matrix(X.pred[i,], 1, 3))
      if(family=="gaussian"){
      y0[s] <- XX1%*%postng.beta + w0 + postng.tausq*rnorm(1)
      }else if (family=="poisson"){
       y0[s] <- rpois(1, exp(XX1%*%postng.beta + w0 ))
     }
      # w0 for k=2
      C12     <- postng.sigmasq2*exp(-postng.phi2*obs.pred.D)
      Cn02    <- postng.sigmasq2*exp(-postng.phi2*obs.D)
      invCn02 <- solve(Cn02)
      m2  <- t(C12)%*%(invCn02 %*%postng.w2[blocki])
      v2  <- postng.sigmasq2 - t(C12)%*%(invCn02 %*%C12)
      w02 <- m2 + sqrt(v2)*rnorm(1)
 # if one covariate     
      XX2 <- t(as.matrix(X2.pred[i,], 1, 3))
      if(family=="gaussian"){
        y2[s] <- XX2%*%postng.beta2 + postng.lambda[1]*w0 + w02 + postng.tausq2*rnorm(1)
      }else if(family=="poisson"){
# if  covariate 
       y2[s] <- rpois(1, exp(XX2%*%postng.beta2 + postng.lambda[1]*w0 + w02) )
      }
    
    if(k==3){
      # w0 for k=3
      postng.w3 <- postng.v3  - postng.lambda[2]*postng.v1 -  postng.lambda[3]*postng.w2 
      C13       <- postng.sigmasq3*exp(-postng.phi3*obs.pred.D)
      Cn03      <- postng.sigmasq3*exp(-postng.phi3*obs.D)
      invCn03   <- solve(Cn03)
      m3  <- t(C13)%*%(invCn03 %*%postng.w3[blocki])
      v3  <- postng.sigmasq3 - t(C13)%*%(invCn03 %*%C13)
      w03 <- m3 + sqrt(v3)*rnorm(1)
# if  covariate 
     XX3 <- t(as.matrix(X3.pred[i,], 1, 3))
      if(family=="gaussian"){
        y3[s] <- XX3%*%postng.beta3 + postng.lambda[2]*w0 + postng.lambda[3]* w02 + w03 + postng.tausq2*rnorm(1)
      }else if(family=="poisson"){
# if  covariate
        y3[s] <- rpois(1, exp(XX3%*%postng.beta3 + postng.lambda[2]*w0 + postng.lambda[3]* w02 + w03) )
      }
     }
    }
 
    mean.y1   <- c(mean.y1, mean(y0))
    mean.y2   <- c(mean.y2, mean(y2))
    if(k==3)  {
      mean.y3   <- c(mean.y3, mean(y3))
  }
  }
  val <- data.frame(Y.pred, Y2.pred, mean.y1, mean.y2)
  if(k==3)  {
    val <- data.frame(Y.pred, Y2.pred, Y3.pred, mean.y1, mean.y2, mean.y3)
}
  return(val)
}

