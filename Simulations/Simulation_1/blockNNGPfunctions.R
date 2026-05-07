## functions to run blockNNGP-GLMC and NNGP-GLMC models.
## @author:  Zaida Quiroz (\email{zquiroz@pucp.edu.pe}).

## Description:
##
## blockNNGP_reg() function  fits  block-NNGP models through REGULAR BLOCKS using INLA
# Arguments:
## case: name to identify regular blockNNGP results
## loc: locations
## y: observed data
## X: covariates
## w: spatial random effects (used to compare possterior estimations with the true values.
## n.blocks: number of  blocks  (regular)
## num.nb: number of neighbor blocks.
## dir.save: directory to save the results
## formula:
##  ‘y ~ 1 + x + f(idx, model = blockNNGP.model)’
##
## For different family distribution chage the inla() function on line 162.
## For a list of possible alternatives and use ‘inla.doc’ for detailed
##  docs for individual families.
##
##  Value:
##      ‘inla’ returns an object of class ‘"inla"’.
##


blockNNGP_reg = function(case= 'regular', data.est, n.blocks, num.nb, name.cov,par.cov=NULL,data.pred){

  all.data <- rbind(data.est,data.pred)
  data.est.orig <- data.est 
  data.est <- all.data
  loc=data.est[,1:2]
  
 

###%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## The next codes implement the Adjacency matrix of blockNNGP and objects that the
## 'inla.rgeneric.blockNNGP.model' function needs, for specific number of
## blocks ( n.blocks) and num.nb neighbor blocks.
## The user does not need to change it.

n.partition <- sqrt(n.blocks)

nloc <- dim(loc)[1]
minlocx<- min(loc[,1])
maxlocx<- max(loc[,1])
minlocy<- min(loc[,2])
maxlocy<- max(loc[,2]) 
hb <- seq(minlocx,maxlocx,length.out=n.partition+1)
vb <- seq(minlocy,maxlocy,length.out=n.partition+1)
#hb <- seq(min(loc[,1]),max(loc[,1]),length.out=n.partition+1)
#vb <- seq(min(loc[,2]),max(loc[,2]),length.out=n.partition+1)

blocks <- matrix(NA, ncol=1, nrow=dim(loc)[1])
loc.blocks <- matrix(NA, ncol=2, nrow=n.blocks)
#file1 <- paste('/case',case,'_',nloc,'_',n.blocks,'_',num.nb,sep="")


#png(paste(dir.save,'/case',case,'_',nloc,'_',n.blocks,"_NNGPblocks1.png",sep=""))

plot(loc,  main= paste('M = ',n.blocks,', nb=',num.nb))
abline(h = hb[2:n.partition],col='red',pch=19)
abline(v = vb[2:n.partition],col='red',pch=19)


k <- 1
for (j in 2:(n.partition+1)){
for (i in 2:(n.partition+1)){
liminfh <- hb[i-1]
limsuph <- hb[i]
liminfv <- vb[j-1]
limsupv <- vb[j]
indblock <- which(loc [,1] >= liminfh & loc[,1] <= limsuph & loc[,2] >= liminfv & loc[,2] <= limsupv)
if(length(indblock)>0){
  blocks[indblock] <- k
  points(loc[blocks==k,1:2],col=k)
  loc.blocks[k,1] <- (liminfh + limsuph)/2
  loc.blocks[k,2] <- (liminfv + limsupv)/2
  text(loc.blocks[k,1],loc.blocks[k,2],k,col='red')
  k <- k + 1
  }
}
}

#dev.off()

n.blocks <- k-1 
loc.blocks <- loc.blocks[1:n.blocks,]

ind <- sort.int(loc.blocks[,2], index.return=TRUE)

x2 <- ind$x
indexsort <- ind$ix
x1 <- loc.blocks[indexsort,1]


sortloc <- cbind(x1, x2)
dist.mat <- rdist(sortloc)


AdjMatrix <-  matrix(0, n.blocks, n.blocks)
AdjMatrix[1,1] <-0


for (j in 2:n.blocks){
if (j <= num.nb+1) AdjMatrix[1:(j-1),j] = 1
if (j > num.nb+1){
 ind1 <-  (sort(dist.mat[,j], index.return=TRUE))$ix
 ind <- (ind1[which(ind1<j)])[1:num.nb]
 AdjMatrix[ind,j] = 1
}
}

#AdjMatrix


#g1 = graph.adjacency(AdjMatrix, mode='directed')
#is.dag(g1)

nloc <- dim(data.est.orig)[1]

blocks.pred <- blocks[(nloc+1):dim(all.data)[1]]  
data.pred$blocks <- blocks.pred

blocks <- blocks[1:nloc]
#data.est <- data.est.orig
loc <- loc[1:nloc,]


## to make things simple reorder loc respect to the blocks
ind1 <- sort.int(blocks, index.return=TRUE)
newdata <- data.est.orig[(ind1$ix),]
loc <- loc[(ind1$ix),]
blocks <- blocks[(ind1$ix)]
newdata$blocks <- blocks
##rnorm_n.obs <- rnorm_n.obs[(ind1$ix)]
#y <- y[(ind1$ix)]
#if(length(data.est$w)>0){
#  w <- w[(ind1$ix)]
#}
#X <- X[(ind1$ix),]

meanloc=mean(table(blocks))
print(paste("approx points per block is" , meanloc))


### needed indexes to built the precision matrix of block_NNGP
  newindex     <- NULL
  nb           <- matrix(NA,n.blocks,1)
  nb[1]        <- length(which(blocks==1))
  for (j in 1:n.blocks){
    ind_obs      <- which(blocks==j)
    newindex     <- c(newindex,ind_obs)
    if(j>1){
    nbj=length(ind_obs)
    nb[j] <- nb[j-1]+nbj
    }
  }
  nloc         <- dim(loc)[1]
  ind_obs1      <- which(blocks==1)
  num1          <- seq(1:length(ind_obs1))


indb <- NULL
for (k in 1:(n.blocks-1)){
indb[[k]] <- util.index(k+1,blocks,AdjMatrix,newindex)
}


## mask for precision-blockNNGP
coords.D 	<- rdist(loc)
if(name.cov==""){name.cov="exponential"}
if(name.cov=="exponential"){
  C1 <-  exp(-0.04*coords.D)
}
if(name.cov=="matern"){
  nu<- par.cov
  C1 <-  (coords.D*0.04)^nu/(2^(nu-1)*gamma(nu))*besselK(x=coords.D*0.04, nu=nu)
diag(C1) <- 1
}
invC   <-  PrecblockNNGP(nloc, n.blocks,C1,nb,ind_obs1,num1,indb)
invCsp <- as.matrix(invC)
invCsp[which(invC>0)] <- 1
invCsp[which(invC<0)] <- 1
invCsp[which(invC==0)] <- 0

W = invCsp
W <- as(W, "sparseMatrix")

###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Xnew <- X[,1:(dim(X)[2])]
#colnames(Xnew) <- 1:(dim(X)[2])
  
# response variable and covariates
data1<- newdata
#if(length(data.est$w)>0){
#  data1<- data.frame(y=y,x=Xnew,w,loc1=loc[,1],loc2=loc[,2])
#}else{
#  data1<- data.frame(y=y,x=Xnew,loc1=loc[,1],loc2=loc[,2])
#}
#data1$idx <- 1:nrow(data1)


datablock <-list(data1=data1, W = W, n= nloc, n.blocks= n.blocks,nb =nb,
            ind_obs1=ind_obs1,num1=num1,indb=indb,coords.D=coords.D,nu=par.cov,data.pred)


return(datablock)
}



blockNNGP_IRREG = function(case="irregular", data.est , n.blocks, num.nb,name.cov,par.cov=NULL,data.pred){
  

 data.est.orig <- data.est 
  data.pred.orig <- data.pred
  all.data <- rbind(data.est,data.pred)
  data.est <- all.data
  loc=data.est[,1:2]

  
  ###%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## The next codes implement the Adjacency matrix of blockNNGP and objects that the
  ## 'inla.rgeneric.blockNNGP.model' function needs, for specific number of
  ## blocks ( n.blocks) and num.nb neighbor blocks.
  ## The user does not need to change it.
  
  nloc <- dim(loc)[1]
  
  
  #file1 <- paste('/case',case,'_',nloc,'_',n.blocks,'_',num.nb,sep="")
 
  ds <-  data.frame(x = loc[ ,1], y = loc[ ,2])  
  fit <- kmeans(ds, centers = n.blocks)
  ds <- cbind(ds, cluster = as.factor(fit$cluster))

  # Call function ---------------------------------------------------------------
  #png(paste(dir.save,'/case',case,'_',nloc,'_',n.blocks,"_NNGPblocks1.png",sep=""))
  VoronoiPlotly(fit, ds, n.sd.x = 3, n.sd.y = 3, print.ggplot = F)
  #dev.off()
  
  
  
  center.b <- fit$centers
  
  par(mfrow=c(1,2))
  plot(center.b)
  text(center.b[,1]+0.05,center.b[,2],rownames(center.b),col='red')
  
  
  #dist.blocks <- rdist(center.b)
  
  ind.block <- sort.int(center.b[,2], index.return=TRUE)$ix
  
  newvindex <- cbind( ind.block,1:n.blocks)
  
  ordered_polygons <- center.b[ind.block, ]
  
  plot(ordered_polygons)
  text(ordered_polygons[,1]+0.05,ordered_polygons[,2],1:n.blocks,col='red')
  rownames(ordered_polygons)<- 1:n.blocks
  
  ds1 <- ds
  for( l in (1:n.blocks)){
    indv <- which(ds$cluster==l)
    indv1 <- newvindex[which(newvindex[,1]==l),2]
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
      AdjMatrix[ind,j] = 1
    }
  }
  
  
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
  nb           <- matrix(NA,n.blocks,1)
  nb[1]        <- length(which(blocks==1))
  for (j in 1:n.blocks){
    ind_obs      <- which(blocks==j)
    newindex     <- c(newindex,ind_obs)
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
    indb[[k]] <- util.index(k+1,blocks,AdjMatrix,newindex)
  }
  
  
  ## mask for precision-blockNNGP
  coords.D 	<- rdist(loc2)
  if(name.cov==""){name.cov="exponential"}
  if(name.cov=="exponential"){
    C1 <-  exp(-0.04*coords.D)
    diag(C1) <- 1
    }
  if(name.cov=="matern"){
    nu<- par.cov
    C1 <-  (coords.D*0.04)^nu/(2^(nu-1)*gamma(nu))*besselK(x=coords.D*0.04, nu=nu)
    diag(C1) <- 1
  }
  invC   <-  PrecblockNNGP(nloc2, n.blocks,C1,nb,ind_obs1,num1,indb)
  invCsp <- as.matrix(invC)
  invCsp[which(invC>0)] <- 1
  invCsp[which(invC<0)] <- 1
  invCsp[which(invC==0)] <- 0
  
  W = invCsp
  W <- as(W, "sparseMatrix")
  
  ###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  data1 <- newdata
    
  datablock <-list(data1=data1, W = W, n= nloc2, n.blocks= n.blocks,nb =nb,
                   ind_obs1=ind_obs1,num1=num1,indb=indb,coords.D=coords.D, nu=par.cov,data.pred)
  
  
  return(datablock)
}





## Description:
##
## blockNNGP_reg() function  fits  NNGP models  using INLA
# Arguments:
## case: name to identify NNGP results
## loc: locations
## y: observed data
## X: covariates
## w: spatial random effects (used to compare possterior estimations with the true values.
## num.nb: number of neighbor observations.
## dir.save: directory to save the results.
## formula:
##  ‘y ~ 1 + x + f(idx, model = NNGP.model)’
##
## For different family distribution chage the inla() function on line 162.
## For a list of possible alternatives and use ‘inla.doc’ for detailed
##  docs for individual families.
##
##  Value:
##      ‘inla’ returns an object of class ‘"inla"’.
##


NNGP = function(case="NNGP", data.est, n.blocks =1, num.nb, name.cov,par.cov=NULL,data.pred ){
  
  
  
  ###%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## The next codes implement the Adjacency matrix of NNGP and objects that the
  ## 'inla.rgeneric.NNGP.model' function needs, for specific number of
  ## neigbors (num.nb).
  ## The user does not need to change it.
  
  n.blocks <- 1
  loc=data.est[,1:2]  
  nloc <- dim(loc)[1]
  
  
  ind <- sort.int(loc[,2], index.return=TRUE)
  indexsort <- ind$ix
  data.est <- data.est[indexsort,]
  
  
  sortloc <- data.est[,1:2]
  dist.mat <- rdist(sortloc)
  
  AdjMatrix <-  matrix(0, nloc, nloc)
  AdjMatrix[1,1] <-0
  
  
  for (j in 2:nloc){
    if (j <= num.nb+1) AdjMatrix[1:(j-1),j] = 1
    if (j > num.nb+1){
      ind1 <-  (sort(dist.mat[,j], index.return=TRUE))$ix
      ind <- (ind1[which(ind1<j)])[1:num.nb] 
      AdjMatrix[ind,j] = 1
    }
  }
  
  
  
  
  
  ## mask for precision-NNGP
  C1 <- exp(-0.04*dist.mat)
  invC   <-  Prec_NNGP(sortloc,AdjMatrix,C1) 
  invCsp <- as.matrix(invC)
  invCsp[which(invC>0)] <- 1
  invCsp[which(invC<0)] <- 1
  invCsp[which(invC==0)] <- 0
  
  W = invCsp
  W <- as(W, "sparseMatrix")
  
   
  datablock <-list(data1=data.est, W = W,coords.D=dist.mat, AdjMatrix=AdjMatrix,data.pred)
  
  return(datablock)
}


# This function creates the structure of neighbor blocks 
# reordering of data
# for blockNNGP regular, block-NNGP irregular and NNGP

blockNNGP.struct <- function(case,data.est, n.blocks, num.nb,name.cov,par.cov,data.pred){
  if(case== 'regular'){
    build.struct <- blockNNGP_reg(case= 'regular', data.est, n.blocks, num.nb,name.cov,par.cov,data.pred)
  }else  if (case== 'irregular'){
      build.struct <- blockNNGP_IRREG(case="irregular", data.est, n.blocks, num.nb,name.cov,par.cov,data.pred)
#      build.struct <- blockNNGP_IRREG (case="irregular", loc,  y, X, w.obs , n.partition, num.nb,name.cov)
  }else if(case=="NNGP"){
    build.struct <- NNGP (case="NNGP", data.est, n.blocks=1, num.nb, name.cov,par.cov,data.pred )
  } else{
    print("Please define in case: regular, irregular or NNGP")
  }
return(build.struct)  
}
  


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




inla.blockNNGP.model <- function(case, name.cov,name.prior,datafill,prior.set) {
  inla.rgeneric.blockNNGP = choose.model(case,
                                         name.cov,
                                         name.prior)
  if(name.prior=="prior1"&& case!= 'NNGP'&& name.cov=="exponential"){
  ## set the blockNNGP as latent effect
  INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                             n= datafill[[3]], n.blocks= datafill[[4]],nb =datafill[[5]],
                             ind_obs1=datafill[[6]],num1=datafill[[7]],indb=datafill[[8]],
                             coords.D=datafill[[9]], a=prior.set[1],b=prior.set[2])
  }else if(name.prior=="prior1"&& case!= 'NNGP'&& name.cov=="matern"){
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               n= datafill[[3]], n.blocks= datafill[[4]],nb =datafill[[5]],
                               ind_obs1=datafill[[6]],num1=datafill[[7]],indb=datafill[[8]],
                               coords.D=datafill[[9]], a=prior.set[1],b=prior.set[2],nu=datafill[[10]])
    }else if(name.prior=="pc.prior" && case!= 'NNGP'&& name.cov=="exponential"){
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1

    ## set the blockNNGP as latent effect
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               n= datafill[[3]], n.blocks= datafill[[4]],nb =datafill[[5]],
                               ind_obs1=datafill[[6]],num1=datafill[[7]],indb=datafill[[8]],
                               coords.D=datafill[[9]], 
                               lam1=lam11, lam2=lam22,
                               initial.range=initial.range,
                               initial.sigma=initial.sigma)
  } else if(name.prior=="pc.prior" && case!= 'NNGP'&& name.cov=="matern"){
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1
    
    ## set the blockNNGP as latent effect
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               n= datafill[[3]], n.blocks= datafill[[4]],nb =datafill[[5]],
                               ind_obs1=datafill[[6]],num1=datafill[[7]],indb=datafill[[8]],
                               coords.D=datafill[[9]], 
                               lam1=lam11, lam2=lam22,
                               initial.range=initial.range,
                               initial.sigma=initial.sigma, nu=datafill[[10]])
  } else if(name.prior=="prior1"&& case== 'NNGP'&& name.cov=="exponential"){
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D= datafill[[3]], AdjMatrix= datafill[[4]],
                               a=prior.set[1],b=prior.set[2])
  }else if(name.prior=="prior1"&& case== 'NNGP'&& name.cov=="matern"){
    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D= datafill[[3]], AdjMatrix= datafill[[4]],
                               a=prior.set[1],b=prior.set[2])
  }else if(name.prior=="pc.prior"&& case== 'NNGP'&& name.cov=="exponential"){ 
    prior.range <- c(prior.set[1:2])
    prior.sigma <- c(prior.set[3:4])
    lam11 <- -log(prior.range[2])*prior.range[1]
    lam22 <- -log(prior.sigma[2])/prior.sigma[1]
    
    initial.range <- log(prior.range[1]) + 1
    initial.sigma <- log(prior.sigma[1]) - 1

    INLA::inla.rgeneric.define(inla.rgeneric.blockNNGP, W = datafill[[2]],
                               coords.D= datafill[[3]], AdjMatrix= datafill[[4]], 
                               lam1=lam11, lam2=lam22,
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
                               coords.D= datafill[[3]], AdjMatrix= datafill[[4]], 
                               lam1=lam11, lam2=lam22,
                               initial.range=initial.range,
                               initial.sigma=initial.sigma)
  }
}





###############################
## 2. function for predictions
##############################


blockNNGP_predLMC = function(case, n.blocks, num.nb, data.est, pred.data,res,n.sample,family="gaussian" ){
  
 
  Y.pred <- pred.data$y1
  n.pred = length(Y.pred)
  w.pred = pred.data$w1
  X.pred = cbind(1,pred.data$X.2)
  Y2.pred = pred.data$y2
  w2.pred = pred.data$w2
  X2.pred = cbind(1,pred.data$X.2)
  pred.coords= cbind(pred.data$loc.1,pred.data$loc.2) 
  
  blocks.pred = pred.data$blocks
  
  ## to make things simple reorder loc respect to the blocks
  ind1 <- sort.int(blocks.pred, index.return=TRUE)
  
  Y.pred <- Y.pred[(ind1$ix)]
  blocks.pred <- blocks.pred[(ind1$ix)]
  X.pred <- X.pred[(ind1$ix),]
  w.pred <- w.pred[(ind1$ix)]
  pred.coords <- pred.coords[(ind1$ix),]
  Y2.pred <- Y2.pred[(ind1$ix)]
  X2.pred <- X2.pred[(ind1$ix),]
  w2.pred <- w2.pred[(ind1$ix)]
  
  # posterior samples
  set.seed(1)
  param.pos <-   inla.posterior.sample( n.sample, res)
  
  ## computing ypred
  n <- length(data.est$y)/2
  k <- 2  # number of response variables
  obs.coords <- cbind(data.est$loc.1, data.est$loc.2)
  indx <- seq(1, k*n, k)
  obs.coords <- obs.coords[indx,]
  blocks <- data.est$blocks[indx]
  mean.y1 <- NULL
  mean.y2 <- NULL
  for(i in 1:n.pred){
    print(i)
    Y.predi <- Y.pred[i]
    blocks.predi <- blocks.pred[i]
    blocki <- which(blocks==blocks.predi)
    if(length(blocki)>1){
      obs.D <- rdist(obs.coords[blocki,])
    }else{
      obs.D <- 0
    }
    obs.pred.D <- rdist( rbind(pred.coords[i,], obs.coords[blocki,]))[1,]
    obs.pred.D <- obs.pred.D[-1]
    y0 <- matrix(NA, n.sample, 1)
    y2 <- matrix(NA, n.sample, 1)
    for(s in 1:n.sample){
      n1 <- seq((2*n+1),(4*n),2)
      n2 <- seq((2*n+2),(4*n),2)
      #n1 <- (2*n+1):(3*n)
      #n2 <- (3*n+1):(4*n)
      ncoef <- (4*n+1):(4*n+4)
      if(family=="gaussian"){
       postng.tausq<- sqrt(1/param.pos[[s]]$hyperpar[1])
        postng.tausq2<- sqrt(1/param.pos[[s]]$hyperpar[3])
        postng.sigmasq <- (exp(param.pos[[s]]$hyperpar[4]))^2
        postng.phi <- 2/exp(param.pos[[s]]$hyperpar[2])
        postng.w1 <- param.pos[[s]]$latent[n1,]
        postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[5]))^2
        postng.phi2 <- 2/exp(param.pos[[s]]$hyperpar[6])
        postng.w2 <- param.pos[[s]]$latent[n2,]
        postng.lambda <- param.pos[[s]]$hyperpar[7]
      }else if(family=="poisson"){
        postng.sigmasq <- (exp(param.pos[[s]]$hyperpar[1]))^2
        postng.phi <- 2/exp(param.pos[[s]]$hyperpar[2])
        postng.v1 <- param.pos[[s]]$latent[n1,]
        postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[3]))^2
        postng.phi2 <- 2/exp(param.pos[[s]]$hyperpar[4])
        postng.v2 <- param.pos[[s]]$latent[n2,]
        postng.lambda <- param.pos[[s]]$hyperpar[5]
      }
          postng.beta<- c(param.pos[[s]]$latent[ncoef[1]], param.pos[[s]]$latent[ncoef[3]])
          postng.beta<-as.matrix(postng.beta, 2, 1)
          postng.beta2<- c(param.pos[[s]]$latent[ncoef[2]], param.pos[[s]]$latent[ncoef[4]])
          postng.beta2<-as.matrix(postng.beta2, 2, 1)
          
       # w0 for k=1
      postng.w2 <- postng.v2  - postng.lambda*postng.v1 
      C1  <- postng.sigmasq*exp(-postng.phi*obs.pred.D)
      Cn0 <- postng.sigmasq*exp(-postng.phi*obs.D)
      invCn0 <- solve(Cn0)
      m  <- t(C1)%*%(invCn0 %*%postng.v1[blocki])
      v <- postng.sigmasq - t(C1)%*%(invCn0 %*%C1)
      w0 <- m + sqrt(v)*rnorm(1)
      XX <- t(as.matrix(X.pred[i,],1,2))
      if(family=="gaussian"){
      y0[s] <- XX%*%postng.beta + w0 + postng.tausq*rnorm(1)
      }else if (family=="poisson"){
        y0[s] <- rpois(1, exp(XX%*%postng.beta + w0 ))
     }
      # w0 for k=2
      C12  <- postng.sigmasq2*exp(-postng.phi2*obs.pred.D)
      Cn02 <- postng.sigmasq2*exp(-postng.phi2*obs.D)
      invCn02 <- solve(Cn02)
      m2  <- t(C12)%*%(invCn02 %*%postng.w2[blocki])
      v2 <- postng.sigmasq2 - t(C12)%*%(invCn02 %*%C12)
      w02 <- m2 + sqrt(v2)*rnorm(1)
      XX <- t(as.matrix(X2.pred[i,],1,2))
      if(family=="gaussian"){
        y2[s] <- XX%*%postng.beta2 + postng.lambda*w0 + w02 + postng.tausq2*rnorm(1)
      }else if(family=="poisson"){
        y2[s] <- rpois(1,exp(XX%*%postng.beta2 + postng.lambda*w0 + w02) )
      }
    }
    
    mean.y1   <- c(mean.y1, mean(y0))
    mean.y2   <- c(mean.y2, mean(y2))
  }
  
  val <- data.frame(Y.pred, Y2.pred, mean.y1, mean.y2)
  
  return(val)
}


NNGP_predLMC = function(case="NNGP", num.nb, data.est, pred.data,res,n.sample,family="gaussian" ){
  
  #pred.data <- list(Y.pred, w.pred, X.pred, Y2.pred,w2.pred, X2.pred, pred.coords)
  
  Y.pred <- pred.data$y1
  n.pred = length(Y.pred)
  w.pred = pred.data$w1
  X.pred = cbind(1,pred.data$X.2)
  Y2.pred = pred.data$y2
  w2.pred = pred.data$w2
  X2.pred = cbind(1,pred.data$X.2)
  pred.coords= cbind(pred.data$loc.1,pred.data$loc.2) 
  
  # posterior samples
  set.seed(1)
  param.pos <-   inla.posterior.sample( n.sample, res)
  
  ## computing ypred
  obs.coords <- cbind(data.est$loc.1, data.est$loc.2)
  k=2
  n = length(data.est$y)/2
  indx <- seq(1, k*n, k)
  obs.coords <- obs.coords[indx,]
  blocks <- data.est$blocks[indx]

  #  blocks <- data.est$blocks
  mean.y1 <- NULL
  mean.y2 <- NULL
  for(i in 1:n.pred){
    print(i)
    Y.predi <- Y.pred[i]
    library(FNN)
    blocki <- (get.knn(data.frame(x=c(pred.coords[i,1],as.numeric(obs.coords[,1])), y=c(pred.coords[i,2],as.numeric(obs.coords[,2]))), num.nb)$nn.index)[1,]-1
    
    if(length(blocki)>1){
      obs.D <- rdist(obs.coords[blocki,])
    }else{
      obs.D <- 0
    }
    obs.pred.D <- rdist( rbind(pred.coords[i,], obs.coords[blocki,]))[1,]
    obs.pred.D <- obs.pred.D[-1]
    y0 <- matrix(NA, n.sample, 1)
    y2 <- matrix(NA, n.sample, 1)
    for(s in 1:n.sample){
      n1 <- seq((2*n+1),(4*n),2)
      n2 <- seq((2*n+2),(4*n),2)
      ncoef <- (4*n+1):(4*n+4)
      if(family=="gaussian"){
      postng.tausq<- sqrt(1/param.pos[[s]]$hyperpar[1])
      postng.tausq2<- sqrt(1/param.pos[[s]]$hyperpar[2])
      postng.sigmasq <- (exp(param.pos[[s]]$hyperpar[3]))^2
      postng.phi <- 2/exp(param.pos[[s]]$hyperpar[4])
      postng.w1 <- param.pos[[s]]$latent[n1,]
      postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[5]))^2
      postng.phi2 <- 2/exp(param.pos[[s]]$hyperpar[6])
      postng.w2 <- param.pos[[s]]$latent[n2,]
      postng.lambda <- param.pos[[s]]$hyperpar[7]
      }else if(family=="poisson"){
        postng.sigmasq <- (exp(param.pos[[s]]$hyperpar[1]))^2
        postng.phi <- 2/exp(param.pos[[s]]$hyperpar[2])
        postng.v1 <- param.pos[[s]]$latent[n1,]
        postng.sigmasq2 <- (exp(param.pos[[s]]$hyperpar[3]))^2
        postng.phi2 <- 2/exp(param.pos[[s]]$hyperpar[4])
        postng.v2 <- param.pos[[s]]$latent[n2,]
        postng.lambda <- param.pos[[s]]$hyperpar[5]
      }
      postng.beta<- c(param.pos[[s]]$latent[ncoef[1]], param.pos[[s]]$latent[ncoef[3]])
      postng.beta<-as.matrix(postng.beta, 2, 1)
      # k=2
      postng.beta2<- c(param.pos[[s]]$latent[ncoef[2]], param.pos[[s]]$latent[ncoef[4]])
      postng.beta2<-as.matrix(postng.beta2, 2, 1)
      postng.w2 <- postng.v2  - postng.lambda*postng.v1 
      # w0 for k=1
      C1  <- postng.sigmasq*exp(-postng.phi*obs.pred.D)
      Cn0 <- postng.sigmasq*exp(-postng.phi*obs.D)
      invCn0 <- solve(Cn0)
      m  <- t(C1)%*%(invCn0 %*%postng.v1[blocki])
      v <- postng.sigmasq - t(C1)%*%(invCn0 %*%C1)
      w0 <- m + sqrt(v)*rnorm(1)
      XX <- t(as.matrix(X.pred[i,],1,2))
      if(family=="gaussian"){
      y0[s] <- XX%*%postng.beta + w0 + postng.tausq*rnorm(1)
      }else if(family=="poisson"){
        y0[s] <- rpois(1,exp(XX%*%postng.beta + w0 ))
      }
      # w0 for k=2
      C12  <- postng.sigmasq2*exp(-postng.phi2*obs.pred.D)
      Cn02 <- postng.sigmasq2*exp(-postng.phi2*obs.D)
      invCn02 <- solve(Cn02)
      m2  <- t(C12)%*%(invCn02 %*%postng.w2[blocki])
      v2 <- postng.sigmasq2 - t(C12)%*%(invCn02 %*%C12)
      w02 <- m2 + sqrt(v2)*rnorm(1)
      XX <- t(as.matrix(X2.pred[i,],1,2))
      if(family=="gaussian"){
      y2[s] <- XX%*%postng.beta2 + postng.lambda*w0 + w02 + postng.tausq2*rnorm(1)
      }else if(family=="poisson"){
        y2[s] <- rpois(1,exp(XX%*%postng.beta2 + postng.lambda*w0 + w02 ))
      }
      }
    mean.y1   <- c(mean.y1, mean(y0))
    mean.y2   <- c(mean.y2, mean(y2))
  }
  
  val <- data.frame(Y.pred, Y2.pred, mean.y1, mean.y2)
  
  return(val)  
}
 
