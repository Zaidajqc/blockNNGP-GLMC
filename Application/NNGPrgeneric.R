# R implementation of NNGP latent effect with rgeneric

# In order to define the NNGP latent effect in INLA, we need to define:
# A variable theta is defined by INLA in the code to store theta=(theta1, theta2)  to provide an internal representation of the hyperparameters (sigma2 and  phi, respectively) to make numerical optimization easier.
# The mean of the latent effects: mu
# The precision of the latent effects: Q(theta)
# A ‘graph’, with a binary representation of the precision matrix: W
# The initial values of the parameters.
# A log-normalizing constant.
# The log-prior of theta.

'inla.rgeneric.NNGP.model' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  # interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
  interpret.theta <- function() {

    return(
      list(sigmasq = exp(-theta[1L]),
           phi = b - (b-a)/ (1 + exp(theta[2L])))
    )
  }
  
  # the graph function represents the entries of the precision matrix that are non-zero. W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    require(Matrix)
    
    return(Diagonal(nrow(W), x = 1) + W)
  }
  
  # meancov_nn function computes the matrices Bk and Fk for each observation k.
  meancov_nn1 = function(i,loc,AdjMatrix,Sigma){
    ind_neig <- which(AdjMatrix[,i]==1)
    invC_nbi <- chol2inv(chol(Sigma[ind_neig,ind_neig])) # inverse of Fbi
    B_bi     <- (Sigma[i,ind_neig] )%*%invC_nbi
    F_bi     <- Sigma[i,i] -(B_bi%*%Sigma[ind_neig,i])
    Bstar_bi <- matrix(0,1,dim(loc)[1])
    Bstar_bi[,ind_neig] <- -B_bi 
    Bstar_bi[i] <- 1 
    val         <- list(B_bi=B_bi, F_bi =F_bi,Bstar_bi=Bstar_bi)
    return(val)
  }
  
  
  # PrecNNGP function computes the precision matrix of NNGP
  Prec_NNGP  = function(loc,AdjMatrix,Sigma){
    
    nloc          <- dim(loc)[1]  
    var           <- matrix(NA,nloc,1)
    Fs_1          <- matrix(0,nloc,nloc)
    Bb            <- matrix(0,nloc,nloc)
    
    for (j in 1:nloc){
      #print(j)
      if(j==1){
        var[j]    <- Sigma[j,j] 
        Fs_1[j,j] <- 1/var[j] 
        Bb[j,j]   <- 1
      }
      if(j>1){
        res       <- meancov_nn1(j,loc,AdjMatrix,Sigma)
        var[j]    <- res$F_bi
        Fs_1[j,j] <- 1/var[j] 
        Bb[j,]    <- res$Bstar_bi
      }
    }
    
    m1 <- as(Bb , "dgCMatrix")
    m2 <- as(Fs_1 , "dgCMatrix")
    invCs        <- crossprod(m1,m2)%*%m1
    
    return(invCs)
  }
  
  
  # Q function defines the precision matrix which is defined in a similar way of W. Here we define the precision matrix of the NNGP latent effect.
  Q <- function() {
    require(Matrix)
    
    param <- interpret.theta()
    R <- exp(-param$phi*coords.D)
    diag(R) <- 1
    C <- param$sigmasq*R

    
    return(  Prec_NNGP(coords.D,AdjMatrix,C) )
  }
  
  # mu function is the mean of the blockNNGP latent effect which is zero.
  mu <- function()
  {
    return(numeric(0))
  }
  
  # log.norm.const function computes the normalizing constant, INLA computes it if numeric is zero.
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  # log.prior function computes the pdf of prior distributions for sigmasq and phi. In particular, for the marginal variance we set a gamma distribution  with parameters 1 and 0.00005,and for phi=2/range we set a uniform distribution on (a,b), where a and b are associated to the minimum and maximum distance between locations, in this case a=1 and b=30. extra terms that appear in the definition of the log-density of the prior are due to the change of variable involved. INLA works with (theta1, theta2)  internally, but the prior is set on (sigmasq, phi).
  log.prior <- function() {
    param <- interpret.theta()
    res   <- dgamma(param$sigmasq, 1, 5e-05, log = TRUE) + log(param$sigmasq) +
      (2*log(1/(b-a))) + log(param$phi-a) + log(b - param$phi) 
    return(res)
  }
  
  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    return(c(0, 0))
  }
  
  # A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}


'inla.rgeneric.NNGP.model.pc.prior' <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
    "log.prior", "quit"),
  theta = NULL) {

# interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
  interpret.theta <- function() {

    return(
      list(sigma = exp(theta[1L]),
           range = exp(theta[2L]))
    )
  }

# the graph function represents the entries of the precision matrix that are non-zero. W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    require(Matrix)

    return(Diagonal(nrow(W), x = 1) + W)
  }

# meancov_nn function computes the matrices Bk and Fk for each observation k.
meancov_nn1 = function(i,loc,AdjMatrix,Sigma){
  ind_neig <- which(AdjMatrix[,i]==1)
  invC_nbi <- chol2inv(chol(Sigma[ind_neig,ind_neig])) # inverse of Fbi
  B_bi     <- (Sigma[i,ind_neig] )%*%invC_nbi
  F_bi     <- Sigma[i,i] -(B_bi%*%Sigma[ind_neig,i])
  Bstar_bi <- matrix(0,1,dim(loc)[1])
  Bstar_bi[,ind_neig] <- -B_bi 
  Bstar_bi[i] <- 1 
  val         <- list(B_bi=B_bi, F_bi =F_bi,Bstar_bi=Bstar_bi)
  return(val)
}


# PrecNNGP function computes the precision matrix of NNGP
Prec_NNGP  = function(loc,AdjMatrix,Sigma){

  nloc          <- dim(loc)[1]  
  var           <- matrix(NA,nloc,1)
  Fs_1          <- matrix(0,nloc,nloc)
  Bb            <- matrix(0,nloc,nloc)
  
  for (j in 1:nloc){
    #print(j)
    if(j==1){
      var[j]    <- Sigma[j,j] 
      Fs_1[j,j] <- 1/var[j] 
      Bb[j,j]   <- 1
    }
    if(j>1){
      res       <- meancov_nn1(j,loc,AdjMatrix,Sigma)
      var[j]    <- res$F_bi
      Fs_1[j,j] <- 1/var[j] 
      Bb[j,]    <- res$Bstar_bi
    }
  }
  
   m1     <- as(Bb , "dgCMatrix")
   m2     <- as(Fs_1 , "dgCMatrix")
   invCs  <- crossprod(m1,m2)%*%m1

  return(invCs)
}


# Q function defines the precision matrix which is defined in a similar way of W. Here we define the precision matrix of the NNGP latent effect.
  Q <- function() {
    require(Matrix)

    param   <- interpret.theta()
    
    sigmasq <- param$sigma^2
    phi     <- 2/param$range
    
    R       <- exp(-phi*coords.D)
    diag(R) <- 1
    C       <- sigmasq*R

    return(  Prec_NNGP(coords.D,AdjMatrix,C) )
  }

  # mu function is the mean of the blockNNGP latent effect which is zero.
 mu <- function()
  {
    return(numeric(0))
  }

  # log.norm.const function computes the normalizing constant, INLA computes it if numeric is zero.
 log.norm.const <- function() {
    return(numeric(0))

  }

 # log.prior function computes the pdf of prior distributions for sigmasq and phi. 
 #In particular, pc.prior for range and sigma
 #INLA works with (theta1, theta2)  internally, but the prior is set on (sigmasq, phi).
 log.prior <- function() {
   param  <- interpret.theta()
   res    <- log(lam1) + log(lam2) -2*log(param$range)- (lam1/(param$range)) - lam2*((param$sigma))
   return(res)
 }

  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    return(c(initial.sigma, initial.range))
  }

# A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }

  res <- do.call(match.arg(cmd), args = list())
  return(res)
}


'inla.rgeneric.NNGP.model.matern' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  # interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
  interpret.theta <- function() {
    return(
      list(sigmasq = exp(-theta[1L]),
           phi = b - (b-a)/ (1 + exp(theta[2L])))
    )
  }
  
  # the graph function represents the entries of the precision matrix that are non-zero. W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    require(Matrix)
    
    return(Diagonal(nrow(W), x = 1) + W)
  }
  
  # meancov_nn function computes the matrices Bk and Fk for each observation k.
  meancov_nn1 = function(i,loc,AdjMatrix,Sigma){
    ind_neig <- which(AdjMatrix[,i]==1)
    invC_nbi <- chol2inv(chol(Sigma[ind_neig,ind_neig])) # inverse of Fbi
    B_bi     <- (Sigma[i,ind_neig] )%*%invC_nbi
    F_bi     <- Sigma[i,i] -(B_bi%*%Sigma[ind_neig,i])
    Bstar_bi <- matrix(0,1,dim(loc)[1])
    Bstar_bi[,ind_neig] <- -B_bi 
    Bstar_bi[i] <- 1 
    val         <- list(B_bi=B_bi, F_bi =F_bi,Bstar_bi=Bstar_bi)
    return(val)
  }
  
  
  # PrecNNGP function computes the precision matrix of NNGP
  Prec_NNGP  = function(loc,AdjMatrix,Sigma){
    
    nloc          <- dim(loc)[1]  
    var           <- matrix(NA,nloc,1)
    Fs_1          <- matrix(0,nloc,nloc)
    Bb            <- matrix(0,nloc,nloc)
    
    for (j in 1:nloc){
      #print(j)
      if(j==1){
        var[j]    <- Sigma[j,j] 
        Fs_1[j,j] <- 1/var[j] 
        Bb[j,j]    <- 1
      }
      if(j>1){
        res       <- meancov_nn1(j,loc,AdjMatrix,Sigma)
        var[j]    <- res$F_bi
        Fs_1[j,j] <- 1/var[j] 
        Bb[j,]    <- res$Bstar_bi
      }
    }
    
    m1 <- as(Bb , "dgCMatrix")
    m2 <- as(Fs_1 , "dgCMatrix")
    invCs        <- crossprod(m1,m2)%*%m1
    
    return(invCs)
  }
  
  
  # Q function defines the precision matrix which is defined in a similar way of W. Here we define the precision matrix of the NNGP latent effect.
  Q <- function() {
    require(Matrix)
    
    param   <- interpret.theta()
    nu      <- 0.5
    R       <- (coords.D*param$phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=coords.D*param$phi, nu=nu)
    diag(R) <- 1
    C       <- param$sigmasq*R
    
    
    return(Prec_NNGP(coords.D,AdjMatrix,C) )
  }
  
  # mu function is the mean of the blockNNGP latent effect which is zero.
  mu <- function()
  {
    return(numeric(0))
  }
  
  # log.norm.const function computes the normalizing constant, INLA computes it if numeric is zero.
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  # log.prior function computes the pdf of prior distributions for sigmasq and phi. In particular, for the marginal variance we set a gamma distribution  with parameters 1 and 0.00005,and for phi=2/range we set a uniform distribution on (a,b), where a and b are associated to the minimum and maximum distance between locations, in this case a=1 and b=30. extra terms that appear in the definition of the log-density of the prior are due to the change of variable involved. INLA works with (theta1, theta2)  internally, but the prior is set on (sigmasq, phi).
  log.prior <- function() {
    param <- interpret.theta()
    res   <- dgamma(param$sigmasq, 1, 5e-05, log = TRUE) + log(param$sigmasq) +
      (2*log(1/(b-a))) + log(param$phi-a) + log(b - param$phi) 
    return(res)
  }
  
  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    return(c(0, 0))
  }
  
  # A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}



'inla.rgeneric.NNGP.model.pc.prior.matern' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  # interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
  interpret.theta <- function() {
    return(
      list(sigma = exp(theta[1L]),
           range = exp(theta[2L]))
    )
  }
  
  # the graph function represents the entries of the precision matrix that are non-zero. W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    require(Matrix)
    
    return(Diagonal(nrow(W), x = 1) + W)
  }
  
  # meancov_nn function computes the matrices Bk and Fk for each observation k.
  meancov_nn1 = function(i,loc,AdjMatrix,Sigma){
    ind_neig <- which(AdjMatrix[,i]==1)
    invC_nbi <- chol2inv(chol(Sigma[ind_neig,ind_neig])) # inverse of Fbi
    B_bi     <- (Sigma[i,ind_neig] )%*%invC_nbi
    F_bi     <- Sigma[i,i] -(B_bi%*%Sigma[ind_neig,i])
    Bstar_bi <- matrix(0,1,dim(loc)[1])
    Bstar_bi[,ind_neig] <- -B_bi 
    Bstar_bi[i] <- 1 
    val         <- list(B_bi=B_bi, F_bi =F_bi,Bstar_bi=Bstar_bi)
    return(val)
  }
  
  
  # PrecNNGP function computes the precision matrix of NNGP
  Prec_NNGP  = function(loc,AdjMatrix,Sigma){
    
    nloc          <- dim(loc)[1]  
    var           <- matrix(NA,nloc,1)
    Fs_1          <- matrix(0,nloc,nloc)
    Bb            <- matrix(0,nloc,nloc)
    
    for (j in 1:nloc){
      if(j==1){
        var[j]    <- Sigma[j,j] 
        Fs_1[j,j] <- 1/var[j] 
        Bb[j,j]   <- 1
      }
      if(j>1){
        res       <- meancov_nn1(j,loc,AdjMatrix,Sigma)
        var[j]    <- res$F_bi
        Fs_1[j,j] <- 1/var[j] 
        Bb[j,]    <- res$Bstar_bi
      }
    }
    
    m1    <- as(Bb , "dgCMatrix")
    m2    <- as(Fs_1 , "dgCMatrix")
    invCs <- crossprod(m1,m2)%*%m1
    
    return(invCs)
  }
  
  
  # Q function defines the precision matrix which is defined in a similar way of W. Here we define the precision matrix of the NNGP latent effect.
  Q <- function() {
    require(Matrix)
    
    param   <- interpret.theta()
    
    sigmasq <- param$sigma^2
    phi     <- 2/param$range
    
    nu      <- 0.5
    R       <- exp(-phi*coords.D)
    diag(R) <- 1
    C       <- sigmasq*R

    return(  Prec_NNGP(coords.D,AdjMatrix,C) )
  }
  
  # mu function is the mean of the blockNNGP latent effect which is zero.
  mu <- function()
  {
    return(numeric(0))
  }
  
  # log.norm.const function computes the normalizing constant, INLA computes it if numeric is zero.
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  # log.prior function computes the pdf of prior distributions for sigmasq and phi. 
  #In particular, pc.prior for range and sigma
  #INLA works with (theta1, theta2)  internally, but the prior is set on (sigmasq, phi).
  log.prior <- function() {
    param <- interpret.theta()
    res   <- log(lam1) + log(lam2) -2*log(param$range)- (lam1/(param$range)) - lam2*((param$sigma))
    return(res)
  }
  
  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    return(c(initial.sigma, initial.range))
  }
  
  # A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

