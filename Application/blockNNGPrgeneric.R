## Description:
##
# R implementation of blockNNGP-GLMC latent effect with rgeneric

# A variable theta is defined by INLA in the code to store theta=(theta1, theta2)  to provide an internal representation of the hyperparameters (sigma2 and  phi, respectively) to make numerical optimization easier.
# In order to define the block-NNGP latent effect in INLA, we need to define the next functions:
# The mean of the latent effects: mu
# The precision of the latent effects: Q(theta)
# A ‘graph’, with a binary representation of the precision matrix: W
# The initial values of the parameters.
# A log-normalizing constant.
# The log-prior of theta.

## Arguments:
##
#blockNNGP.model <- inla.rgeneric.define(inla.rgeneric.blockNNGP.model, W = W, n= n, n.blocks= n.blocks,nb =nb,ind_obs1=ind_obs1,num1=num1,indb=indb,coords.D=coords.D)


# Define previous variables as global to avoid warnings()
utils::globalVariables(c("k", "W", "lam1", "lam2",
                         "n" , "n.blocks" ,"nb" ,
                         "ind_obs1" ,"num1" ,"indb" ,
                         "coords.D" , "a" ,"b" ,
                         "initial.range", "initial.sigma"))


'inla.rgeneric.LMCk2_QblockNNGP.model.pc.prior' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  # interpret.theta function will take the parameters in the internal scale and return the marginal 
  # variance and phi parameters:
  interpret.theta <- function() {
    return(
      list(sigma = sapply(theta[as.integer(1:k)], function(x) { exp(x) }),
           range = sapply(theta[as.integer((k+1):(2*k))], function(x) { exp(x) }),
           lambda =  theta[-as.integer(1:(2*k))] )
    )
  }
  
  # the graph function represents the entries of the precision matrix that are non-zero. 
  #W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    
    require(Matrix)
    
    return( W)
  }
  
  
  # meancov_nn function computes the matrices Bbk and Fbk for each block k.
  meancov_nn <- function(Sigma, ind_obs, ind_neigblocks, indnum){
    invC_nbi <- chol2inv(chol( Sigma[ind_neigblocks, ind_neigblocks] ))
    B_bi     <- (Sigma[ind_obs, ind_neigblocks])%*%invC_nbi
    F_bi     <- Sigma[ind_obs, ind_obs] -(B_bi%*%Sigma[ind_neigblocks, ind_obs])
    invFbi   <- chol2inv(chol(F_bi))
    Bstar_bi <- matrix(0,length(ind_obs), dim(Sigma)[1])
    
    Bstar_bi[,ind_neigblocks] <- -B_bi
    Bstar_bi[indnum]       <- 1
    
    val     <- list( invFbi = invFbi, Bstar_bi = Bstar_bi)
    return(val)
  }
  
  # PrecNNGP function computes the precision matrix of blockNNGP
  PrecblockNNGP  <- function(nloc, n.blocks, Sigma, nb, ind_obs1, num1, indb){
    
    Fs_1         <- matrix(0,nloc,nloc)
    Bb           <- matrix(0,nloc,nloc)
    
    Fs_1[1:nb[1],1:nb[1]]   <- chol2inv(chol(Sigma[ind_obs1,ind_obs1]))
    
    Bstar_bi                 <- matrix(0,length(ind_obs1),nloc)
    diag(Bstar_bi[num1,num1]) <- 1
    Bb[1:nb[1],]             <- Bstar_bi
    
    for (j in 2:n.blocks){
      ress <- meancov_nn(Sigma,indb[[j-1]][[1]],indb[[j-1]][[2]],indb[[j-1]][[4]])
      Fs_1[(nb[j-1]+1):(nb[j]),(nb[j-1]+1):(nb[j])] <- ress$invFbi
      Bb[(nb[j-1]+1):(nb[j]),] <- ress$Bstar_bi
    }
    
    Bbb <- as(Bb , "dgCMatrix")
    Fs_11 <- as(Fs_1 , "dgCMatrix")
    
    invCs        <- crossprod(Bbb,Fs_11)%*%Bbb
    
    return(invCs)
  }
  
  
  
  # Q function defines the precision matrix which is defined in a similar way of W. 
  # Here we define the precision matrix of the GLMC-blockNNGP latent effect.
  Q <- function() {
    require(Matrix)
    
    param <- interpret.theta()
    
    sigmasq <-  (param$sigma)^2 
    phi <- 2/(param$range) 
    
    M <-  diag(1, k)
    M[upper.tri(M)] <- param$lambda #fill by rows 
    invgamma <- solve(M)
 

 
    Q <- 0
    for(j in 1:k){
      gamma1inv <- matrix(invgamma[,j], k,1)
      Cj 	 <- sigmasq[j] * exp(-phi[j]*coords.D) 
      Qblockj<-  PrecblockNNGP(n, n.blocks,Cj,nb,ind_obs1,num1,indb)
      rm(Cj)
      Qnewblockj <- kronecker(Qblockj,(gamma1inv)%*%t(gamma1inv))
      rm(Qblockj)
      Q <- Q + Qnewblockj
      rm(Qnewblockj)
    }
   
    return( Q)
  }
  
  # mu function is the mean of the blockNNGP latent effect which is zero.
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  #In particular, pc.prior for range and sigma
  #INLA works with (theta1, theta2)  internally, but the prior is set on (sigmasq, phi).
  log.prior <- function() {
    param <- interpret.theta()
    # log-Prior for  lambdaij parameters 
    res <- sum(dnorm(param$lambda, 0, 10, log = TRUE)) +
      k*(log(lam1) + log(lam2)) - (2*sum(log(param$range))) - sum(lam1/param$range) - (lam2*sum(param$sigma)) 
    return(res)
  }
  
  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    # The Initial values form a diagonal matrix
    return ( c( rep(0.1, k),  rep(0.1, k), rep(0.1,k*(k-1)/2)))
  }
  
  # A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}


