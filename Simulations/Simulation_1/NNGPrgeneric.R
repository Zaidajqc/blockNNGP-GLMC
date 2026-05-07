# R implementation of NNGP latent effect for NNGP-GLM model
## with rgeneric


# Define previous variables as global to avoid warnings()
utils::globalVariables(c("k", "W", "lam1", "lam2",
                         "n", "coords.D" ,  "a" ,"b",
                         "initial.range", "initial.sigma"))


'inla.rgeneric.NNGP.model.pc.prior2' <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
    "log.prior", "quit"),
  theta = NULL) {

# interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
  interpret.theta <- function() {
    return(
      list(sigma1 = exp(theta[1L]),
           range1 = exp(theta[2L]),
           sigma2 = exp(theta[3L]),
           range2 = exp(theta[4L]),
           lambda = (theta[5L]) )
    )
  }

# the graph function represents the entries of the precision matrix that are non-zero. W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    return( W)
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

    param <- interpret.theta()
    
    sigmasq1 <- param$sigma1^2
    phi1<- 2/param$range1
    sigmasq2 <- param$sigma2^2
    phi2<- 2/param$range2
    
    M <-  diag(1, 2)
    M[upper.tri(M)] <- param$lambda
    invgamma <- solve(M)
    
    gamma1inv <- matrix(invgamma[,1], k,1)
    
    R1 <-   exp(-phi1 *coords.D)
    diag(R1) <- 1
    C1 	 <-  sigmasq1* R1
    Qblock1<-  Prec_NNGP(coords.D,AdjMatrix,C1) 
   Qnewblock1 <- kronecker(Qblock1,(gamma1inv)%*%t(gamma1inv))
    
    gamma2inv <- matrix(invgamma[,2], k,1)
    
    R2 <-   exp(-phi2 *coords.D)
    diag(R2) <- 1
    C2 	 <- sigmasq2 * R2
    Qblock2<-   Prec_NNGP(coords.D,AdjMatrix,C2) 
    Qnewblock2 <- kronecker(Qblock2,(gamma2inv)%*%t(gamma2inv))
    
    Q <- Qnewblock1 + Qnewblock2
    
    return( Q)
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
   param = interpret.theta()
   res <- log(lam1) + log(lam2) - (2*(log(param$range1))) - (lam1/param$range1) - (lam2*param$sigma1)+
     log(lam1) + log(lam2) - (2*(log(param$range2))) - (lam1/param$range2) - (lam2*param$sigma2)+
     dnorm(param$lambda, 0, 10, log = TRUE) 
   return(res)
 }

  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    return(c(initial.sigma, initial.range,initial.sigma, initial.range, 0.1))
  }

# A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }

  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

