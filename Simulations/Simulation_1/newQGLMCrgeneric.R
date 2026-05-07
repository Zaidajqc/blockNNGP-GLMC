#' @name inla.rgeneric.LMCk2_QblockNNGP.model.pc.prior
#'
#' @description blockNNGP-GLMC
#'
## @author:  Zaida Quiroz (\email{zquiroz@pucp.edu.pe}).

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

  #  interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
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


  # Q function defines the precision matrix which is defined in a similar way of W. Here we define the precision matrix of the blockNNGP latent effect.
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
    Qblock1<-  PrecblockNNGP(n, n.blocks,C1,nb,ind_obs1,num1,indb)
    Qnewblock1 <- kronecker(Qblock1,(gamma1inv)%*%t(gamma1inv))

    gamma2inv <- matrix(invgamma[,2], k,1)

    R2 <-   exp(-phi2 *coords.D)
    diag(R2) <- 1
    C2 	 <- sigmasq2 * R2
    Qblock2<-  PrecblockNNGP(n, n.blocks,C2,nb,ind_obs1,num1,indb)
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
    param <- interpret.theta()
    res <- log(lam1) + log(lam2) - (2*(log(param$range1))) - (lam1/param$range1) - (lam2*param$sigma1)+
      log(lam1) + log(lam2) - (2*(log(param$range2))) - (lam1/param$range2) - (lam2*param$sigma2)+
      dnorm(param$lambda, 0, 10, log = TRUE) 
      return(res)
  }

  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    # The Initial values form a diagonal matrix
     return(c(initial.sigma, initial.range,initial.sigma, initial.range, 0.1))
  }

  # A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }

  res <- do.call(match.arg(cmd), args = list())
  return(res)
}


