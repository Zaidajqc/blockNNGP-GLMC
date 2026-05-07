#' @title summary.blockNNGP()
#' @description This function returns summary of posterior parameters for GLMC-blockNNGP models.
## @author:  Zaida Quiroz (\email{zquiroz@pucp.edu.pe}).

#summary for blockNNGP 
summary.blockNNGP_LMC <- function(name.prior, resf,n.blocks, 
                                  num.nb,prior.set=NULL,family="gaussian"){
  
if(name.prior=="pc.prior")  {
  if (!is.null(resf$summary.hyperpar) && length(resf$summary.hyperpar) > 0) {
   
    if ( length(names(resf$marginals.hyperpar)) > 0) {
      if ( length(resf$marginals.hyperpar$"Precision for the Gaussian observations") > 0) {
        marg.tau1 <- inla.tmarginal(function(x) { 1/x},
                                    resf$marginals.hyperpar$"Precision for the Gaussian observations")
        res.tau1 <- inla.zmarginal(marg.tau1, TRUE)
        res.tau1$mode <- inla.mmarginal(marg.tau1)
        #marg.tau2 <- inla.tmarginal(function(x) { 1/x},
        #                            resf$marginals.hyperpar$"Precision for the Gaussian observations[2]")
        #res.tau2 <- inla.zmarginal(marg.tau2, TRUE)
        #res.tau2$mode <- inla.mmarginal(marg.tau2)
      }else{
        if(family=="gamma") res.hyper <- summary(resf)$hyperpar[1,]
      }
    }
    res.sigma <- NULL
    resphi <- NULL
    namesj <- names(resf$marginals.hyperpar)
    for(j in(1:k)){
      marg.sigmasqj <- inla.tmarginal(function(x) { (exp(x))^2},
                                     resf$marginals.hyperpar[j][[1]])
      res.sigmaj <- inla.zmarginal(marg.sigmasqj, TRUE)
      res.sigmaj$mode <- inla.mmarginal(marg.sigmasqj)

      
      marg.phij <- inla.tmarginal(function(x) { 2/exp(x)},
                                  resf$marginals.hyperpar[j+k][[1]])
      resphij <- inla.zmarginal(marg.phij, TRUE)
      resphij$mode <- inla.mmarginal(marg.phij)
      
      res.sigma <- rbind(res.sigma, unlist(res.sigmaj))
      resphi <- rbind(resphi, unlist(resphij))
    }
    

    
    if ( length(resf$marginals.hyperpar$"Precision for the Gaussian observations") > 0) {
      res.transf <- rbind(unlist(res.tau1), unlist(res.sigma1),unlist(resphi), unlist(res.sigma2),unlist(resphi2))
      res.transf <- data.frame(res.transf)
      row.names(res.transf) <- c("var.nugget1", "var.marg1", "phi1","var.marg2", "phi2")
      #row.names(res.transf) <- c("var.nugget1","var.nugget2", "var.marg", "phi","var.marg2", "phi2")
      res.transf <- res.transf[,c(-4,-6)]

      
    } else{
      res.transf <- rbind( res.sigma ,resphi)
      res.transf <- data.frame(res.transf)
      row.names(res.transf) <- paste(c( rep("var.margj",k), rep("phij",k)),1:k)
      res.transf <- res.transf[,c(-4,-6)]
      if(family=="gamma") {
        names(res.hyper) <- names(res.transf)
        res.transf <- rbind(res.hyper, res.transf)
      }
    }
  }

  }



   ## provides a summary for a inla object
   ret <- list()
   
   if (!is.null(resf$summary.fixed) && length(resf$summary.fixed) > 0) {
     ret <- c(ret, list(fixed = round(as.matrix(resf$summary.fixed), digits = 3)))
   }
   
   
   if (!is.null(resf$summary.hyperpar) && length(resf$summary.hyperpar) > 0) {
     ret <- c(ret, list(hyperpar = round(res.transf, digits = 3)))
   }

   indjk <- seq((2*k+1),(2*k+k*(k-1)/2))
   hyperlambda <- list(hyperparLMC = round(as.matrix(resf$summary.hyperpar[indjk,]), digits = 3))
     rownames(hyperlambda$hyperparLMC) <- paste(rep("lambdajk",length(indjk)), indjk)
     ret <- c(ret,  hyperlambda)
   
   return(ret)
} 




