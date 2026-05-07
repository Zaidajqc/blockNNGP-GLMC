#summary for blockNNGP 

summary.blockNNGP_LMC <- function(name.prior, resf, data1,n.blocks, 
                                  num.nb,prior.set=NULL,family="gaussian"){
  
if(name.prior=="pc.prior")  {
  if (!is.null(resf$summary.hyperpar) && length(resf$summary.hyperpar) > 0) {
    if ( length(names(resf$marginals.hyperpar)) > 0) {
      if ( length(resf$marginals.hyperpar$"Precision for the Gaussian observations") > 0) {
        marg.tau1 <- inla.tmarginal(function(x) { 1/x},
                                    resf$marginals.hyperpar$"Precision for the Gaussian observations")
        res.tau1 <- inla.zmarginal(marg.tau1, TRUE)
        res.tau1$mode <- inla.mmarginal(marg.tau1)
      }else{
        if(family=="gamma") res.hyper <- summary(resf)$hyperpar[1,]
      }
    }
    marg.sigmasq1 <- inla.tmarginal(function(x) { (exp(x))^2},
                                    resf$marginals.hyperpar$"Theta1 for spatial.field.y")
    res.sigma1 <- inla.zmarginal(marg.sigmasq1, TRUE)
    res.sigma1$mode <- inla.mmarginal(marg.sigmasq1)
    marg.phi1 <- inla.tmarginal(function(x) { 2/exp(x)},
                                resf$marginals.hyperpar$"Theta2 for spatial.field.y")
    resphi <- inla.zmarginal(marg.phi1, TRUE)
    resphi$mode <- inla.mmarginal(marg.phi1)
    
    marg.sigmasq2 <- inla.tmarginal(function(x) { (exp(x))^2},
                                    resf$marginals.hyperpar$"Theta3 for spatial.field.y")
    res.sigma2 <- inla.zmarginal(marg.sigmasq2, TRUE)
    res.sigma2$mode <- inla.mmarginal(marg.sigmasq2)
    marg.phi2 <- inla.tmarginal(function(x) { 2/exp(x)},
                                resf$marginals.hyperpar$"Theta4 for spatial.field.y")
    resphi2 <- inla.zmarginal(marg.phi2, TRUE)
    resphi2$mode <- inla.mmarginal(marg.phi2)
    
    
    
    if ( length(resf$marginals.hyperpar$"Precision for the Gaussian observations") > 0) {
      res.transf <- rbind(unlist(res.tau1), unlist(res.sigma1),unlist(resphi), unlist(res.sigma2),unlist(resphi2))
      res.transf <- data.frame(res.transf)
      row.names(res.transf) <- c("var.nugget1", "var.marg1", "phi1","var.marg2", "phi2")
      res.transf <- res.transf[,c(-4,-6)]

    } else{
      res.transf <- rbind( unlist(res.sigma1),unlist(resphi), unlist(res.sigma2),unlist(resphi2))
      res.transf <- data.frame(res.transf)
      row.names(res.transf) <- c( "var.marg1", "phi1","var.marg2", "phi2")
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

   # now just for poisson, binomial.
   if (!is.null(resf$marginals.hyperpar$"Theta5 for spatial.field.y") && length(resf$summary.hyperpar) > 0) {
     hyperlambda <- list(hyperparLMC = round(as.matrix(resf$summary.hyperpar[5,]), digits = 3))
     rownames(hyperlambda$hyperparLMC) <- "lambda"
     ret <- c(ret,  hyperlambda)

     }   
   
   return(ret)
} 


