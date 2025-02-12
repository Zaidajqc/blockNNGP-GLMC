#' @title summary.blockNNGP()
#' @description This function returns summary of posterior parameters for blockNNGP and NNGP models.

summary.blockNNGP <- function(name.prior, resf, data1,n.blocks, num.nb,prior.set=NULL){

if(name.prior=="pc.prior")  {
  if (!is.null(resf$summary.hyperpar) && length(resf$summary.hyperpar) > 0) {
  
    if ( length(resf$marginals.hyperpar$"Precision for the Gaussian observations") > 0) {
      marg.tau1     <- inla.tmarginal(function(x) { 1/x},
                                  resf$marginals.hyperpar[[1]])
      res.tau2      <- inla.zmarginal(marg.tau1, TRUE)
      res.tau2$mode <- inla.mmarginal(marg.tau1)
    }
    marg.sigmasq1   <- inla.tmarginal(function(x) { (exp(x))^2},
                                    resf$marginals.hyperpar$"Theta1 for spatial.field.y1")
    res.sqsigma1      <- inla.zmarginal(marg.sigmasq1, TRUE)
    res.sqsigma1$mode <- inla.mmarginal(marg.sigmasq1)
    marg.phi1         <- inla.tmarginal(function(x) { 2/exp(x)},
                                resf$marginals.hyperpar$"Theta2 for spatial.field.y1")
    resphi1       <- inla.zmarginal(marg.phi1, TRUE)
    resphi1$mode  <- inla.mmarginal(marg.phi1)
    marg.sigmasq2 <- inla.tmarginal(function(x) { (exp(x))^2},
                                    resf$marginals.hyperpar$"Theta1 for spatial.field.y2")
    res.sqsigma2      <- inla.zmarginal(marg.sigmasq2, TRUE)
    res.sqsigma2$mode <- inla.mmarginal(marg.sigmasq2)
    marg.phi2         <- inla.tmarginal(function(x) { 2/exp(x)},
                                resf$marginals.hyperpar$"Theta2 for spatial.field.y2")
    resphi2       <- inla.zmarginal(marg.phi2, TRUE)
    resphi2$mode  <- inla.mmarginal(marg.phi2)
    
    if ( length(resf$marginals.hyperpar$"Precision for the Gaussian observations") > 0) {
      res.transf  <- rbind(unlist(res.tau2), unlist(res.sqsigma1),unlist(resphi1), unlist(res.sqsigma2),unlist(resphi2))
      res.transf  <- data.frame(res.transf)
      row.names(res.transf) <- c("var.nugget", "var.marg1", "phi1", "var.marg2", "phi2")
      res.transf  <- res.transf[,c(-4,-6)]
    } else{
      res.transf <- as.matrix(rbind( unlist(res.sqsigma1),unlist(resphi1),unlist(res.sqsigma2),unlist(resphi2)))
      res.transf <- data.frame(res.transf)
      row.names(res.transf) <- c( "var.marg1", "phi1","var.marg2", "phi2")
      res.transf <- res.transf[,c(-4,-6)]
    }
  }  
}

   ## provides a summary for a inla object
   ret <- list()
   
   if (!is.null(resf$summary.fixed) && length(resf$summary.fixed) > 0) {
     ret  <- c(ret, list(fixed = round(as.matrix(resf$summary.fixed), digits = 3)))
   }
   
   
   if (!is.null(resf$summary.hyperpar) && length(resf$summary.hyperpar) > 0) {
     ret  <- c(ret, list(hyperpar = round(res.transf, digits = 3)))
   }
   
   return(ret)
} 

plot_birds = function(resf, data1,variable.est,n.blocks, num.nb ){
  library(rworldmap )
  
  worldmap    <- getMap()
  study_area  <- subset(worldmap, NAME %in% c('United States', "Canada"))
  
  loc <- cbind(data1$loc1, data1$loc2)
  int.elev2 <- mba.surf(cbind(loc,variable.est), 100, 100, extend=TRUE)$xyz.est
  
  image.plot(int.elev2  , main= paste('M = ',n.blocks,', nb=',num.nb),
             xaxs = 'r', yaxs = 'r',
             xlim = range(loc[,1]),
             ylim = range(loc[,2]),
             xlab='Longitude', ylab='Latitude')
  lines(study_area)
}


