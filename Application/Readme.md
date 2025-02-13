# GLMC-blockNNGP 

## Zaida Quiroz, Marcos Prates and Carlos Pizango.

We provide the R code needed to run  the application of "Fast Bayesian inference of generalized linear blockNNGP model of coregionalization through INLA". 

You have to add "INLA"  from https://www.r-inla.org/download-install.  

Download BSS data: data.est.Rdata and data.pred.Rdata

Main code (with example usage) is denoted by [main]. 

- [main]runblockNNGP-GLMC.R: fit the  models to the BSS data.
- blockNNGPfunctions.R: auxiliary code that contains blockNNGP functions to run blockNNGP-GLMC and NNGP-GLMC models. 
- blockNNGPfunctions2.R: auxiliary code that contains blockNNGP functions to run Palmí-Perales model using blockNNGP and NNGP. 
- blockNNGPrgeneric.R: INLA-rgeneric code for  blockNNGP spatial random effect.
- NNGPrgeneric.R: INLA-rgeneric code for  NNGP spatial random effect.
- summary_functions.R: to print summary (similar to the INLA summary function) of reparameterized parameters for blockNNGP-GLMC and NNGP-GLMC models.
- summary_functions2.R: to print summary (similar to the INLA summary function) of reparameterized parameters for Palmí-Perales(2023)  model  using blockNNGP and NNGP.
- utils.R: auxiliary functions.
