# GLMC-blockNNGP 

## Z. Quiroz, M. Prates, C. Gonzáles and H. Rue.

We provide the R code needed to run simulation 2 of "Fast Bayesian inference of generalized linear blockNNGP model of coregionalization through INLA". 

You have to add "INLA"  from https://www.r-inla.org/download-install.  

Main code (with example usage) is denoted by [main]. 

- [main]run_all.R: fit the  models.
- newQGLMCrgeneric.R: auxiliary code that contains blockNNGP functions to run blockNNGP-GLMC model.  
- blockNNGPrgeneric.R: INLA-rgeneric code for  blockNNGP-GLMC spatial random effect.
- NNGPrgeneric.R: INLA-rgeneric code for  NNGP-GLMC spatial random effect.
- summary_functions.R: to print summary (similar to the INLA summary function) of reparameterized parameters for blockNNGP-GLMC model.
- utils.R: auxiliary functions.
