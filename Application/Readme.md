# GLMC-blockNNGP 

## Z. Quiroz, M. Prates, C. Pizango and H. Rue.

We provide the R code needed to run  the application of "Fast Bayesian inference of generalized linear blockNNGP model of coregionalization through INLA". 

You have to add "INLA"  from https://www.r-inla.org/download-install.  

Download BSS data: data.est.Rdata and data.pred.Rdata

Main code (with example usage) is denoted by [main]. 

- [main]runblockNNGP-GLMC.R: fit the  models to the BSS data.
- blockNNGPfunctions.R: auxiliary code that contains blockNNGP functions to run blockNNGP-GLMC model.  
- blockNNGPrgeneric.R: INLA-rgeneric code for  blockNNGP spatial random effect.
- summary_functions.R: to print summary (similar to the INLA summary function) of reparameterized parameters for blockNNGP-GLMC model.
- utils.R: auxiliary functions.

```
Zaida Quiroz
Science Department
School of Mathematics
Pontificia Universidad Católica del Perú
Lima, Perú
E-mail: zquiroz@pucp.edu.pe
```
