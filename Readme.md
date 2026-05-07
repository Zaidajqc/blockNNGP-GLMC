
# Fast  blockNNGP-GLMC Modeling through INLA

R code needed to run simulations and  application study of generalized linear model of coregionalization (blockNNGP-GLMC)  model  using Integrated Nested Laplace Approximation (INLA) as proposed in "Fast Bayesian inference of generalized linear blockNNGP model of coregionalization through INLA" by   Z. Quiroz, M.Prates, C. Gonzáles and H. Rue (https://....). 

We extend the reformulation of the LMC presented by Krainski et al. (2018) to define the generalized LMC model (GLMC), in which the precision matrix of all N spatial random effects is efficiently computed as a sum of Kronecker products, as proposed by Alie et al. (2024). However, we remark that the GLMC is defined as a linear combination of spatial processes, where each univariate kth Gaussian spatial process has its own covariance function, that is, its own range and spatial variance. Due to this feature, we can present a computationally efficient extension of GLMC for large data, merging it with the blockNNGP. The new model, called blockNNGP-GLMC, fits multivariate Gaussian and non-Gaussian spatially dependent response variables, and the spatial random effects for each kth response variable are defined via a blockNNGP prior. We implement fast Bayesian inference for blockNNGP-GLMC models using INLA, which  particularly explores the sparsity of the precision matrix of the Gaussian random effects.

In the application, specifically, we estimate and predict the abundance of the American Robin, Chipping Sparrow, and European Starling, one of the most abundant species in North America.

