# BLSM

R package allowing the computation of a Bayesian latent space model for complex networks, either weighted or unweighted.

Given an observed input graph, the estimates for the latent coordinates of the nodes are obtained through a Bayesian MCMC algorithm. 

The overall likelihood of the graph depends on a fundamental probability equation, which is defined so that ties are more likely to exist between nodes whose latent space coordinates are close. 

The package is mainly based on the model by Hoff, Raftery and Handcock (JASA, 2002) and contains some extra features (e.g., removal of the Procrustean step, weights implemented as coefficients of the latent distances, 3D plots). The original code related to the above model was retrieved from https://www.stat.washington.edu/people/pdhoff/Code/hoff_raftery_handcock_2002_jasa/.

Users can inspect the MCMC simulation, create and customize insightful graphical representations or apply clustering techniques. 