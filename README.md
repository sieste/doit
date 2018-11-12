# Implementation of the DoIt method 

## Papers

**Bayesian Computation Using Design of Experiments-Based Interpolation
Technique**, Joseph, 2012, Technometrics, 10.1080/00401706.2012.680399

improvement:

**A Note on Nonnegative {DoIt} Approximation**, Joseph, 2012, Technometrics,
10.1080/00401706.2012.759154

## The method

- essentially an improvement of the Laplace approximation 
- approximate target density by mixture of gaussian kernels at appropriately chosen design points
- can calculate normalisation constant, expectation, variance, marginals
- use design-of-experiments technique to select new design points



## Files

- [doit1-1d.md](doit1-1d.md): 1d example from paper 1
- [doit2-1d.md](doit2-1d.md): 1d example from paper 2, solving the negativity problem
- [doit2-2d.md](doit2-2d.md): 2d example from paper 2
- [doit2-sequential.md](doit2-sequential.md): sequential updating
- [doit_update.md](doit_update.md): new implementation, now with posterior summary measures


