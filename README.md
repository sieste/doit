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



## Demo Files

- [doit1-1d.md](old_stuff/doit1-1d.md): 1d example from paper 1
- [doit2-1d.md](old_stuff/doit2-1d.md): 1d example from paper 2, solving the negativity problem
- [doit2-2d.md](old_stuff/doit2-2d.md): 2d example from paper 2
- [doit2-sequential.md](old_stuff/doit2-sequential.md): sequential updating
- [doit_update.md](old_stuff/doit_update.md): new implementation, now with posterior summary measures



## Todo

- k in doit_marginal should accept parameter name 
- test that doit_marginal behaves as expected with different settings for
  theta_eval
- make nu_ij part of doit object
- add @seealso tags
- add links to roxygen blocks (had/r-pkgs/man)
- add example data
- replace `crossing` by `expand`
- in vignette plot truth and approximation side-by-side instead of contour over raster
- add argument `near` to specify starting point in `doit_propose_new`
- write vignettes for 1d and 2d function, using the package functions
- make `GGfun` more efficient by accounting for symmetry
- in `doit_fit` and `doit_update` check for duplicate design points



