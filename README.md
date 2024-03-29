# doit: Bayesian Computation using Design of Experiments-based Interpolation Technique in R

## Overview

This R package implements the design of experiments-based interpolation
technique (DoIt, Joseph 2012) for approximate Bayesian computations. 

The method uses evaluations of an unnormalised density at a space-filling
design of parameter values. Normalisation is achieved by approximating the
target density by a weighted sum of Gaussian kernels centered on the design
points. 

DoIt approximates the joint density, marginal densities, as well as expecations
and variances. The package contains functions to optimally choose additional
design points, and to calculate the optimal kernel width by cross validation.

![Example plot of 2d DoIt approximation](fig/doit-2d.png)

_Figure: DoIt approximation of a complicated 2-dimensional density. See `vignette('doit_2d')` for details._


## Installation

```r
remotes::install_github('sieste/doit', build_opts=NULL)
```

To install the package without using `remotes`, run the following shell
commands:

```sh
git clone git@github.com:sieste/doit.git
cd doit
R CMD build .
R CMD INSTALL doit_*.tar.gz
```


## Vignettes

The usage of the package is documented in 2 vignettes, where results from the
original papers are reproduced.

```r
vignette('doit_1d') # 1d example from Joseph (2012)
vignette('doit_2d') # 2d example from Joseph (2012)
```


## References

Joseph (2012) _Bayesian Computation Using Design of Experiments-Based
Interpolation Technique_, Technometrics,
[10.1080/00401706.2012.680399](http://dx.doi.org/10.1080/00401706.2012.680399)

Joseph (2012) _A Note on Nonnegative DoIt Approximation_, Technometrics,
[10.1080/00401706.2012.759154](http://dx.doi.org/10.1080/00401706.2012.759154)



