% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doit.R
\docType{package}
\name{doit}
\alias{doit}
\alias{doit-package}
\title{doit: A package for approximate Bayesian computation using
design-of-experiments based interpolation technique}
\description{
The package uses as input a data frame of parameter values theta and
corresponding evaluations of an unnormalised probability density function
f(theta), as is a common situation in Bayesian inference problems. By
interpolating the target density using Gaussian kernels, the normalisation
constant, as well as marginal densities and expectations can be
approximated.
}
\section{Fitting methods}{


doit_estimate_w obtains the optimal width of the Gaussian interpolation
kernels by cross validation.

doit_fit estimates the parameters for the DoIt approximation

doit_update updates a fitted DoIt approximation by a new parameter value and
function evaluation.
}

\section{Approximation methods}{


doit_approx evaluates the DoIt approximation at different parameter values

doit_marginal approximates the marginal distribution

doit_marginal_A approximates the marginal distribution of a linear
transformation of the inputs

doit_expectation and doit_variance approximate expectation and variance

doit_integral approximates the integral under the unnormalised density
}

