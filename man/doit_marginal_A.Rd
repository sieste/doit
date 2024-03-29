% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doit.R
\name{doit_marginal_A}
\alias{doit_marginal_A}
\title{DoIt approximation of the marginal of a linear transformation}
\usage{
doit_marginal_A(doit, A = NULL, theta_eval = NULL)
}
\arguments{
\item{doit}{An object of class `doit`, see function `doit_fit`.}

\item{A}{The transformation matrix.}

\item{theta_eval}{Evaluation points at which to apply the linear
transformation and then approximate the marginal distribution. If `NULL`
(the default) the original design points are used.}
}
\value{
A data frame of the transformed evaluation points and the corresponding
DoIt approximation of the marginal density.
}
\description{
Approximate the marginal density of a linear transformation `A %*% theta`.
}
\examples{
design = data.frame(x=rnorm(10), y=rnorm(10))
design$f = with(design, exp(-0.5*(x+y)^2))
fit = doit_fit(design)
# distribution of x + y
mar_xpy = doit_marginal_A(fit, A=matrix(c(1,1), 1, 2))

}
