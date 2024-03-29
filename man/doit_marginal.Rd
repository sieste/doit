% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doit.R
\name{doit_marginal}
\alias{doit_marginal}
\title{DoIt approximation of the marginal density}
\usage{
doit_marginal(doit, param, theta_eval = NULL)
}
\arguments{
\item{doit}{An object of class `doit`, see function `doit_fit`.}

\item{param}{Parameter name whose marginal density is calculated. The
original design data frame must have a column with that name.  If integer,
`param` is interpreted as the column index in the *original* design.}

\item{theta_eval}{Data frame or matrix of parameter values at which to
approximate the marginal distribution. Should have a column named `param`.
If `NULL` (the default) the original design is used.}
}
\value{
A data frame of the provided evaluation points and the corresponding
DoIt approximation of the marginal density.
}
\description{
Approximate the marginal density of one element of `theta`.
}
\examples{
design   = data.frame(x=rnorm(10), y=rnorm(10))
design$f = with(design, exp(-0.5*(x+y)^2))
fit      = doit_fit(design)
mar_x    = doit_marginal(fit, 'x')

}
