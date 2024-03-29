% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doit.R
\name{doit_fit}
\alias{doit_fit}
\title{Fit a DoIt object}
\usage{
doit_fit(design, w = NULL)
}
\arguments{
\item{design}{A data frame of design points and corresponding function
evaluations. Must contain a column named `f` with function values. The other
columns are treated as design points.}

\item{w}{vector of variances used for the Gaussian kernels. If `NULL`
(the default), `w` is calculated by `doit_estimate_w`.}
}
\value{
An object of class `doit` used in further calculations.
}
\description{
Fit the parameters of the DoIt approximation.
}
\examples{
design = data.frame(x=rnorm(10), y=rnorm(10))
design$f = with(design, exp(-0.5*(x+y)^2))
fit = doit_fit(design)

}
