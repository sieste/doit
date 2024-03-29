#' doit: A package for approximate Bayesian computation using
#' design-of-experiments based interpolation technique
#'
#' The package uses as input a data frame of parameter values theta and
#' corresponding evaluations of an unnormalised probability density function
#' f(theta), as is a common situation in Bayesian inference problems. By
#' interpolating the target density using Gaussian kernels, the normalisation
#' constant, as well as marginal densities and expectations can be
#' approximated.
#'
#' @section Fitting methods:
#'
#' doit_estimate_w obtains the optimal width of the Gaussian interpolation
#' kernels by cross validation.
#'
#' doit_fit estimates the parameters for the DoIt approximation
#'
#' doit_update updates a fitted DoIt approximation by a new parameter value and
#' function evaluation.
#'
#' @section Approximation methods:
#'
#' doit_approx evaluates the DoIt approximation at different parameter values
#'
#' doit_marginal approximates the marginal distribution
#'
#' doit_marginal_A approximates the marginal distribution of a linear
#' transformation of the inputs
#'
#' doit_expectation and doit_variance approximate expectation and variance
#'
#' doit_integral approximates the integral under the unnormalised density
#'
#' @docType package
#' @name doit
NULL



#' Estimate kernel widths for DoIt
#'
#' Estimate the optimal kernel bandwith for the approximation by minimising the
#' weighted mean squared cross validation error.
#' 
#' @param design A data frame of design points and corresponding function
#' evaluations. Must contain a column named `f` with function values. The other
#' columns are treated as design points.
#' @param w_0 Starting point for the optimisation routine. Either NULL
#' (default) in wwhich case an educated guess is made, or otherwise a vector of
#' the same dimension as the parameter space. 
#' @param optim_control A list passed to `optimize` (for 1d problems) or
#' `optim` (if the parameter has dimension 2 or more) 
#' @return A vector of optimal kernel widths.
#'
#' @examples 
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' doit_estimate_w(design)
#'
#' @importFrom stats optim optimize sd
#' @export
#'
doit_estimate_w = function(design, w_0=NULL, optim_control=NULL) {
  stopifnot(nrow(design) > 1)
  m     = nrow(design)
  theta = as.matrix(design[, names(design) != 'f', drop=FALSE])
  dd    = ncol(theta)
  ff    = design$f

  theta2 = lapply(1:dd, function(ii) outer(theta[,ii], theta[,ii], '-')^2 / 2)

  GGfun = function(w) {
    mat_ = matrix(0, nrow=m, ncol=m)
    for (ii in seq_len(dd)) {    
      mat_ = mat_ - theta2[[ii]] / w[ii]
    }
    return(drop(exp(mat_)))
  }

  # leave-one-out MSE as function of the kernel width
  wmscv = function(log_w) {
    GGinv_ = solve(GGfun(exp(log_w)))
    ee_    = drop(1/diag(GGinv_) * GGinv_ %*% sqrt(ff))
    wmscv  = sum(ee_ * diag(GGinv_) * ee_) / m
    return(wmscv)
 }

  # minimise wmscv wrt w, use `optimize` for 1d and `optim` for >1d
  if (dd == 1) {
    opt = optimize(wmscv, c(log(1e-6), log(diff(range(theta)))))
    w = exp(opt$minimum)
  } else {
    if (is.null(w_0)) { # use component-wise silvermans rule as starting point 
      w_0 = 1.06 * m^(-.2) * apply(theta, 2, sd)
    }
    opt = optim(log(w_0), wmscv, control=optim_control)
    w   = exp(opt$par)
  }
  names(w) = colnames(theta)
  return(w)
}


#' Fit a DoIt object
#'
#' Fit the parameters of the DoIt approximation.
#'
#' @param design A data frame of design points and corresponding function
#' evaluations. Must contain a column named `f` with function values. The other
#' columns are treated as design points.
#' @param w vector of variances used for the Gaussian kernels. If `NULL`
#' (the default), `w` is calculated by `doit_estimate_w`.
#' @return An object of class `doit` used in further calculations.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#'
#' @export
#'
doit_fit = function(design, w=NULL) {
  m     = nrow(design)
  theta = as.matrix(design[, names(design) != 'f', drop=FALSE])
  dd    = ncol(theta)
  ff    = design$f

  GGfun = function(xx, yy, w) {
    if (is.null(dim(xx))) xx = matrix(xx, nrow=1)
    if (is.null(dim(yy))) yy = matrix(yy, nrow=1)
    stopifnot(ncol(xx) == ncol(yy), ncol(xx) == length(w))
    mat_ = matrix(0, nrow=nrow(xx), ncol=nrow(yy))
    for (ii in seq_along(w)) {    
      mat_ = mat_ - outer(xx[,ii], yy[,ii], '-')^2 / (2 * w[ii])
    }
    return(drop(exp(mat_)))
  }
  
  if(is.null(w)) {
    w = doit_estimate_w(design)
  }

  # calculate parameters
  GG       = GGfun(theta, theta, w)
  GGinv    = solve(GG)
  GG2      = sqrt(GG)
  bb       = drop(solve(GG, sqrt(ff)))
  ee       = drop(1/diag(GGinv) * GGinv %*% sqrt(ff))

  ij = expand.grid(i=1:m, j=1:m)
  ij = ij[ with(ij, i>=j), ]
  diag_w = diag(x=1/(2*w), nrow=dd)
  theta_imj2 = (theta[ij$i, , drop=FALSE] - theta[ij$j, , drop=FALSE])^2
  d_ij = bb[ij$i] * bb[ij$j] * exp(-0.5*rowSums(theta_imj2 %*% diag_w))
  d_ij = ifelse(ij$i == ij$j, d_ij, 2*d_ij)
  nu_ij = 0.5 * (theta[ij$i, , drop=FALSE] + theta[ij$j, , drop=FALSE])
  sum_d_ij = sum(d_ij)

  ans = list(GGfun=GGfun, theta=theta, w=w, 
             ff=ff, bb=bb, GG=GG, GG2=GG2, ee=ee, GGinv=GGinv,
             d_ij=d_ij, sum_d_ij=sum_d_ij, nu_ij=nu_ij, ij=ij)

  class(ans) = c('doit', class(ans))
  return(ans)
}


#' DoIt approximation of the target density
#'
#' Evaluate mean and variance of the DoIt approximation of the target density
#' at a number of evaluation points.
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param theta_eval A data frame of evaluation points.
#' @return A data frame of the provided evaluation points and the corresponding
#' mean and variance of the DoIt approximation. 
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' theta_eval = expand.grid(x=seq(-2,2,.1), y=seq(-2,2,.1))
#' approx = doit_approx(fit, theta_eval)
#'
#' @export
#'
doit_approx = function(doit, theta_eval) with(doit, {
  theta_eval = as.matrix(theta_eval)
  bGG2b_ = drop(bb %*% GG2 %*% bb)
  gg_    = GGfun(theta_eval, theta, w)
  ggbb2_ = drop(gg_ %*% bb)^2
  ee_    = ggbb2_ / (sqrt(prod(pi*w)) * bGG2b_)
  vv_    = ggbb2_ * drop(1 - rowSums(gg_ * (gg_ %*% GGinv)))
  return(as.data.frame(cbind(theta_eval, dens_approx=ee_, variance=pmax(0, vv_))))
})


#' Propose a new design point for the DoIt approximation
#'
#' Find the point with largest estimation variance by numerical optimisation,
#' using the current design point with largest leave-one-out prediction error
#' as the starting point.
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param theta_0 Data frame or vector. Initial value used by the optimisation
#' routine. Ignored for 1d problems.
#' @return A parameter value.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' theta_new = doit_propose_new(fit)
#'
#' @export
#'
doit_propose_new = function(doit, theta_0=NULL) with(doit, {
  fn = function(r) {
    gg = GGfun(r, theta, w)
    -1 * sum(bb * gg)^2 * drop(1 - gg %*% GGinv %*% gg)
  }
  if (is.null(theta_0)) {
    vv =  (sqrt(ff) - ee)^2 / diag(GGinv)
    theta_0 = drop(theta[which.max(vv), ])
  } else {
    theta_0 = unlist(theta_0)
  }
  stopifnot(length(theta_0) == ncol(theta))
  if (ncol(theta) == 1) {
    theta_new = optimize(fn, range(theta)+c(-1,1)*diff(range(theta)))$minimum
  } else {
    theta_new = optim(theta_0, fn)$par
  }
  names(theta_new) = colnames(theta)
  return(data.frame(as.list(theta_new)))
})


#' Update a DoIt approximation
#'
#' Update a `doit` object by adding a new design point, but without changing
#' the kernel width.
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param design_new The new design point.
#' @return The updated `doit` object.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' theta_new = doit_propose_new(fit)
#' theta_new$f = with(theta_new, exp(-0.5*(x+y)^2))
#' fit = doit_update(fit, theta_new)
#'
#' @export
#'
doit_update = function(doit, design_new) with(doit, {
  design_new = design_new[1, ] 
  theta_new  = as.matrix(design_new[, names(design_new) != 'f', drop=FALSE])
  theta_up   = rbind(theta, theta_new)
  gg_        = GGfun(theta, theta_new, w)
  gg2_       = sqrt(gg_)
  GG_up      = rbind(cbind(GG, gg_), cbind(t(gg_), 1))
  GG2_up     = rbind(cbind(GG2, gg2_), cbind(t(gg2_), 1))
  dd_        = 1 / drop(1 - crossprod(gg_, GGinv %*% gg_))
  hh_        = GGinv %*% gg_
  GGinv_up   = rbind(cbind(GGinv + dd_*tcrossprod(hh_), - dd_*hh_),
                     cbind(-dd_ * t(hh_), dd_))
  ff_up   = c(ff, design_new$f)
  bb_up   = drop(GGinv_up %*% sqrt(ff_up))
  ee_up   = drop(1/diag(GGinv_up) * GGinv_up %*% sqrt(ff_up))

  m = nrow(theta_up)
  dd = ncol(theta_up)
  ij_up = expand.grid(i=1:m, j=1:m)
  ij_up = ij_up[ with(ij_up, i>=j), ]
  diag_w = diag(x=1/(2*w), nrow=dd)
  theta_up_imj2 = (theta_up[ij_up$i, , drop=FALSE] - theta_up[ij_up$j, , drop=FALSE])^2
  d_ij_up = bb_up[ij_up$i] * bb_up[ij_up$j] * exp(-0.5*rowSums(theta_up_imj2 %*% diag_w))
  d_ij_up = ifelse(ij_up$i == ij_up$j, d_ij_up, 2*d_ij_up)
  nu_ij_up = 0.5 * (theta_up[ij_up$i, , drop=FALSE] + theta_up[ij_up$j, , drop=FALSE])

  # write new object
  doit$theta    = theta_up
  doit$ff       = ff_up
  doit$bb       = bb_up
  doit$GG       = GG_up
  doit$GG2      = GG2_up
  doit$ee       = ee_up
  doit$GGinv    = GGinv_up
  doit$ij       = ij_up
  doit$d_ij     = d_ij_up
  doit$sum_d_ij = sum(d_ij_up)
  doit$nu_ij    = nu_ij_up

  return(doit)
})


#' DoIt approximation of the expectation
#'
#' Approximate the expectation of the vector `theta`.
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return Vector of expected values of the target density.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' expec = doit_expectation(fit)
#'
#' @export
#'
doit_expectation = function(doit) with(doit, {
  e_theta = colSums(d_ij * nu_ij) / sum_d_ij
  return(e_theta)
})


#' DoIt approximation of the marginal density
#'
#' Approximate the marginal density of one element of `theta`.
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param param Parameter name whose marginal density is calculated. The
#' original design data frame must have a column with that name.  If integer,
#' `param` is interpreted as the column index in the *original* design.
#' @param theta_eval Data frame or matrix of parameter values at which to
#' approximate the marginal distribution. Should have a column named `param`.
#' If `NULL` (the default) the original design is used.
#' @return A data frame of the provided evaluation points and the corresponding
#' DoIt approximation of the marginal density.
#'
#' @examples
#' design   = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit      = doit_fit(design)
#' mar_x    = doit_marginal(fit, 'x')
#'
#' @export
#'
doit_marginal = function(doit, param, theta_eval=NULL) with(doit, {
  stopifnot('doit' %in% class(doit), 
            length(param) == 1)
  if (is.numeric(param)) {
    stopifnot(param > 0, param <= ncol(theta))
    param = colnames(theta)[param]
  }
  if (is.character(param)) {
    stopifnot(param %in% colnames(theta))
  }
  if (is.null(theta_eval)) theta_eval = theta
  theta_eval = unique(theta_eval[, param])
  sd_ = sqrt(w[param]/2)
  dens_ = sapply(theta_eval, function(tt) {
    phi_ij = dnorm(tt, nu_ij[,param], sd_)
    return(sum(d_ij * phi_ij) / sum_d_ij)
  })
  return(data.frame(par=param, theta=theta_eval, dens_approx=dens_))
})


#' DoIt approximation of all marginal densities
#'
#' Approximate the marginal densities of all elements of the vector `theta`.
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param theta_eval Evaluation points at which to approximate the marginal
#' distribution. If `NULL` (the default) the original design points are used.
#' @return A data frame of the provided evaluation points and the corresponding
#' DoIt approximation of the marginal density.
#'
#' @examples
#' design   = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit      = doit_fit(design)
#' mar      = doit_marginals(fit)
#'
#' @export
#'
doit_marginals = function(doit, theta_eval=NULL) {
  if (!is.null(theta_eval)) theta_eval = as.matrix(theta_eval)
  margs = lapply(colnames(doit$theta), function(param) 
                 doit_marginal(doit, param, theta_eval))
  margs = do.call(rbind, margs)
  return(margs)
}


#' DoIt approximation of the variance
#'
#' Approximate the covariance matrix of the vector `theta`.
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return (Co-)variance matrix of the target density.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' vari = doit_variance(fit)
#'
#' @export
#'
doit_variance = function(doit) with(doit, {
  dd = ncol(theta)
  expec = colSums(d_ij * nu_ij) / sum_d_ij
  var_ij = lapply(1:nrow(nu_ij), function(kk) d_ij[kk] * tcrossprod(nu_ij[kk, ]))
  v_theta = Reduce(`+`, var_ij) / sum_d_ij + diag(x=w/2, nrow=dd) - tcrossprod(expec)
  rownames(v_theta) = colnames(v_theta) = colnames(theta)
  return(v_theta)
})


#' DoIt approximation of the normalisation constant
#'
#' Calculate the integral of the approximated function.
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return The integral of the approximated function over all `theta`.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' Z = doit_integral(fit)
#'
#' @export
#'
doit_integral = function(doit) with(doit, {
  ans = sqrt(prod(pi*w)) * drop(bb %*% GG2 %*% bb)
  return(ans)
})


#' DoIt approximation of the marginal of a linear transformation
#'
#' Approximate the marginal density of a linear transformation `A %*% theta`.
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param A The transformation matrix.
#' @param theta_eval Evaluation points at which to apply the linear
#' transformation and then approximate the marginal distribution. If `NULL`
#' (the default) the original design points are used.
#' @return A data frame of the transformed evaluation points and the corresponding
#' DoIt approximation of the marginal density.
#'
#' @examples
#' design = data.frame(x=rnorm(10), y=rnorm(10))
#' design$f = with(design, exp(-0.5*(x+y)^2))
#' fit = doit_fit(design)
#' # distribution of x + y
#' mar_xpy = doit_marginal_A(fit, A=matrix(c(1,1), 1, 2))
#'
#' @export
#' 
doit_marginal_A = function(doit, A=NULL, theta_eval=NULL) with(doit, {
  if (is.null(theta_eval)) theta_eval = theta
  dd = ncol(theta)
  if (is.null(A)) A = diag(ncol=dd)
  if (is.null(dim(A))) A = matrix(A, nrow=1)
  stopifnot(ncol(A) == dd)
  theta_eval = matrix(theta_eval, ncol=dd)
  tau_eval   = theta_eval %*% t(A)
  colnames(tau_eval) = paste('tau', 1:nrow(A), sep='')
  mu_ij      = nu_ij %*% t(A)
  Sigma      = A %*% diag(x=w/2, nrow=dd) %*% t(A)
  phi_ij = t(apply(mu_ij, 1, function(mu) {
    mvtnorm::dmvnorm(tau_eval, mu, Sigma)
  }))
  dens_approx = colSums(d_ij * phi_ij) / sum_d_ij
  return(data.frame(tau_eval, dens_approx=dens_approx))
})


