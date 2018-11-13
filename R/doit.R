#' Estimate kernel widths for DoIt
#'
#' Estimate the optimal kernel bandwith by minimising the weighted
#' mean squared cross validation error.
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
#' @export
#'
doit_estimate_w = function(design, w_0=NULL, optim_control=NULL) {
  stopifnot(nrow(design) > 1)
  m     = nrow(design)
  theta = as.matrix(design[, names(design) != 'f', drop=FALSE])
  dd    = ncol(theta)
  ff    = design$f

  GGfun = function(w) {
    mat_ = matrix(0, nrow=m, ncol=m)
    for (ii in seq_len(dd)) {    
      mat_ = mat_ - outer(theta[,ii], theta[,ii], '-')^2 / (2 * w[ii])
    }
    return(drop(exp(mat_)))
  }

  # leave-one-out MSE as function of the kernel width
  wmscv = function(w) {
    GGinv_ = solve(GGfun(w))
    ee_    = drop(1/diag(GGinv_) * GGinv_ %*% sqrt(ff))
    wmscv  = sum(ee_ * diag(GGinv_) * ee_) / m
    return(wmscv)
 }

  # minimise wmscv wrt w, use `optimize` for 1d and `optim` for >1d
  if (dd == 1) {
    w = optimize(wmscv, c(1e-6, diff(range(theta))))$minimum
  } else {
    if (is.null(w_0)) { # use component-wise silvermans rule as starting point 
      w_0 = 1.06 * m^(-.2) * apply(theta, 2, sd)
    }
    w = optim(w_0, wmscv, control=optim_control)$par
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
  GG    = GGfun(theta, theta, w)
  GGinv = solve(GG)
  GG2   = sqrt(GG)
  bb    = drop(solve(GG, sqrt(ff)))
  ee    = drop(1/diag(GGinv) * GGinv %*% sqrt(ff))
  d_ij  = tcrossprod(bb) * GG2
  sum_d_ij = sum(d_ij)

  ans = list(GGfun=GGfun, theta=theta, w=w, 
             ff=ff, bb=bb, GG=GG, GG2=GG2, ee=ee, GGinv=GGinv,
             d_ij=d_ij, sum_d_ij = sum_d_ij)

  class(ans) = c('doit', class(ans))
  return(ans)
}


#' DoIt approximation of the posterior density
#'
#' Evaluate mean and variance of the DoIt approximation of the target density
#' at a number of evaluation points.
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param theta_eval A data frame of evaluation points.
#' @return A data frame of the provided evaluation points and the corresponding
#' mean and variance of the DoIt approximation. 
#' @export
#'
doit_approx = function(doit, theta_eval) with(doit, {
  bGG2b_ = drop(bb %*% GG2 %*% bb)
  gg_    = GGfun(theta_eval, theta, w)
  ggbb2_ = drop(gg_ %*% bb)^2
  ee_    = ggbb2_ / (sqrt(prod(pi*w)) * bGG2b_)
  vv_    = ggbb2_ * drop(1 - rowSums(gg_ * (gg_ %*% GGinv)))
  return(as.data.frame(cbind(theta_eval, ee=ee_, vv=pmax(0, vv_))))
})


#' Propose a new design point for the DoIt approximation
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return A parameter value.
#' @export
#'
doit_propose_new = function(doit) with(doit, {
  fn = function(r) {
    gg = GGfun(r, theta, w)
    -1 * sum(bb * gg)^2 * drop(1 - gg %*% GGinv %*% gg)
  }
  vv =  (sqrt(ff) - ee)^2 / diag(GGinv)
  r0 = drop(theta[which.max(vv), ])
  if (ncol(theta) == 1) {
    theta_new = optimize(fn, range(theta)+c(-1,1)*diff(range(theta)))$minimum
  } else {
    theta_new = optim(r0, fn)$par
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
  d_ij_up = tcrossprod(bb_up) * GG2_up
  # write new object
  doit$theta    = theta_up
  doit$ff       = ff_up
  doit$bb       = bb_up
  doit$GG       = GG_up
  doit$GG2      = GG2_up
  doit$ee       = ee_up
  doit$GGinv    = GGinv_up
  doit$d_ij     = d_ij_up
  doit$sum_d_ij = sum(d_ij_up)
  return(doit)
})


#' DoIt approximation of the marginal posterior distribution
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param k parameter dimension
#' @param theta_eval Evaluation points at which to approximate the marginal
#' distribution. If `NULL` (the default) the original design points are used.
#' @return A data frame of the provided evaluation points and the corresponding
#' DoIt approximation of the marginal posterior.
#' @export
#'
doit_marginal = function(doit, k, theta_eval=NULL) with(doit, {
  if (is.numeric(k)) stopifnot(k > 0, k <= ncol(theta))
  if (is.character(k)) stopifnot(k %in% colnames(theta))
  if (is.null(theta_eval)) theta_eval = theta[, k]
  nu_ij = 0.5 * outer(theta[,k], theta[,k], '+')
  sd_ = sqrt(w[k]/2)
  ans = sapply(theta_eval, function(tt) {
    phi_ij = dnorm(tt, nu_ij, sd_)
    return(sum(d_ij * phi_ij) / sum_d_ij)
  })
  return(data_frame(par=colnames(theta)[k], theta=theta_eval, posterior=ans))
})


#' DoIt approximation of all marginal posterior distributions.
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @param theta_eval Evaluation points at which to approximate the marginal
#' distribution. If `NULL` (the default) the original design points are used.
#' @return A data frame of the provided evaluation points and the corresponding
#' DoIt approximation of the marginal posterior.
#' @export
#'
doit_marginals = function(doit, theta_eval=NULL) {
  margs = lapply(1:ncol(doit$theta), function(kk) 
                 doit_marginal(doit, kk, theta_eval))
  margs = do.call(rbind, margs)
  return(margs)
}



#' DoIt approximation of posterior expectation
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return Vector of posterior expectations
#' @export
#'
doit_expectation = function(doit) with(doit, {
  nu_ij = lapply(as.data.frame(theta), function(tt) 0.5 * outer(tt, tt, '+'))
  e_theta = sapply(nu_ij, function(nu) sum(d_ij * nu) / sum_d_ij)
  return(e_theta)
})


#' DoIt approximation of posterior variance
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return Posterior (co-)variance matrix.
#' @export
#'
doit_variance = function(doit) with(doit, {
  nu_ij = lapply(as.data.frame(theta), function(tt) 0.5 * outer(tt, tt, '+'))
  dd = ncol(theta)
  kk = expand.grid(ii = 1:dd, jj = 1:dd) 
  kk = kk[with(kk, jj>=ii), ]
  vv = sapply(1:nrow(kk), function(ll) 
         sum(d_ij * nu_ij[[ kk[ll, 'ii'] ]] * nu_ij[[ kk[ll, 'jj'] ]]) / 
         sum_d_ij)
  v_theta = matrix(0, dd, dd)
  for (l in 1:nrow(kk)) {
    v_theta[kk[l,'ii'], kk[l,'jj']] = v_theta[kk[l,'jj'], kk[l,'ii']] = vv[l]
  }
  rownames(v_theta) = colnames(v_theta) = colnames(theta)
  return(v_theta)
})


