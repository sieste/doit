#' Estimate kernel widths for DoIt
#'
#' Estimate the optimal kernel bandwith by minimising the weighted
#' mean squared cross validation error.
#' 
#' @param design A data frame of design points and corresponding function
#' evaluations. Must contain a column named `f` with function values. The other
#' columns are treated as design points.
#' @param sigma_0 Starting point for the optimisation routine. Either NULL
#' (default) in wwhich case an educated guess is made, or otherwise a vector of
#' the same dimension as the parameter space. 
#' @param optim_control A list passed to `optimize` (for 1d problems) or
#' `optim` (if the parameter has dimension 2 or more) 
#' @return A vector of optimal kernel widths.
#' @export
#'
doit_estimate_sigma = function(design, sigma_0=NULL, optim_control=NULL) {
  stopifnot(nrow(design) > 1)
  m     = nrow(design)
  theta = as.matrix(dplyr::select(design, -f))
  dd    = ncol(theta)
  ff    = design$f

  GGfun = function(sigma2) {
    mat_ = matrix(0, nrow=m, ncol=m)
    for (ii in seq_len(dd)) {    
      mat_ = mat_ - outer(theta[,ii], theta[,ii], '-')^2 / (2 * sigma2[ii])
    }
    return(drop(exp(mat_)))
  }

  # leave-one-out MSE as function of the kernel width
  wmscv = function(sigma2) {
    GGinv_ = solve(GGfun(sigma2))
    ee_    = drop(1/diag(GGinv_) * GGinv_ %*% sqrt(ff))
    wmscv  = sum(ee_ * diag(GGinv_) * ee_) / m
    return(wmscv)
 }

  # minimise wmscv wrt sigma2, use `optimize` for 1d and `optim` for >1d
  if (dd == 1) {
    sigma2 = optimize(wmscv, c(1e-6, diff(range(design$xx))))$minimum
  } else {
    if (is.null(sigma_0)) { # use component-wise silvermans rule as starting point 
      sigma_0 = 1.06 * m^(-.2) * apply(theta, 2, sd)
    }
    sigma2 = optim(sigma_0, wmscv)$par
  }
  return(sigma2)
}



#' Fit a DoIt object
#'
#' Fit the parameters of the DoIt approximation.
#'
#' @param design A data frame of design points and corresponding function
#' evaluations. Must contain a column named `f` with function values. The other
#' columns are treated as design points.
#' @param sigma2 vector of variances used for the Gaussian kernels. If `NULL`
#' (the default), sigma2 is calculated by `doit_estimate_sigma`.
#' @return An object of class `doit` used in further calculations.
#' @export
#'
doit_fit = function(design, sigma2=NULL) {
  m     = nrow(design)
  theta = as.matrix(dplyr::select(design, -f))
  dd    = ncol(theta)
  ff    = design$f

  GGfun = function(xx, yy, sigma2) {
    if (is.null(dim(xx))) xx = matrix(xx, nrow=1)
    if (is.null(dim(yy))) yy = matrix(yy, nrow=1)
    stopifnot(ncol(xx) == ncol(yy), ncol(xx) == length(sigma2))
    mat_ = matrix(0, nrow=nrow(xx), ncol=nrow(yy))
    for (ii in seq_along(sigma2)) {    
      mat_ = mat_ - outer(xx[,ii], yy[,ii], '-')^2 / (2 * sigma2[ii])
    }
    return(drop(exp(mat_)))
  }
  
  if(is.null(sigma2)) {
    sigma2 = doit_estimate_sigma(design)
  }

  # calculate parameters
  GG    = GGfun(theta, theta, sigma2)
  GGinv = solve(GG)
  GG2   = sqrt(GG)
  bb    = drop(solve(GG, sqrt(ff)))
  ee    = drop(1/diag(GGinv) * GGinv %*% sqrt(ff))
  d_ij  = tcrossprod(bb) * GG2
  sum_d_ij = sum(d_ij)

  ans = list(GGfun=GGfun, theta=theta, sigma2=sigma2, 
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
  gg_    = GGfun(theta_eval, theta, sigma2)
  ggbb2_ = drop(gg_ %*% bb)^2
  ee_    = ggbb2_ / (sqrt(prod(pi*sigma2)) * bGG2b_)
  vv_    = ggbb2_ * drop(1 - rowSums(gg_ * (gg_ %*% GGinv)))
  return(as_data_frame(cbind(theta_eval, ee=ee_, vv=vv_)))
})


#' Propose a new design point for the DoIt approximation
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return A parameter value.
#' @export
#'
doit_propose_new = function(doit) with(doit, {
  fn = function(r) {
    gg = GGfun(r, theta, sigma2)
    sum(bb * gg)^2 * drop(1 - gg %*% GGinv %*% gg)
  }
  vv =  (sqrt(ff) - ee)^2 / diag(GGinv)
  r0 = theta[which.max(vv), ]
  optim(r0, fn, control=list(fnscale=-1))$par
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
  theta_new  = as.matrix(dplyr::select(design_new, -f))
  theta_up   = rbind(theta, theta_new)
  gg_        = GGfun(theta, theta_new, sigma2)
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
  if (is.null(theta_eval)) theta_eval = theta[, k]
  nu_ij = 0.5 * outer(theta[,k], theta[,k], '+')
  sd_ = sqrt(sigma2[k]/2)
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
  purrr::map_df(1:ncol(doit$theta), ~doit_marginal(doit, ., theta_eval))
}



#' DoIt approximation of posterior expectation
#'
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return Vector of posterior expectations
#' @export
#'
doit_expectation = function(doit) with(doit, {
  nu_ij = purrr::map(colnames(theta), ~0.5 * outer(theta[,.], theta[,.], '+'))
  return(purrr::map_dbl(nu_ij, ~sum(d_ij * .) / sum_d_ij))
})


#' DoIt approximation of posterior variance
#' 
#' @param doit An object of class `doit`, see function `doit_fit`.
#' @return Posterior (co-)variance matrix.
#' @export
#'
doit_variance = function(doit) with(doit, {
  nu_ij = purrr::map(colnames(theta), ~0.5 * outer(theta[,.], theta[,.], '+'))
  dd = ncol(theta)
  kk = expand.grid(ii = 1:dd, jj = 1:dd) 
  kk = dplyr::filter(kk, jj >= ii) 
  kk = dplyr::mutate(kk, v = purrr::map2_dbl(ii, jj, 
              ~ sum(d_ij * nu_ij[[.x]] * nu_ij[[.y]]) / sum_d_ij))
  v_theta = matrix(0, dd, dd)
  for (l in 1:nrow(kk)) {
    v_theta[kk[l,'ii'], kk[l,'jj']] = v_theta[kk[l,'jj'], kk[l,'ii']] = kk[l,'v']
  }
  rownames(v_theta) = colnames(v_theta) = colnames(theta)
  return(v_theta)
})


