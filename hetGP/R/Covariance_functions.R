###############################################################################
### Covariance functions and their derivatives
### Alternatives so far: Squared-exponential and Matern with nu = 5/2 and 3/2
###############################################################################

### 0) Generic covariance function

##' Correlation function of selected type, supporting both isotropic and product forms
##' @param X1 matrix of design locations, one point per row
##' @param X2 matrix of design locations if correlation is calculated between \code{X1} and \code{X2} (otherwise calculated between \code{X1} and itself)
##' @param theta vector of lengthscale parameters (either of size one if isotropic or of size d if anisotropic)
##' @param type one of "\code{Gaussian}", "\code{Matern5_2}", "\code{Matern3_2}"
##' @export
##' @details 
##' Definition of univariate correlation function and hyperparameters:
##' \itemize{
##' \item "\code{Gaussian}": \eqn{c(x, y) = exp(-(x-y)^2/theta)}
##' \item "\code{Matern5_2}": \eqn{c(x, y) = (1+sqrt(5)/theta * abs(x-y) + 5/(3*theta^2)(x-y)^2) * exp(-sqrt(5)*abs(x-y)/theta)}
##' \item "\code{Matern3_2}": \eqn{c(x, y) = (1+sqrt(3)/theta * abs(x-y)) * exp(-sqrt(3)*abs(x-y)/theta)}
##' }
##' Multivariate correlations are product of univariate ones.
##' @useDynLib hetGP
##' @importFrom  Rcpp evalCpp
cov_gen <- function(X1, X2 = NULL, theta, type = c("Gaussian", "Matern5_2", "Matern3_2")){
  type <- match.arg(type)
  if(type == "Gaussian")
    return(cov_Gaussian(X1 = X1, X2 = X2, theta = theta))
  if(type == "Matern5_2")
    return(cov_Matern5_2(X1 = X1, X2 = X2, theta = theta))
  if(type == "Matern3_2")
    return(cov_Matern3_2(X1 = X1, X2 = X2, theta = theta))
}

## ' Partial derivative of covariance function
## ' @details for compatibility and efficiency, the results is to be multiplied by the matrix it is derived from
## ' @param X1 matrix of design locations
## ' @param X2 matrix of design locations if covariance is calculated between X1 and X2
## ' @param theta vector of lengthscale parameters (either of size one if isotropic or of size d if anisotropic)
## ' @param type one of "Gaussian", "Matern5_2"
## ' @param arg argument corresponding to the required partial derivative: "theta_k", "k_theta_g", "X_i_j"
## ' @param ... additional arguments to be passed for the derivative computation
## ' @export
partial_cov_gen <- function(X1, theta, type = "Gaussian", arg, ..., X2 = NULL){
  if(is.null(X2)){
    if(type == "Gaussian"){
      if(arg == "theta_k")
        return(partial_d_C_Gaussian_dtheta_k(X1 = X1, theta = theta, ...))
      if(arg == "k_theta_g")
        return(partial_d_Cg_Gaussian_d_k_theta_g(X1 = X1, theta = theta, ... ))
      if(arg == "X_i_j")
        return(partial_d_C_Gaussian_dX_i_j(X1 = X1, theta = theta, ... ))
    }
    if(type == "Matern5_2"){
      if(arg == "theta_k")
        return(partial_d_C_Matern5_2_dtheta_k(X1 = X1, theta = theta, ...))
      if(arg == "k_theta_g")
        return(partial_d_Cg_Matern5_2_d_k_theta_g(X1 = X1, theta = theta, ... ))
      if(arg == "X_i_j")
        return(partial_d_C_Matern5_2_dX_i_j(X1 = X1, theta = theta, ... ))
    }
    if(type == "Matern3_2"){
      if(arg == "theta_k")
        return(partial_d_C_Matern3_2_dtheta_k(X1 = X1, theta = theta, ...))
      if(arg == "k_theta_g")
        return(partial_d_Cg_Matern3_2_d_k_theta_g(X1 = X1, theta = theta, ... ))
    }
  }else{
    if(type == "Gaussian"){
      if(arg == "theta_k")
        return(partial_d_k_Gaussian_dtheta_k(X1 = X1, X2 = X2, theta = theta, ...))
      if(arg == "k_theta_g")
        return(partial_d_kg_Gaussian_d_k_theta_g(X1 = X1, X2 = X2, theta = theta, ... ))
      if(arg == "X_i_j")
        return(partial_d_k_Gaussian_dX_i_j(X1 = X1, X2 = X2, theta = theta, ... ))
    }
    if(type == "Matern5_2"){
      if(arg == "theta_k")
        return(partial_d_k_Matern5_2_dtheta_k(X1 = X1, X2 = X2, theta = theta, ...))
      if(arg == "k_theta_g")
        return(partial_d_kg_Matern5_2_d_k_theta_g(X1 = X1, X2 = X2, theta = theta, ... ))
      if(arg == "X_i_j")
        return(partial_d_k_Matern5_2_dX_i_j(X1 = X1, X2 = X2, theta = theta, ... ))
    }
    if(type == "Matern3_2"){
      if(arg == "theta_k")
        return(partial_d_k_Matern3_2_dtheta_k(X1 = X1, X2 = X2, theta = theta, ...))
      if(arg == "k_theta_g")
        return(partial_d_kg_Matern3_2_d_k_theta_g(X1 = X1, X2 = X2, theta = theta, ... ))
    }
  }
}


### A) Gaussian covariance

## Gaussian covariance function : K = exp(-(X1 - X2)^2/theta)
## @param X1 matrix of design locations
## @param X2 matrix of design locations if covariance is calculated between X1 and X2
## @param theta vector of lengthscale parameters (either of size one if isotropic or of size d if anisotropic)
cov_Gaussian <- function(X1, X2 = NULL, theta){
  if(length(theta) == 1){
    return(exp(-distance_cpp(X1, X2)/theta))
  }
  if(is.null(X2)){
    return(exp(-distance_cpp(X1 * rep(1/sqrt(theta), rep(nrow(X1),length(theta))))))
  }
  return(exp(-distance_cpp(X1 * rep(1/sqrt(theta), rep(nrow(X1),length(theta))), X2 * rep(1/sqrt(theta), rep(nrow(X2),length(theta))))))
}

## Partial derivative of the covariance matrix with respect to theta[k] (to be multiplied by the covariance matrix)
partial_d_C_Gaussian_dtheta_k <- function(X1, theta){
  return(distance_cpp(X1)/theta^2)
}

## Partial derivative of the covariance vector with respect to theta[k] (to be multiplied by the covariance vector)
partial_d_k_Gaussian_dtheta_k <- function(X1, X2, theta){
  return(distance_cpp(X1, X2)/theta^2)
}

## Partial derivative of the covariance matrix of the noise process with respect to k_theta_g (to be multiplied by the covariance matrix)
partial_d_Cg_Gaussian_d_k_theta_g <- function(X1, theta, k_theta_g){
  ## 1-dimensional/isotropic case
  if(length(theta) == 1)
    return(distance_cpp(X1)/(theta*k_theta_g^2))

  return(distance_cpp(X1 * rep(1/sqrt(theta), rep(nrow(X1),length(theta))))/k_theta_g^2)
}

## Partial derivative of the covariance vector of the noise process with respect to k_theta_g (to be multiplied by the covariance vector)
partial_d_kg_Gaussian_d_k_theta_g <- function(X1, X2, theta, k_theta_g){
  ## 1-dimensional/isotropic case
  if(length(theta) == 1)
    return(distance_cpp(X1, X2)/(theta*k_theta_g^2))

  return(distance_cpp(X1 * rep(1/sqrt(theta), rep(nrow(X1),length(theta))), X2 * rep(1/sqrt(theta), rep(nrow(X2),length(theta))))/k_theta_g^2)
}


## Derivative with respect to X[i,j]. Useful for pseudo inputs, to be multiplied by the covariance matrix
## @param i1 row
## @param i2 column
partial_d_C_Gaussian_dX_i_j <- function(X1, theta, i1, i2){
  tmp <- partial_d_dist_dX_i1_i2(X1, i1, i2)

  ## 1-dimensional/isotropic case
  if(length(theta) == 1)
    return(tmp/theta)

  return(tmp / theta[i2])
}

## Derivative with respect to X[i,j]. Useful for pseudo inputs, to be multiplied by the covariance matrix
## @param i1 row
## @param i2 column
partial_d_k_Gaussian_dX_i_j <- function(X1, X2, theta, i1, i2){
  tmp <- partial_d_dist_dX1_i1_i2_X2(X1, X2, i1, i2)

  ## 1-dimensional/isotropic case
  if(length(theta) == 1)
    return(tmp/theta)

  return(tmp / theta[i2])
}


### B) Matern 5/2 covariance

## Matern 5/2 covariance function (tensor product)
## @param X1 matrix of design locations
## @param X2 matrix of design locations if covariance is calculated between X1 and X2
## @param theta vector of lengthscale parameters (either of size one if isotropic or of size d if anisotropic)
cov_Matern5_2 <- function(X1, X2 = NULL, theta){
  if(length(theta) == 1){
    if(is.null(X2)){
      return(matern5_2_1args(X1/theta))
    }else{
      return(matern5_2_2args(X1/theta, X2/theta))
    }
  }
  
  if(is.null(X2)){
    return(matern5_2_1args(X1 * rep(1/theta, rep(nrow(X1),length(theta)))))
  }else{
    return(matern5_2_2args(X1 * rep(1/theta, rep(nrow(X1),length(theta))), X2 * rep(1/theta, rep(nrow(X2),length(theta)))))
  }
}

## Partial derivative of the covariance matrix with respect to theta[k] (to be multiplied by the covariance matrix)
partial_d_C_Matern5_2_dtheta_k <- function(X1, theta){
    if(ncol(X1) == 1){
      tmp <- d_matern5_2_1args_theta_k(X1 = X1, theta = theta)
    }else{
      tmp <- d_matern5_2_1args_theta_k_iso(X1 = X1, theta = theta)
    }
    return(tmp)
}

## Partial derivative of the covariance vector with respect to theta[k] (to be multiplied by the covariance vector)
partial_d_k_Matern5_2_dtheta_k <- function(X1, X2, theta){
    tmp <- d_matern5_2_2args_theta_k_iso(X1 = X1, X2 = X2, theta = theta)
    return(tmp)
}

## Partial derivative of the covariance matrix of the noise process with respect to k_theta_g (to be multiplied by the covariance matrix)
partial_d_Cg_Matern5_2_d_k_theta_g <- function(X1, theta, k_theta_g){
  if(length(theta) == 1) return(d_matern5_2_1args_kthetag(X1/theta, k_theta_g))
  return(d_matern5_2_1args_kthetag(X1 * rep(1/theta, rep(nrow(X1),length(theta))), k_theta_g))
}

## Partial derivative of the covariance vector of the noise process with respect to k_theta_g (to be multiplied by the covariance vector)
partial_d_kg_Matern5_2_d_k_theta_g <- function(X1, X2, theta, k_theta_g){
  if(length(theta) == 1) return(d_matern5_2_2args_kthetag(X1/theta, X2/theta, k_theta_g))
  return(d_matern5_2_2args_kthetag(X1 * rep(1/theta, rep(nrow(X1),length(theta))), X2 * rep(1/theta, rep(nrow(X2),length(theta))), k_theta_g))
}

## Derivative with respect to X[i,j]. Useful for pseudo inputs, to be multiplied by the covariance matrix
## @param i1 row
## @param i2 column
partial_d_C_Matern5_2_dX_i_j <- function(X1, theta, i1, i2){


  ## 1-dimensional/isotropic case
  if(length(theta) == 1){
    tmp <- partial_d_dist_abs_dX_i1_i2(X1/theta, i1, i2)
    return(tmp/theta)
  }

  tmp <- partial_d_dist_abs_dX_i1_i2(X1/theta[i2], i1, i2)
  return(tmp / theta[i2])
}

## Derivative with respect to X[i,j]. Useful for pseudo inputs, to be multiplied by the covariance matrix
## @param i1 row
## @param i2 column
## @param theta lengthscales
partial_d_k_Matern5_2_dX_i_j <- function(X1, X2, theta, i1, i2){

  ## 1-dimensional/isotropic case
  if(length(theta) == 1){
    tmp <- partial_d_dist_abs_dX1_i1_i2_X2(X1/theta, X2/theta, i1, i2)
    return(tmp/theta)
  }

  tmp <- partial_d_dist_abs_dX1_i1_i2_X2(X1/theta[i2], X2/theta[i2], i1, i2)
  return(tmp / theta[i2])
}


### C) Matern 3/2 covariance

## Matern 3/2 covariance function (tensor product)
## @param X1 matrix of design locations
## @param X2 matrix of design locations if covariance is calculated between X1 and X2
## @param theta vector of lengthscale parameters (either of size one if isotropic or of size d if anisotropic)
cov_Matern3_2 <- function(X1, X2 = NULL, theta){
  if(length(theta) == 1){
    if(is.null(X2)){
      return(matern3_2_1args(X1/theta))
    }else{
      return(matern3_2_2args(X1/theta, X2/theta))
    }
  }
  
  if(is.null(X2)){
    return(matern3_2_1args(X1 * rep(1/theta, rep(nrow(X1),length(theta)))))
  }else{
    return(matern3_2_2args(X1 * rep(1/theta, rep(nrow(X1),length(theta))), X2 * rep(1/theta, rep(nrow(X2),length(theta)))))
  }
}

## Partial derivative of the covariance matrix with respect to theta[k] (to be multiplied by the covariance matrix) 
partial_d_C_Matern3_2_dtheta_k <- function(X1, theta){
  if(ncol(X1) == 1) return(d_matern3_2_1args_theta_k(X1 = X1, theta = theta))
  return(d_matern3_2_1args_theta_k_iso(X1 = X1, theta = theta))
}

## Partial derivative of the covariance vector with respect to theta[k] (to be multiplied by the covariance vector)
partial_d_k_Matern3_2_dtheta_k <- function(X1, X2, theta){
  tmp <- d_matern3_2_2args_theta_k_iso(X1 = X1, X2 = X2, theta = theta)
  return(tmp)
}

## Partial derivative of the covariance matrix of the noise process with respect to k_theta_g (to be multiplied by the covariance matrix)
partial_d_Cg_Matern3_2_d_k_theta_g <- function(X1, theta, k_theta_g){
  if(length(theta) == 1) return(d_matern3_2_1args_kthetag(X1/theta, k_theta_g))
  return(d_matern3_2_1args_kthetag(X1 * rep(1/theta, rep(nrow(X1),length(theta))), k_theta_g))
}

## Partial derivative of the covariance vector of the noise process with respect to k_theta_g (to be multiplied by the covariance vector)
partial_d_kg_Matern3_2_d_k_theta_g <- function(X1, X2, theta, k_theta_g){
  if(length(theta) == 1) return(d_matern3_2_2args_kthetag(X1/theta, X2/theta, k_theta_g))
  return(d_matern3_2_2args_kthetag(X1 * rep(1/theta, rep(nrow(X1),length(theta))), X2 * rep(1/theta, rep(nrow(X2),length(theta))), k_theta_g))
}


