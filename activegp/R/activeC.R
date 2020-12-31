#' @title C matrix closed form expression for a GP.
#' @description 
#' Computes the integral (uniform measure) over the [0,1] box domain of the outer product of the gradients of a Gaussian process. 
#' The corresponding matrix is the C matrix central in the active subspace methodology.
#' @param model \code{homGP} or \code{hetGP} GP model, see \code{\link[hetGP]{hetGP-package}} containing, e.g.,
#'  a vector of \code{theta}s, type of covariance \code{ct}, an inverse covariance matrix \code{Ki},
#' a design matrix \code{X0}, and response vector \code{Z0}.
#' @return a \code{const_C} object with elements
#' \itemize{
#' \item \code{model}: GP model provided;
#' \item \code{mat}: C matrix estimated;
#' \item \code{Wij}: list of W matrices, of size number of variables;
#' \item \code{ct}: covariance type (1 for "Gaussian", 2 for "Matern3_2", 3 for "Matern5_2").
#' }
#' @references 
#' N. Wycoff, M. Binois, S. Wild (2019+), Sequential Learning of Active Subspaces, preprint.\cr
#' 
#' P. Constantine (2015), Active Subspaces, Philadelphia, PA: SIAM.
#' @export
#' @useDynLib activegp
#' @importFrom  Rcpp evalCpp
#' @import hetGP
#' @seealso \code{\link[activegp]{print.const_C}}, \code{\link[activegp]{plot.const_C}}
#' @examples 
#' ################################################################################
#' ### Active subspace of a Gaussian process
#' ################################################################################
#' library(hetGP); library(lhs)
#' set.seed(42)
#' 
#' nvar <- 2
#' n <- 100
#' 
#' # theta gives the subspace direction
#' f <- function(x, theta, nugget = 1e-3){
#'   if(is.null(dim(x))) x <- matrix(x, 1)
#'   xact <- cos(theta) * x[,1] - sin(theta) * x[,2]
#'   return(hetGP::f1d(xact) + rnorm(n = nrow(x), sd = rep(nugget, nrow(x))))
#' }
#' 
#' theta_dir <- pi/6
#' act_dir <- c(cos(theta_dir), -sin(theta_dir))
#' 
#' # Create design of experiments and initial GP model
#' design <- X <- matrix(signif(maximinLHS(n, nvar), 2), ncol = nvar)
#' response <- Y <- apply(design, 1, f, theta = theta_dir)
#' model <- mleHomGP(design, response, known = list(beta0 = 0))
#' 
#' C_hat <- C_GP(model)
#' 
#' # Subspace distance to true subspace:
#' print(subspace_dist(C_hat, matrix(act_dir, nrow = nvar), r = 1))
#' plot(design %*% eigen(C_hat$mat)$vectors[,1], response, 
#'   main = "Projection along estimated active direction")
#' plot(design %*% eigen(C_hat$mat)$vectors[,2], response, 
#'   main = "Projection along estimated inactive direction")
#'   
#' # For other plots:
#' # par(mfrow = c(1, 3)) # uncomment to have all plots together
#' plot(C_hat)
#' # par(mfrow = c(1, 1)) # restore graphical window
#' 
C_GP <- function(model){
  # Extract model design and response 
  design <- model$X0
  response <- model$Z0
  
  # isotropic case:
  if(length(model$theta) < ncol(model$X0)) model$theta <- rep(model$theta, ncol(model$X0))
  
  # Discover covariance type
  if (model$covtype == "Gaussian") {
    M_num <- 1
    ct <- 1
    theta <- sqrt(model$theta/2)
  } else if (model$covtype == "Matern3_2") {
    M_num <- 3
    ct <- 2
    theta <- model$theta
  } else if (model$covtype == "Matern5_2") {
    M_num <- 5/3
    ct <- 3
    theta <- model$theta
  } else {
    stop("Unsupported Covariance type in model passed to C_GP.")
  }
  
  nvar <- ncol(design)
  if (nrow(design) != length(response)) {
    stop("microbenchmark")
  }
  if(max(design) > 1 + sqrt(.Machine$double.eps)|| min(design) < 0 - sqrt(.Machine$double.eps)) {
      warning("Designs are supposed to be in [0,1]^d; please rescale.\n Extreme values of ", min(design), ", ", max(design), " detected in design. ")
  }
  
  # Create a const_C object to return
  C <- list()
  C$mat <- matrix(0, nvar, nvar)
  class(C) <- "const_C"
  C$Wij <- lapply(1:nvar, function(i) list())
  C$model <- model
  # C$response <- response
  # C$design <- design
  C$ct <- ct
  
  # Populate its elements
  # Ki <- C$model$Ki 
  Kir <- crossprod(model$Ki, response)
  for(i in 1:nvar)
    for(j in i:nvar){
      Wij <- W_kappa_ij(design = design, theta = theta, i1 = i - 1, i2 = j - 1, ct = ct)
      
      M <- M_num/theta[i]^2 * (i == j) - sum(model$Ki * Wij) + crossprod(Kir, Wij %*% Kir) 
      C$mat[i,j] <- C$mat[j,i] <- M
      C$Wij[[i]][[j]] <- Wij
    }
  return(C)
}

if(!isGeneric("update")) {
  setGeneric(name = "update",
             def = function(object, ...) standardGeneric("update")
  )
}

#' C update with new observations
#'
#' Update Constantine's C with new point(s) for a GP
#' @param object A const_C object, the result of a call to the C_GP function. 
#' @param Xnew matrix (one point per row) corresponding to the new designs
#' @param Znew vector of size \code{nrow(Xnew)} for the new responses at \code{Xnew}
#' @param ... not used (for consistency of update method)
#' @return The updated const_C object originally provided. 
#' @seealso \code{\link[activegp]{C_GP}} to generate const_C objects from \code{\link[hetGP]{mleHomGP}} objects; \code{\link[activegp]{update_C2}} for an update using faster expressions.  
#' @export
#' @useDynLib activegp
#' @importFrom  Rcpp evalCpp
#' @examples
#' \donttest{ 
#' ################################################################################
#' ### Active subspace of a Gaussian process
#' ################################################################################
#' library(hetGP); library(lhs)
#' set.seed(42)
#' 
#' nvar <- 2
#' n <- 100
#' 
#' # theta gives the subspace direction
#' f <- function(x, theta, nugget = 1e-3){
#'   if(is.null(dim(x))) x <- matrix(x, 1)
#'   xact <- cos(theta) * x[,1] - sin(theta) * x[,2]
#'   return(hetGP::f1d(xact) + 
#'     rnorm(n = nrow(x), sd = rep(nugget, nrow(x))))
#' }
#' 
#' theta_dir <- pi/6
#' act_dir <- c(cos(theta_dir), -sin(theta_dir))
#' 
#' # Create design of experiments and initial GP model
#' design <- X <- matrix(signif(maximinLHS(n, nvar), 2), ncol = nvar)
#' response <- Y <- apply(design, 1, f, theta = theta_dir)
#' model <- mleHomGP(design, response, known = list(beta0 = 0))
#' 
#' C_hat <- C_GP(model)
#' 
#' print(C_hat)
#' print(subspace_dist(C_hat, matrix(act_dir, nrow = nvar), r = 1))
#' 
#' # New designs
#' Xnew <- matrix(runif(2), 1)
#' Znew <- f(Xnew, theta_dir)
#' 
#' C_new <- update(C_hat, Xnew, Znew)
#' print(C_new)
#' subspace_dist(C_new, matrix(act_dir, nrow = nvar), r = 1)
#' }
update.const_C <- function(object, Xnew, Znew, ...){
  # Evaluate new quantities
  if(is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  n0 <- nrow(object$model$X0) # to identify replicates
  nvar <- ncol(Xnew)
  
  # Update ancillary quantities.
  object$model <- update(object$model, Xnew = Xnew, Znew = Znew, maxit = 0)
  
  theta <- object$model$theta
  # isotropic case:
  if(length(theta) < ncol(object$model$X0)) theta <- rep(theta, nvar)
  if(object$ct == 1){
    theta <- sqrt(theta/2)
    M_num <- 1
  }else{
    if(object$ct == 2) M_num <- 3 else M_num <- 5/3
  } 
  
  n <- nrow(object$model$X0)
  Kir <- crossprod(object$model$Ki, object$model$Z0)
  
  for(i in 1:nvar) {
    for(j in i:nvar) {
      # If xnew is a replicate, n = n0 and Wijs are unchanged
      if(n > n0){
        Wij <- rbind(cbind(object$Wij[[i]][[j]], matrix(NA, n0, n - n0)), matrix(NA, n - n0, n))
        W_kappa_ij_up(W = Wij, design = object$model$X0, theta, i - 1, j - 1, start = n0, ct = object$ct)
        object$Wij[[i]][[j]] <- Wij
      }else{
        Wij <- object$Wij[[i]][[j]]
      }
      
      # Cov(dYi(X),dYj(X))
      M <- M_num/theta[i]^2 * (i == j) - sum(object$model$Ki * Wij) + crossprod(Kir, Wij) %*% Kir 
      object$mat[i,j] <- object$mat[j,i] <- M
    }
  }
  return(object)
  
}

#' Update Constantine's C, using update formula
#'
#' @param C A const_C object, the result of a call to \code{\link[activegp]{C_GP}}.
#' @param xnew The new design point
#' @param ynew The new response
#' @importFrom stats update predict
#' @references 
#' N. Wycoff, M. Binois, S. Wild (2019+), Sequential Learning of Active Subspaces, preprint.\cr
#' @return Updated C matrix, a const_C object.
#' @export
update_C2 <- function(C, xnew, ynew){
  if(is.null(nrow(xnew))) xnew <- matrix(xnew, nrow = 1)
  nvar <- ncol(xnew)
  
  Cup <- C$mat
  
  kn1 <- cov_gen(xnew, C$model$X0, theta = C$model$theta, type = C$model$covtype)
  
  # for shorter expressions
  if (C$ct == 1) {
    theta <- sqrt(C$model$theta/2)
  } else {
    theta <- C$model$theta
  }
  
  new_lambda <- predict(object = C$model, x = xnew, nugs.only = TRUE)$nugs/C$model$nu_hat
  vn <- drop(1 - kn1 %*% tcrossprod(C$model$Ki, kn1)) + new_lambda + C$model$eps
  
  # precomputations
  Kikn <- tcrossprod(C$model$Ki, kn1)
  gn <- - Kikn / vn
  Kiyn <- C$model$Ki %*% C$model$Z0 # Ki yn
  gyn <- crossprod(gn, C$model$Z0)
  
  for(i in 1:nvar) {
    for(j in i:nvar){
      wa <- drop(W_kappa_ij2(C$model$X0, xnew, theta = theta, i - 1, j - 1, ct = C$ct))  # w(X, xnew)
      wb <- drop(W_kappa_ij2(xnew, rbind(C$model$X0, xnew), theta = theta, i - 1, j - 1, ct = C$ct))  # c(w(xnew, X), w(xnew, xnew))
      w <-  wb[length(wb)]# w(xnew, xnew)
      wb <- wb[-length(wb)]
      kntKiWij <- crossprod(Kikn, C$Wij[[i]][[j]])
      
      tmp <- - crossprod(wa + wb, gn)
      tmp <- tmp - (gyn + ynew/vn) * (kntKiWij %*% Kiyn + crossprod(Kiyn, C$Wij[[i]][[j]] %*% Kikn))
      tmp <- tmp + (gyn + ynew/vn) * crossprod(wa + wb, Kiyn + gn * ynew - gn * drop(kn1 %*% Kiyn))
      tmp <- tmp + ((gyn + ynew/vn)^2 - 1/vn) * (w + kntKiWij %*% Kikn)
      
      Cup[i, j] <- Cup[j, i] <- C$mat[i, j] + tmp 
    }
  }
  return(Cup)
}

#' Expected variance of trace of C 
#' 
#' @param C A const_C object, the result of a call to \code{\link[activegp]{C_GP}}.
#' @param xnew The new design point
#' @param grad If \code{FALSE}, calculate variance of trace after update. If \code{TRUE}, returns the gradient.
#' @return A real number giving the expected variance of the trace of C given the current design.
#' @export
#' @references 
#' N. Wycoff, M. Binois, S. Wild (2019+), Sequential Learning of Active Subspaces, preprint.\cr
#' @examples 
#' \donttest{
#' ################################################################################
#' ### Variance of trace criterion landscape
#' ################################################################################
#'     library(hetGP)
#'     set.seed(42)
#'     nvar <- 2
#'     n <- 20
#' 
#'     # theta gives the subspace direction
#'     f <- function(x, theta = pi/6, nugget = 1e-6){
#'      if(is.null(dim(x))) x <- matrix(x, 1)
#'      xact <- cos(theta) * x[,1] - sin(theta) * x[,2]
#'      return(hetGP::f1d(xact) + 
#'        rnorm(n = nrow(x), sd = rep(nugget, nrow(x))))
#'     }
#' 
#'     design <- matrix(signif(runif(nvar*n), 2), ncol = nvar)
#'     response <- apply(design, 1, f)
#'     model <- mleHomGP(design, response, lower = rep(1e-4, nvar),
#'                       upper = rep(0.5,nvar), known = list(g = 1e-4))
#'                       
#'     C_hat <- C_GP(model)
#' 
#'     ngrid <- 101
#'     xgrid <- seq(0, 1,, ngrid)
#'     Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#'     filled.contour(matrix(f(Xgrid), ngrid))
#' 
#'     Ctr_grid <- apply(Xgrid, 1, C_tr, C = C_hat)
#'     filled.contour(matrix(Ctr_grid, ngrid), color.palette = terrain.colors,
#'                    plot.axes = {axis(1); axis(2); points(design, pch = 20)})
#' }
C_tr <- function(C, xnew, grad = FALSE){
  if(is.null(nrow(xnew))) xnew <- matrix(xnew, nrow = 1)
  nvar <- ncol(xnew)
  #for(i in 1:nvar) {
  #  ret <- get_betagamma(C, xnew, i, i, kn1, Kikn, Kiyn, vn, grad = grad)
  #  beta <- beta + ret$beta
  #  gamma <- gamma + ret$gamma
  #  if (grad) {
  #    dbeta <- dbeta + ret$dbeta
  #    dgamma <- dgamma + ret$dgamma
  #  }
  #}
  #if (grad) {
  #  return(2 * dbeta * beta + 4 * dgamma * gamma)
  #} else {
  #  return(drop(beta^2 + 2*gamma^2))
  #}
  ret <- get_betagamma(C, xnew, grad)
  if (grad) {
    return(sapply(1:nvar, function(d) 2 * sum(diag(ret$dBETA[,,d])) * sum(diag(ret$BETA)) + 4 * sum(diag(ret$dGAMMA[,,d])) * sum(diag(ret$GAMMA))))
  } else {
    return(drop(sum(diag(ret$BETA))^2 + 2*sum(diag(ret$GAMMA))^2))
  }
}

#' Element-wise Cn+1 variance
#'
#' @param C A const_C object, the result of a call to \code{\link[activegp]{C_GP}}.
#' @param xnew The new design point
#' @param grad If \code{FALSE}, calculate variance of update. If \code{TRUE}, returns the gradient.
#' @return A real number giving the expected elementwise variance of C given the current design.
#' @references 
#' N. Wycoff, M. Binois, S. Wild (2019+), Sequential Learning of Active Subspaces, preprint.\cr
#' @export
#' @examples 
#' ################################################################################
#' ### Norm of the variance of C criterion landscape
#' ################################################################################
#' library(hetGP)
#' set.seed(42)
#' nvar <- 2
#' n <- 20
#' 
#' # theta gives the subspace direction
#' f <- function(x, theta = pi/6, nugget = 1e-6){
#'  if(is.null(dim(x))) x <- matrix(x, 1)
#'  xact <- cos(theta) * x[,1] - sin(theta) * x[,2]
#'  return(hetGP::f1d(xact) 
#'    + rnorm(n = nrow(x), sd = rep(nugget, nrow(x))))
#' }
#' 
#' design <- matrix(signif(runif(nvar*n), 2), ncol = nvar)
#' response <- apply(design, 1, f)
#' model <- mleHomGP(design, response, lower = rep(1e-4, nvar),
#'                   upper = rep(0.5,nvar), known = list(g = 1e-4))
#'                   
#' C_hat <- C_GP(model)
#' 
#' ngrid <- 51
#' xgrid <- seq(0, 1,, ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' filled.contour(matrix(f(Xgrid), ngrid))
#' 
#' cvar_crit <- function(C, xnew){
#'  return(sqrt(sum(C_var(C, xnew)^2)))
#' }
#' 
#' Cvar_grid <- apply(Xgrid, 1, cvar_crit, C = C_hat)
#' filled.contour(matrix(Cvar_grid, ngrid), color.palette = terrain.colors,
#'                plot.axes = {axis(1); axis(2); points(design, pch = 20)})
C_var <- function(C, xnew, grad = FALSE){
  if(is.null(nrow(xnew))) xnew <- matrix(xnew, nrow = 1)
  nvar <- ncol(xnew)
  #for(i in 1:nvar) {
  #  for(j in i:nvar) {
  #    #beta <- (crossprod(Kiyn, C$Wij[[i]][[j]] %*% Kikn) + kntKiWij %*% Kiyn - crossprod(wa + wb, Kiyn))/sqrt(vn)
  #    #gamma <- (w + kntKiWij %*% Kikn - crossprod(wa + wb, Kikn))/vn
  #    ret <- get_betagamma(C, xnew, i, j, kn1, Kikn, Kiyn, vn, grad = grad)
  #    if (grad) {
  #      dCvar[i, j,] <- dCvar[j, i,] <- drop(2*ret$beta*ret$dbeta + 4*ret$gamma*ret$dgamma)
  #    } 
  #    Cvar[i, j] <- Cvar[j, i] <- drop(ret$beta^2 + 2*ret$gamma^2)
  #  }
  #}
  #if (grad) {
  #  normgrads <- 2*sapply(1:nvar, function(d) t(as.numeric(dCvar[,,d])) %*% as.numeric(Cvar))
  #  return(normgrads)
  #  #return(dCvar)
  #} else {
  #  return(norm(Cvar, 'F')^2)
  #  #return(Cvar)
  #}
  ret <- get_betagamma(C, xnew, grad = grad)
  Cvar <- ret$BETA^2 + 2 * ret$GAMMA^2
  if (grad) {
    #return(2*sapply(1:nvar, function(d) t(as.numeric(dCvar[,,d])) %*% as.numeric(Cvar)))
    return(2*sapply(1:nvar, function(d) t(as.numeric(2*ret$BETA*ret$dBETA[,,d] + 4 * ret$GAMMA*ret$dGAMMA[,,d])) %*% as.numeric(Cvar)))
  } else {
    return(norm(Cvar, 'F')^2)
  }
}

#' Alternative Variance of Update
#'
#' Defined as E[(C - E[C])^2], where A^2 = AA (not elementwise multiplication).
#'
#' @param C A const_C object, the result of a call to \code{\link[activegp]{C_GP}}.
#' @param xnew The new design point
#' @param grad If \code{FALSE}, calculate variance of update. If \code{TRUE}, returns the gradient.
#' @return A real number giving the expected variance of C defined via matrix multiplication given the current design.
#' @references 
#' N. Wycoff, M. Binois, S. Wild (2019+), Sequential Learning of Active Subspaces, preprint.\cr
#' @export
#' @examples 
#' ################################################################################
#' ### Norm of the variance of C criterion landscape
#' ################################################################################
#' library(hetGP)
#' set.seed(42)
#' nvar <- 2
#' n <- 20
#' 
#' # theta gives the subspace direction
#' f <- function(x, theta = pi/6, nugget = 1e-6){
#'  if(is.null(dim(x))) x <- matrix(x, 1)
#'  xact <- cos(theta) * x[,1] - sin(theta) * x[,2]
#'  return(hetGP::f1d(xact) + rnorm(n = nrow(x), sd = rep(nugget, nrow(x))))
#' }
#' 
#' design <- matrix(signif(runif(nvar*n), 2), ncol = nvar)
#' response <- apply(design, 1, f)
#' model <- mleHomGP(design, response, lower = rep(1e-4, nvar),
#'                   upper = rep(0.5,nvar), known = list(g = 1e-4))
#'                   
#' C_hat <- C_GP(model)
#' 
#' ngrid <- 51
#' xgrid <- seq(0, 1,, ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' filled.contour(matrix(f(Xgrid), ngrid))
#' 
#' cvar_crit <- function(C, xnew){
#'  return(sqrt(sum(C_var(C, xnew)^2)))
#' }
#' 
#' Cvar_grid <- apply(Xgrid, 1, cvar_crit, C = C_hat)
#' filled.contour(matrix(Cvar_grid, ngrid), color.palette = terrain.colors,
#'                plot.axes = {axis(1); axis(2); points(design, pch = 20)})
C_var2 <- function(C, xnew, grad = FALSE){
  if(is.null(nrow(xnew))) xnew <- matrix(xnew, nrow = 1)
  nvar <- ncol(xnew)
  ret <- get_betagamma(C, xnew, grad = grad)
  Cvar2 <- ret$BETA %*% ret$BETA + 2 * ret$GAMMA %*% ret$GAMMA
  if (grad) {
    Cvar2d <- array(NA, dim = c(nvar, nvar, nvar))
    for (p in 1:nvar) {
      Cvar2d[,,p] <- ret$dBETA[,,p] %*% ret$BETA + ret$BETA %*% ret$dBETA[,,p] + 2 * (ret$dGAMMA[,,p] %*% ret$GAMMA + ret$GAMMA %*% ret$dGAMMA[,,p])
    }
    normgrads <- 2*sapply(1:nvar, function(d) t(as.numeric(Cvar2d[,,d])) %*% as.numeric(Cvar2))
    return(normgrads)
  } else {
    return(norm(Cvar2, 'F')^2)
  }
}

#' Quantities for Acquisition Functions
#'
#' Create a single element of the BETA/GAMMA matrix. Used to compute acquisition functions and their gradients.
#'
#' @param C A const_C object, the result of a call to C_GP
#' @param xnew The new design point
#' @param grad If \code{FALSE}, calculate beta and gamma only. If \code{TRUE}, calculate their gradient too.
#' @return If \code{grad == FALSE}, A numeric vector of length 2, whose first element of beta_ij and the second gamma_ij. 
#' If \code{grad == TRUE}, a list with 3 numeric vector elements, the first giving the gradient for beta_ij, and the second for gamma_ij,
#' and the third is the same vector as would have been returned if grad was \code{FALSE}: simply the values of beta and gamma.
#' @keywords internal
get_betagamma <- function(C, xnew, grad = FALSE) {
  if(is.null(nrow(xnew))) xnew <- matrix(xnew, nrow = 1)
  nvar <- ncol(xnew)
  kn1 <- cov_gen(xnew, C$model$X0, theta = C$model$theta, type = C$model$covtype)
  Ki <- C$model$Ki
  
  new_lambda <- predict(object = C$model, x = xnew, nugs.only = TRUE)$nugs/C$model$nu_hat
  vn <- drop(1 - kn1 %*% tcrossprod(Ki, kn1)) + new_lambda + C$model$eps
  
  # precomputations
  Kikn <- tcrossprod(Ki, kn1)
  Kiyn <- Ki %*% C$model$Z0 # Ki yn
  
  # h <- 1e-6
  dkn1 <- matrix(NA, nvar, nrow(Ki))
  for (ih in 1:nvar) {
    # xh <- rep(0, nvar)
    # xh[ih] <- h
    # kn1h <- cov_gen(xnew + xh, C$model$X0, theta = C$model$theta, type = C$model$covtype)
    # dkn1 <- rbind(dkn1, (kn1h - kn1) / h)
    dkn1[ih,] <- d1(C$model$X0[, ih], x = xnew[ih], sigma = C$model$theta[ih], type = C$model$covtype) * kn1
  }
  
  dvn <- t(-2 * Ki %*% t(kn1)) %*% t(dkn1)
  
  Cvar <- C$mat
  
  if (C$ct == 1) {
    theta <- sqrt(C$model$theta/2)
  } else {
    theta <- C$model$theta
  }
  
  if (grad) {
    dBETA <- dGAMMA <- array(NA, dim = c(nvar, nvar, nvar))
  }
  BETA <- GAMMA <- matrix(NA, nrow = nvar, ncol = nvar)
  
  for (i in 1:nvar) {
    for (j in i:nvar) {
      wa <- drop(W_kappa_ij2(C$model$X0, xnew, theta = theta, i - 1, j - 1, ct = C$ct))  # w(X, xnew)
      wb <- drop(W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct))  # c(w(xnew, X), w(xnew, xnew))
      Wij <- C$Wij[[i]][[j]]
      kntKiWij <- crossprod(Kikn, Wij)
      WijKiKn <- Wij %*% Kikn
      
      w <- drop(W_kappa_ij2(xnew, xnew, theta = theta, i - 1, j - 1, ct = C$ct))
      
      betanum <- drop((crossprod(Kiyn, WijKiKn) + kntKiWij %*% Kiyn - crossprod(wa + wb, Kiyn)))
      beta <- betanum / sqrt(vn)
      gammanum <- drop((w + kntKiWij %*% Kikn - crossprod(wa + wb, Kikn)))
      gamma <- gammanum / vn
      
      BETA[i,j] <- BETA[j,i] <- beta
      GAMMA[i,j] <- GAMMA[j,i] <- gamma
      
      if (grad) {
        # Get W's derivative via finite differencing
        dWa <- grad_W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)
        dWb <- grad_W_kappa_ij2_w2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)
        
        dwa <- drop(grad_W_kappa_ij2(xnew, xnew, theta = theta, i - 1, j - 1, ct = C$ct))
        dwb <- drop(grad_W_kappa_ij2_w2(xnew, xnew, theta = theta, i - 1, j - 1, ct = C$ct))
        dw <- dwa + dwb
        #h <- 1e-6
        #dw <- rep(NA, nvar)
        #for (ih in 1:nvar) {
        #  xh <- rep(0, nvar)
        #  xh[ih] <- h
        #  wh <- drop(W_kappa_ij2(xnew + xh, xnew + xh, theta = theta, i - 1, j - 1, ct = C$ct))
        #  dw[ih] <- (wh - w) / h
        #}
        
        AA <- Ki %*% (Wij %*% Kiyn + t(crossprod(Kiyn, Wij)))
        dbeta <- ((t(dkn1 %*% AA) - t((dWa + dWb) %*% Kiyn)) * sqrt(vn) - 
                    drop(betanum * 0.5 * vn^(-0.5)) * dvn) / (vn)
        BB <- t(Ki %*% (WijKiKn + t(kntKiWij))) %*% t(dkn1)
        CC <- t((dWa + dWb) %*% Kikn) + t(dkn1 %*% Ki %*% (wa + wb))
        dgamma <- ((t(dw) + BB - CC) * vn - drop(gammanum * dvn)) / (vn)^2
        dBETA[i,j,] <- dBETA[j,i,] <- dbeta
        dGAMMA[i,j,] <- dGAMMA[j,i,] <- dgamma
      } 
    }
  }
  
  if (grad) {
    return(list(BETA = BETA, GAMMA = GAMMA, dBETA = dBETA, dGAMMA = dGAMMA))
  } else {
    return(list(BETA = BETA, GAMMA = GAMMA))
  }
}

#' Print const_C objects
#' @param x A const_C object, the result of a call to C_GP
#' @param ... Additional parameters. Not used. 
#' @export
print.const_C <- function(x, ...) {
  cts <- c("Gaussian", "Matern3_2", "Matern5_2")
  cat(paste(cts[x$ct], "kernel GP Estimate of Constantine's C:\n"))
  print(x$mat)
}

#' Plot const_C objectc
#' @param x A const_C object, the result of a call to C_GP
#' @param output one of \code{"image"} (image of the C matrix), \code{"logvals"} (log-eigen values), 
#' \code{"projfn"} projected function on first eigen vector or all plots at once (default).
#' @param ... Additional parameters. Not used. 
#' @importFrom graphics image plot
#' @export
plot.const_C <- function(x, output = c("all", "matrix", "logvals", "projfn"), ...) {
  output <- match.arg(output)
  if(output %in% c("all", "matrix")) image(x$mat, main = "C matrix values heatmap")
  if(output %in% c("all", "logvals")) plot(log(eigen(x$mat)$values), main = "log eigen values of C", xlab = "index", ylab = "")
  if(output %in% c("all", "projfn")) plot(x$model$X0 %*% eigen(x$mat)$vectors[,1], x$model$Y0, xlab = "First AS direction", ylab = "Function values")

}


#' Active subspace for second order linear model
#' @param design A matrix of design points, one in each row, in [-1,1]^d
#' @param response A vector of observations at each design point.
#' @return A matrix corresponding to the active subspace C matrix. 
#' @importFrom stats lm reformulate
#' @keywords internal
#' @export
#' @examples
#' set.seed(42) 
#' A <- matrix(c(1, -1, 0, -1, 2, -1.5, 0, -1.5, 4), nrow = 3, byrow = TRUE)
#' b <- c(1, 4, 9)
#'
#' # Quadratic function
#' ftest <- function(x, sd = 1e-6){
#'    if(is.null(dim(x))) x <- matrix(x, nrow = 1)
#'    return(3 + drop(diag(x %*% A %*% t(x)) + x %*% b) + 
#'      rnorm(nrow(x), sd = sd))
#' }
#' 
#' ntrain <- 10000
#' design <- 2 * matrix(runif(ntrain * 3), ntrain) - 1
#' response <- ftest(design)
#' 
#' C_hat <- C_Q(design, response)
#' 
#' plot(design %*% eigen(C_hat)$vectors[,1], response)
#' 
#' # Test 
#' gfun <- function(x){2 * A %*% t(x) + matrix(b, nrow = nrow(A), ncol = nrow(x))}
#' grads <- gfun(design)
#' C_MC <- tcrossprod(grads)/ntrain
#' C_true <- 4/3 * A %*% A + tcrossprod(b)
#' subspace_dist(eigen(C_MC)$vectors[,1:2], eigen(C_true)$vectors[,1:2]) 
C_Q <- function(design, response){
  d <- ncol(design)
  
  # create second order formula
  formulatmp <- reformulate(c(".^2", paste0("I(X", 1:ncol(design), "^2)")), response = "y") 
  
  model <- lm(formulatmp, data = data.frame(design, y = response))
  b <- model$coefficients[2:(d + 1)]
  A <- matrix(0, d, d)
  A[lower.tri(A)] <- 1/2*model$coefficients[(2*d + 2):length(model$coefficients)]
  A <- (A + t(A))
  diag(A) <- model$coefficients[(d+2):(2 * d + 1)]
  
  return(4/3 * A %*% A + tcrossprod(b))
}

