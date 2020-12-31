#' Hessian of the log-likelihood with respect to lengthscales hyperparameters
#' Works for homGP and hetGP models from the hetGP package for now. 
#' @param model homGP model
#' @return A matrix giving the Hessian of the GP loglikelihood.
#' @importFrom numDeriv hessian
#' @export
logLikHessian <- function(model){
  hess <- hessian(func = logLikH, x = model$theta, g = model$g, X0 = model$X0, Z0 = model$Z0,
                  Z = model$Z, mult = model$mult, beta0 = model$beta0,
                  Delta = model$Delta, k_theta_g = model$k_theta_g, theta_g = model$theta_g,
                  logN = model$logN)
  if(!isSymmetric(hess)) hess <- (hess + t(hess))/2
  return(hess)
}


#' CI on Eigenvalues via Monte Carlo/GP
#' @param model A homGP model
#' @param B Monte Carlo iterates
#' @return A list with elements ci giving 95% quantile based intervals, as well as eigen_draws giving the raw draws
#' @examples 
#' ################################################################################
#' ## Example of uncertainty quantification on C estimate
#' ################################################################################
#' library(hetGP); library(lhs)
#' set.seed(42)
#' 
#' nvar <- 2
#' n <- 20
#' nits <- 20
#' 
#' # theta gives the subspace direction
#' f <- function(x, theta, nugget = 1e-6){
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
#' model <- mleHomGP(design, response, known = list(beta0 =  0))
#' 
#' res <- C_GP_ci(model)
#' 
#' plot(c(1, 2), log(c(mean(res$eigen_draws[,1]), mean(res$eigen_draws[,2]))),
#'   ylim = range(log(res$eigen_draws)), ylab = "Eigenvalue", xlab = "Index")
#'   segments(1, log(res$ci[1,1]), 1, log(res$ci[2,1]))
#'   segments(2, log(res$ci[1,2]), 2, log(res$ci[2,2]))
#' 
#' @export
#' @importFrom stats quantile
C_GP_ci <- function(model, B = 100) {
  
  ll_func <- function(logtheta) logLikH(X0 = model$X0, Z0 = model$Z0, Z = model$Z, mult = model$mult, theta = exp(logtheta), g = model$g, 
                                        covtype = model$covtype, logN = TRUE)
  hess <- hessian(ll_func, log(model$theta))
  
  if(!isSymmetric(hess)) hess <- (hess + t(hess))/2
  
  SIG_THETA <- -solve(hess)
  #t_draws <- rmvnorm(n = B, mean = model$theta, sigma = SIG_THETA)
  R <- chol(SIG_THETA)
  #R <- matrix(0, ncol = 2, nrow = 2)
  sn_draws <- matrix(rnorm(B*ncol(SIG_THETA)), ncol = B)
  t_draws <- t(t(R) %*% sn_draws) 
  lmodeltheta <- log(model$theta)
  t_draws <- t(apply(t_draws, 1, function(x) x + lmodeltheta))
  eigen_draws <- matrix(NA, nrow = B, ncol = length(model$theta))
  
  for (b in 1:B) {
    modeli <- model
    modeli$theta <- exp(t_draws[b,])
    #TODO: rebuild is doing weird stuff.
    modeli <- rebuild(modeli)
    C <- C_GP(modeli)
    
    if (sum(is.nan(C$mat)) > 0) {
      warning("NaN Encountered in C_GP.")
    } else {
      eigen_draws[b,] <- eigen(C$mat)$values
    }
  }
  
  return(list(ci = apply(eigen_draws, 2, quantile, probs = c(0.025, 0.975)),
              eigen_draws = eigen_draws))
}
