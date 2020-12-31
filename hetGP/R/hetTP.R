#  dmt(Z, mean = rep(0,133), S = (nu - 2)/nu*(model$sigma2 * cov_gen(X, theta = model$theta) + diag(model$g, 133)), df = model$nu, log = T)
loglik_HomTP <- function(X0, Z0, Z, mult, theta, g, nu, sigma2, beta0 = 0,
                         covtype = "Gaussian", eps = sqrt(.Machine$double.eps), env = NULL){
  n <- nrow(X0)
  N <- length(Z)
  
  C <- cov_gen(X1 = X0, theta = theta, type = covtype)
  # Temporarily store Cholesky transform of K in Ki
  Ki <- chol(add_diag(sigma2 * C, eps + g/mult))
  ldetKi <- - 2 * sum(log(diag(Ki))) # log determinant from Cholesky
  Ki <- chol2inv(Ki)
  
  if(!is.null(env)){
    env$C <- C
    env$Ki <- Ki
  } 
  
  psi_0 <- drop(crossprod(Z0 - beta0, Ki) %*% (Z0 - beta0))
  
  psi <- (crossprod(Z - beta0) - crossprod((Z0 - beta0) * mult, Z0 - beta0))/g + psi_0
  
  return(-N/2 * log((nu - 2) * pi) + ldetKi/2 - (N - n)/2 * log(g) - 1/2 * sum(log(mult)) + lgamma((nu + N)/2) - lgamma(nu/2) - (nu + N)/2 * log(1 + psi/(nu - 2)))
}


# derivative of log-likelihood for logLikHom with respect to theta with all observations (Gaussian kernel)
## Model: noisy observations with unknown homoskedastic noise
## K = tau^2 * C + g * I
# X0  design matrix (no replicates)
# Z0 averaged observations
# Z observations vector (all observations)
# mult number of replicates per unique design point
# theta vector of lengthscale hyperparameters (or one for isotropy)
# g noise variance for the process
# tau2 scale
# nu degrees of freedom
# beta0 trend
## @return gradient with respect to theta and g
dlogLikHomTP <- function(X0, Z0, Z, mult, theta, g, nu, sigma2, beta0 = 0, covtype = "Gaussian", eps = sqrt(.Machine$double.eps),
                         components = c("theta", "g", "nu", "sigma2"), env = NULL){
  N <- length(Z)
  n <- nrow(X0)
  
  if(!is.null(env)){
    C <- env$C
    Ki <- env$Ki
  }else{
    C <- cov_gen(X1 = X0, theta = theta, type = covtype)
    Ki <- chol2inv(chol(sigma2*C + diag(eps + g / mult)))
  }
  
  Z0 <- Z0 - beta0
  Z <- Z - beta0
  
  KiZ0 <- Ki %*% Z0 ## to avoid recomputing  
  
  psi_0 <- drop(crossprod(Z0, KiZ0))
  
  psi <- (crossprod(Z) - crossprod(Z0 * mult, Z0))/g + psi_0
  
  
  tmp1 <- tmp2 <- tmp3 <- tmp4 <- NULL
  
  # First component, derivative with respect to theta
  if("theta" %in% components){
    tmp1 <- rep(NA, length(theta))
    if(length(theta)==1){
      dC_dthetak <- partial_cov_gen(X1 = X0, theta = theta, type = covtype, arg = "theta_k") * C
      tmp1 <- sigma2 * (nu + N)/(2 * (nu + psi - 2)) * crossprod(KiZ0, dC_dthetak) %*% KiZ0  - sigma2/2 * trace_sym(Ki, dC_dthetak)
    }else{
      for(i in 1:length(theta)){
        dC_dthetak <- partial_cov_gen(X1 = X0[,i, drop = F], theta = theta[i], type = covtype, arg = "theta_k") * C
        tmp1[i] <- sigma2 * (nu + N)/(2 * (nu + psi - 2)) * crossprod(KiZ0, dC_dthetak) %*% KiZ0 - sigma2/2 * trace_sym(Ki, dC_dthetak)
      }
    } 
  }
  
  if("nu" %in% components) 
    tmp2 <- -N / (2*(nu - 2)) + 1/2*digamma((nu + N)/2) - 1/2*digamma(nu/2) - 1/2 * log(1 + psi/(nu - 2)) + psi * (nu + N)/(2*(nu - 2)^2 + 2 * psi * (nu - 2))
  
  if("sigma2" %in% components) 
    tmp3 <- (nu + N)/(2 * (nu + psi - 2)) * crossprod(KiZ0, C) %*% KiZ0 -1/2 * trace_sym(Ki, C)
  
  # Fourth component derivative with respect to g
  if("g" %in% components)
    tmp4 <- (nu + N)/(2 * (nu + psi - 2)) * ((crossprod(Z) - crossprod(Z0 * mult, Z0))/g^2 + sum(KiZ0^2/mult)) - (N - n)/ (2*g) - 1/2 * sum(diag(Ki)/mult)
  
  return(c(tmp1, tmp2, tmp3, tmp4))
}



##' Student-t process regression under homoskedastic noise based on maximum likelihood estimation of the 
##' hyperparameters. This function is enhanced to deal with replicated observations.
##' @title Student-T process modeling with homoskedastic noise
##' @param X matrix of all designs, one per row, or list with elements:
##' \itemize{
##'   \item \code{X0} matrix of unique design locations, one point per row
##'   \item \code{Z0} vector of averaged observations, of length \code{nrow(X0)}
##'   \item \code{mult} number of replicates at designs in \code{X0}, of length \code{nrow(X0)}
##' } 
##' @param Z vector of all observations. If using a list with \code{X}, \code{Z} has to be ordered with respect to \code{X0}, and of length \code{sum(mult)}
##' @param lower,upper bounds for the \code{theta} parameter (see \code{\link[hetGP]{cov_gen}} for the exact parameterization).
##' In the multivariate case, it is possible to give vectors for bounds (resp. scalars) for anisotropy (resp. isotropy) 
##' @param known optional list of known parameters, e.g., \code{beta0} (default to \code{0}), \code{theta}, \code{g}, \code{sigma2} or \code{nu}
##' @param covtype covariance kernel type, either 'Gaussian', 'Matern5_2' or 'Matern3_2', see \code{\link[hetGP]{cov_gen}}
##' @param noiseControl list with element,
##' \itemize{
##' \item \code{g_bound}, vector providing minimal and maximal noise variance
##' \item \code{sigma2_bounds}, vector providing minimal and maximal signal variance
##' \item \code{nu_bounds}, vector providing minimal and maximal values for the degrees of freedom. 
##' The mininal value has to be stricly greater than 2. If the mle optimization gives a large value, e.g., 30,
##' considering a GP with \code{\link[hetGP]{mleHomGP}} may be better. 
##' } 
##' @param init list specifying starting values for MLE optimization, with elements:
##' \itemize{
##'  \item \code{theta_init} initial value of the theta parameters to be optimized over (default to 10\% of the range determined with \code{lower} and \code{upper})
##'  \item \code{g_init} initial value of the nugget parameter to be optimized over (based on the variance at replicates if there are any, else 10\% of the variance)
##'  \item \code{sigma2} initial value of the variance paramter (default to \code{1})
##'  \item \code{nu} initial value of the degrees of freedom parameter (default to \code{3})
##' }
##' @param maxit maximum number of iteration for L-BFGS-B of \code{\link[stats]{optim}}
##' @param eps jitter used in the inversion of the covariance matrix for numerical stability
##' @param settings list with argument \code{return.Ki}, to include the inverse covariance matrix in the object for further use (e.g., prediction).
##' Arguments \code{factr} (default to 1e9) and \code{pgtol} are available to be passed to \code{control} for L-BFGS-B in \code{\link[stats]{optim}}.
##' @return a list which is given the S3 class "\code{homGP}", with elements:
##' \itemize{
##' \item \code{theta}: maximum likelihood estimate of the lengthscale parameter(s),
##' \item \code{g}: maximum likelihood estimate of the nugget variance,
##' \item \code{trendtype}: either "\code{SK}" if \code{beta0} is given, else "\code{OK}" 
##' \item \code{beta0}: estimated trend unless given in input,
##' \item \code{sigma2}:  maximum likelihood estimate of the scale variance,
##' \item \code{nu2}:  maximum likelihood estimate of the degrees of freedom parameter,
##' \item \code{ll}: log-likelihood value,
##' \item \code{X0}, \code{Z0}, \code{Z}, \code{mult}, \code{eps}, \code{covtype}: values given in input,
##' \item \code{call}: user call of the function
##' \item \code{used_args}: list with arguments provided in the call
##' \item \code{nit_opt}, \code{msg}: \code{counts} and {msg} returned by \code{\link[stats]{optim}}
##' \item \code{Ki}, inverse covariance matrix (if \code{return.Ki} is \code{TRUE} in \code{settings})
##' \item \code{time}: time to train the model, in seconds.
##' 
##'}
##' @details
##' The global covariance matrix of the model is parameterized as \code{K = sigma2 * C + g * diag(1/mult)},
##' with \code{C} the correlation matrix between unique designs, depending on the family of kernel used (see \code{\link[hetGP]{cov_gen}} for available choices).
##' 
##' It is generally recommended to use \code{\link[hetGP]{find_reps}} to pre-process the data, to rescale the inputs to the unit cube and to normalize the outputs.
##' 
##' @seealso \code{\link[hetGP]{predict.homTP}} for predictions.
##' \code{summary} and \code{plot} functions are available as well. 
##' @references 
##' M. Binois, Robert B. Gramacy, M. Ludkovski (2017+), Practical heteroskedastic Gaussian process modeling for large simulation experiments, arXiv preprint arXiv:1611.05902. \cr \cr
##' 
##' A. Shah, A. Wilson, Z. Ghahramani (2014), Student-t processes as alternatives to Gaussian processes, Artificial Intelligence and Statistics, 877--885. \cr \cr
##' 
##' M. Chung, M. Binois, RB Gramacy, DJ Moquin, AP Smith, AM Smith (2018+). Parameter and Uncertainty Estimation for Dynamical Systems Using Surrogate Stochastic Processes. Preprint available on arXiv:1802.00852. 
##' @examples
##' ##------------------------------------------------------------
##' ## Example 1: Homoskedastic Student-t modeling on the motorcycle data
##' ##------------------------------------------------------------
##' set.seed(32)
##' 
##' ## motorcycle data
##' library(MASS)
##' X <- matrix(mcycle$times, ncol = 1)
##' Z <- mcycle$accel
##' plot(X, Z, ylim = c(-160, 90), ylab = 'acceleration', xlab = "time")
##' 
##' noiseControl = list(g_bounds = c(1e-3, 1e4))
##' model <- mleHomTP(X = X, Z = Z, lower = 0.01, upper = 100, noiseControl = noiseControl)
##' summary(model)
##'   
##' ## Display averaged observations
##' points(model$X0, model$Z0, pch = 20) 
##' xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 
##' preds <- predict(x = xgrid, object =  model)
##' 
##' ## Display mean prediction
##' lines(xgrid, preds$mean, col = 'red', lwd = 2)
##' ## Display 95% confidence intervals
##' lines(xgrid, preds$mean + sqrt(preds$sd2) * qt(0.05, df = model$nu + nrow(X)), col = 2, lty = 2)
##' lines(xgrid, preds$mean + sqrt(preds$sd2) * qt(0.95, df = model$nu + nrow(X)), col = 2, lty = 2)
##' ## Display 95% prediction intervals
##' lines(xgrid, preds$mean + sqrt(preds$sd2 + preds$nugs) * qt(0.05, df = model$nu + nrow(X)), 
##'   col = 3, lty = 2)
##' lines(xgrid, preds$mean + sqrt(preds$sd2 + preds$nugs) * qt(0.95, df = model$nu + nrow(X)), 
##'   col = 3, lty = 2)
##' @export
mleHomTP <- function(X, Z, lower = NULL, upper = NULL, known = list(beta0 = 0),
                     noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), 1e4),
                                         nu_bounds = c(2 + 1e-3, 30),
                                         sigma2_bounds = c(sqrt(.Machine$double.eps), 1e4)),
                     init = list(nu = 3), covtype = c("Gaussian", "Matern5_2", "Matern3_2"), maxit = 100,
                     eps = sqrt(.Machine$double.eps),
                     settings = list(return.Ki = TRUE, factr = 1e9)){
  
  if(typeof(X)=="list"){
    X0 <- X$X0
    Z0 <- X$Z0
    mult <- X$mult
    if(sum(mult) != length(Z)) stop("Length(Z) should be equal to sum(mult)")
    if(is.null(dim(X0))) X0 <- matrix(X0, ncol = 1)
    if(length(Z0) != nrow(X0)) stop("Dimension mismatch between Z0 and X0")
  }else{
    if(is.null(dim(X))) X <- matrix(X, ncol = 1)
    if(nrow(X) != length(Z)) stop("Dimension mismatch between Z and X")
    elem <- find_reps(X, Z, return.Zlist = F)
    X0 <- elem$X0
    Z0 <- elem$Z0
    Z <- elem$Z
    mult <- elem$mult
  }
  
  covtype <- match.arg(covtype)
  
  if(is.null(lower) || is.null(upper)){
    auto_thetas <- auto_bounds(X = X0, covtype = covtype)
    if(is.null(lower)) lower <- auto_thetas$lower
    if(is.null(upper)) upper <- auto_thetas$upper
    if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- sqrt(upper * lower)
  }
  
  if(length(lower) != length(upper)) stop("upper and lower should have the same size")
  
  ## Save time to train model
  tic <- proc.time()[3]
  
  if(is.null(settings$return.Ki)) settings$return.Ki <- TRUE
  if(is.null(settings$factr)) settings$factr <- 1e9
  if(is.null(noiseControl$g_bounds)) noiseControl$g_bounds <- c(sqrt(.Machine$double.eps), 1e4)
  if(is.null(noiseControl$sigma2_bounds)) noiseControl$sigma2_bounds <- c(sqrt(.Machine$double.eps), 1e4)
  if(is.null(noiseControl$nu_bounds)) noiseControl$nu_bounds <- c(2 + 1e-3, 30)
  
  if(is.null(known$beta0)) known$beta0 <- 0
  beta0 <- known$beta0
  
  N <- length(Z)
  n <- nrow(X0)
  
  if(is.null(n))
    stop("X0 should be a matrix. \n")
  
  if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- 0.1 * lower + 0.9 * upper
  if(is.null(known$g) && is.null(init$g)){
    if(any(mult > 2)) init$g <- mean((fast_tUY2(mult, (Z - rep(Z0, times = mult))^2)/mult)[which(mult > 2)]) else init$g <- 0.1 * var(Z0)
  } 
  if(is.null(known$nu) && is.null(init$nu)) init$nu <- 3
  if(is.null(known$sigma2) && is.null(init$sigma2)) init$sigma2 <- var(Z0)

  parinit <- c(init$theta, init$nu, init$sigma2, init$g)
  
  # trendtype <- 'OK'
  # if(!is.null(beta0))
  #   trendtype <- 'SK'
  
  ## General definition of fn and gr
  fn <- function(par, X0, Z0, Z, mult, beta0, theta, sigma2, nu, g, env){
    idx <- 1 # to store the first non used element of par
    
    if(is.null(theta)){
      idx <- length(init$theta)
      theta <- par[1:idx]
      idx <- idx + 1
    }
    if(is.null(nu)){
      nu <- par[idx]
      idx <- idx + 1
    }
    if(is.null(sigma2)){
      sigma2 <- par[idx]
      idx <- idx + 1
    }
    if(is.null(g)){
      g <- par[idx]
    }
    
    loglik <- loglik_HomTP(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, nu = nu, sigma2 = sigma2,
                           beta0 = beta0, covtype = covtype, eps = eps, env = env)
    
    if(!is.null(env) && !is.na(loglik)){
      if(is.null(env$max_loglik) || loglik > env$max_loglik){
        env$max_loglik <- loglik
        env$arg_max <- par
      }
    } 
    
    return(loglik)
  }
  
  gr <- function(par, X0, Z0, Z, mult, beta0, theta, sigma2, nu, g, env){
    idx <- 1
    components <- NULL
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
      components <- "theta"
    }
    if(is.null(nu)){
      nu <- par[idx]
      idx <- idx + 1
      components <- c(components, "nu")
    }
    if(is.null(sigma2)){
      sigma2 <- par[idx]
      idx <- idx + 1
      components <- c(components, "sigma2")
    }
    if(is.null(g)){
      g <- par[idx]
      components <- c(components, "g")
    }
    
    return(dlogLikHomTP(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, nu = nu, sigma2 = sigma2,
                        beta0 = beta0, covtype = covtype, eps = eps, components = components, env = env))
  }
  
  ## All known
  if(!is.null(known$g) && !is.null(known[["theta"]]) && !is.null(known$nu) && !is.null(known$sigma2)){
    theta_out <- known[["theta"]]
    nu_out <- known$nu
    sigma2_out <- known$sigma2
    g_out <- known$g
    out <- list(value = loglik_HomTP(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta_out, g = g_out, nu = nu_out, sigma2 = sigma2_out, 
                                     beta0 = beta0, covtype = covtype, eps = eps),
                message = "All hyperparameters given", counts = 0, time = proc.time()[3] - tic)
  }else{
    parinit <- lowerOpt <- upperOpt <- NULL
    if(is.null(known[["theta"]])){
      parinit <- init$theta
      lowerOpt <- c(lower)
      upperOpt <- c(upper)
    }
    if(is.null(known$nu)){
      parinit <- c(parinit, init$nu)
      lowerOpt <- c(lowerOpt, noiseControl$nu_bounds[1])
      upperOpt <- c(upperOpt, noiseControl$nu_bounds[2])
    }
    if(is.null(known$sigma2)){
      parinit <- c(parinit, init$sigma2)
      lowerOpt <- c(lowerOpt, noiseControl$sigma2_bounds[1])
      upperOpt <- c(upperOpt, noiseControl$sigma2_bounds[2])
    }
    if(is.null(known$g)){
      parinit <- c(parinit, init$g)
      lowerOpt <- c(lowerOpt, noiseControl$g_bounds[1])
      upperOpt <- c(upperOpt, noiseControl$g_bounds[2])
    }
    
    
    # environment storing values from log to pass to dlog
    envtmp <- environment()
    out <- try(optim(par = parinit, fn = fn, gr = gr, method = "L-BFGS-B", lower = lowerOpt, upper = upperOpt, theta = known[["theta"]], 
                 nu = known$nu, sigma2 = known$sigma2, g = known$g, env = envtmp,
                 X0 = X0, Z0 = Z0, Z = Z, mult = mult, beta0 = beta0,
                 control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol)))
    ## Catch errors when at least one likelihood evaluation worked
    if(class(out) == "try-error")
      out <- list(par = envtmp$arg_max, value = envtmp$max_loglik, counts = NA,
                  message = "Optimization stopped due to NAs, use best value so far")
    
    ## Post-processing
    idx <- 1
    if(is.null(known[["theta"]])){
      theta_out <- out$par[1:length(init$theta)]
      idx <- idx + length(init$theta)
    }else{
      theta_out <- known$theta
    }
    if(is.null(known$nu)){
      nu_out <-  out$par[idx]
      idx <- idx + 1
    }else{
      nu_out <- known$nu
    }
    if(is.null(known$sigma2)){
      sigma2_out <-  out$par[idx]
      idx <- idx + 1
    }else{
      sigma2_out <- known$sigma2
    }
    if(is.null(known$g)){
      g_out <-  out$par[idx]
    }else{
      g_out <- known$g
    }
    
  }
  
  Ki <- chol2inv(chol(add_diag(sigma2_out * cov_gen(X1 = X0, theta = theta_out, type = covtype), eps + g_out/mult)))
  
  # if(is.null(beta0))
  #   beta0 <- drop(colSums(Ki) %*% Z0 / sum(Ki))
  
  psi_0 <- drop(crossprod(Z0 - beta0, Ki) %*% (Z0 - beta0))
  
  psi <- drop(crossprod(Z - beta0) - crossprod((Z0 - beta0) * mult, Z0 - beta0))/g_out + psi_0
  
  res <- list(theta = theta_out, g = g_out, nu = nu_out, sigma2 = sigma2_out, mult = mult,
              # trendtype = trendtype,
              ll = out$value, nit_opt = out$counts,
              psi = psi,
              beta0 = beta0,  covtype = covtype, msg = out$message, eps = eps,
              X0 = X0, Z = Z, mult = mult, Z0 = Z0, call = match.call(),
              used_args = list(lower = lower, upper = upper, known = known, noiseControl = noiseControl), 
              time = proc.time()[3] - tic)
  
  if(settings$return.Ki) res <- c(res, list(Ki = Ki))
  
  class(res) <- "homTP"
  return(res)
  
}

##' @method summary homTP
##' @export
summary.homTP <- function(object,...){
  ans <- object
  class(ans) <- "summary.homTP"
  ans
}

##' @export
print.summary.homTP <- function(x, ...){
  
  cat("N = ", length(x$Z), " n = ", length(x$Z0), " d = ", ncol(x$X0), "\n")
  cat(x$covtype, " covariance lengthscale values: ", x$theta, "\n")
  
  cat("Homoskedastic nugget value: ", x$g, "\n")
  
  cat("Variance/scale hyperparameter: ", x$sigma2, "\n")
  cat("Degree of freedom nu:", x$nu, "\n")
  
  cat("Given constant trend value: ", x$beta0, "\n")
  
  cat("MLE optimization: \n", "Log-likelihood = ", x$ll, "; Nb of evaluations (obj, gradient) by L-BFGS-B: ", x$nit_opt, "; message: ", x$msg, "\n")
}

##' @method print homTP
##' @export
print.homTP <- function(x, ...){
  print(summary(x))
  
}

##' @method plot homTP
##' @export
plot.homTP <- function(x, ...){
  LOOpreds <- LOO_preds(x)
  plot(x$Z, LOOpreds$mean[rep(1:nrow(x$X0), times = x$mult)], xlab = "Observed values", ylab = "Predicted values",
       main = "Leave-one-out predictions")
  arrows(x0 = LOOpreds$mean + sqrt(LOOpreds$sd2) * qt(0.05, df = x$nu + nrow(x$X0)),
         x1 = LOOpreds$mean + sqrt(LOOpreds$sd2) * qt(0.95, df = x$nu + nrow(x$X0)),
         y0 = LOOpreds$mean, length = 0, col = "blue")
  points(x$Z0[which(x$mult > 1)], LOOpreds$mean[which(x$mult > 1)], pch = 20, col = 2)
  abline(a = 0, b = 1, lty = 3)
  legend("topleft", pch = c(1, 20, NA), lty = c(NA, NA, 1), col = c(1, 2, 4),
         legend = c("observations", "averages (if > 1 observation)", "LOO prediction interval"))
}

##' Student-t process predictions using a homoskedastic noise GP object (of class \code{homGP})
##' @param x matrix of designs locations to predict at
##' @param object an object of class \code{homGP}; e.g., as returned by \code{\link[hetGP]{mleHomTP}}
##' @param xprime optional second matrix of predictive locations to obtain the predictive covariance matrix between \code{x} and \code{xprime}
##' @param ... no other argument for this method
##' @return list with elements
##' \itemize{
##' \item \code{mean}: kriging mean;
##' \item \code{sd2}: kriging variance (filtered, e.g. without the nugget value)
##' \item \code{cov}: predictive covariance matrix between \code{x} and \code{xprime}
##' \item \code{nugs}: nugget value at each prediction location
##' }
##' @importFrom MASS ginv
##' @details The full predictive variance corresponds to the sum of \code{sd2} and \code{nugs}.
##' @method predict homTP
##' @export
predict.homTP <- function(object, x, xprime = NULL, ...){
  
  if(is.null(dim(x))){
    x <- matrix(x, nrow = 1)
    if(ncol(x) != ncol(object$X0)) stop("x is not a matrix")
  }
  
  if(!is.null(xprime) && is.null(dim(xprime))){
    xprime <- matrix(xprime, nrow = 1)
    if(ncol(xprime) != ncol(object$X0)) stop("xprime is not a matrix")
  }
  
  if(is.null(object$Ki))
    object$Ki <- chol2inv(chol(object$sigma2 * cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype) + diag(object$g/object$mult + object$eps)))
  
  
  n1 <- length(object$Z)
  kx <- object$sigma2 * cov_gen(X1 = x, X2 = object$X0, theta = object$theta, type = object$covtype)
  
  mean <- as.vector(object$beta0 + kx %*% (object$Ki %*% (object$Z0 - object$beta0)))
  
  # if(object$trendtype == 'SK'){
  sd2 <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * as.vector(object$sigma2 - fast_diag(kx, tcrossprod(object$Ki, kx)))
  # }else{
  # sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)) + (1 - tcrossprod(rowSums(object$Ki), kx))^2/sum(object$Ki))
  # }
  
  ## In case of numerical errors, some sd2 values may become negative
  if(any(sd2 < 0)){
    # object$Ki <- ginv(add_diag(cov_gen(X1 = object$X, theta = object$theta, type = object$covtype), rep(object$g + object$eps, nrow(object$X0))))/object$sigma2
    # mean <- as.vector(object$beta0 + kx %*% (object$Ki %*% (object$Z0 - object$beta0)))
    # 
    # sd2 <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * as.vector(object$sigma2 - fast_diag(kx, tcrossprod(object$Ki, kx)))
    sd2 <- pmax(0, sd2)
    warning("Numerical errors caused some negative predictive variances to be thresholded to zero. Consider using ginv via rebuild.homTP")
  }
  
  if(!is.null(xprime)){
    kxprime <- object$sigma2 * cov_gen(X1 = object$X0, X2 = xprime, theta = object$theta, type = object$covtype)
    # if(object$trendtype == 'SK'){
    if(nrow(x) > nrow(xprime)){
      cov <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * (object$sigma2 * cov_gen(X1 = x, X2 = xprime, theta = object$theta, type = object$covtype) - kx %*% (object$Ki %*% kxprime))
    }else{
      cov <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * (object$sigma2 * cov_gen(X1 = x, X2 = xprime, theta = object$theta, type = object$covtype) - (kx %*% object$Ki) %*% kxprime)
    }
    # }else{
    # cov <- object$sigma2 * cov_gen(X1 = x, X2 = xprime, theta = object$theta, type = object$covtype) - kx %*% tcrossprod(object$Ki, kxprime) + crossprod(1 - tcrossprod(rowSums(object$Ki), kx), 1 - tcrossprod(rowSums(object$Ki), kxprime))/sum(object$Ki)
    # }
  }else{
    cov = NULL
  }
  
  return(list(mean = mean, sd2 = sd2, cov = cov, nugs = rep(object$g, length(mean))))
}


## ' Rebuild inverse covariance matrix of \code{homTP} (e.g., if exported without \code{Ki})
## ' @param object \code{homTP} model without slot \code{Ki} (inverse covariance matrix)
## ' @param robust if \code{TRUE} \code{\link[MASS]{ginv}} is used for matrix inversion, otherwise it is done via Cholesky.
##' @rdname ExpImp
##' @method rebuild homTP
##' @export
rebuild.homTP <- function(object, robust = FALSE){
  if(robust){
    object$Ki <- ginv(add_diag(cov_gen(X1 = object$X, theta = object$theta, type = object$covtype), rep(object$g + object$eps, nrow(object$X0))))/object$sigma2
  }else{
    object$Ki <- chol2inv(chol(object$sigma2 * cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype) + diag(object$g/object$mult + object$eps)))
  }
  
  return(object)
}

###############################################################################
## Part II: Heterogeneous TP with all options for the fit
###############################################################################


# dmt(Z, mean = rep(0,length(Z)), S = (model$nu - 2)/model$nu*(model$sigma2 * cov_gen(X, theta = model$theta, type = "Matern5_2") + diag(rep(model$Lambda, times = model$mult) + 1e-8)), df = model$nu,  log = T)
logLik_HetTP <- function(X0, Z0, Z, mult, Delta, theta, g, nu, sigma2, k_theta_g = NULL, theta_g = NULL, logN = TRUE,
                         beta0 = 0, eps = sqrt(.Machine$double.eps), covtype = "Gaussian",
                         penalty = T, hom_ll = NULL, env = NULL, trace = 0){
  n <- nrow(X0)
  N <- length(Z)
  
  if(is.null(theta_g))
    theta_g <- k_theta_g * theta
  
  Cg <- cov_gen(X1 = X0, theta = theta_g, type = covtype)
  Kg_c <- chol(Cg + diag(eps + g/mult))
  Kgi <- chol2inv(Kg_c)
  
  nmean <- drop(rowSums(Kgi) %*% Delta / sum(Kgi)) ## ordinary kriging mean
  
  M <- Cg %*% (Kgi %*% (Delta - nmean))
  
  
  Lambda <- drop(nmean + M)
  
  if(logN){
    Lambda <- exp(Lambda)
  }
  else{
    Lambda[Lambda <= 0] <- eps
  }
  
  LambdaN <- rep(Lambda, times = mult)
  
  C <- cov_gen(X1 = X0, theta = theta, type = covtype)
  Ki <- chol(add_diag(sigma2 * C, Lambda/mult + eps))
  ldetKi <- - 2 * sum(log(diag(Ki))) # log determinant from Cholesky
  Ki <- chol2inv(Ki)
  
  if(!is.null(env)){
    env$C <- C
    env$Cg <- Cg
    env$Kg_c <- Kg_c
    env$Kgi <- Kgi
    env$ldetKi <- ldetKi
    env$Ki <- Ki
  }
  
  psi_0 <- drop(crossprod(Z0 - beta0, Ki) %*% (Z0 - beta0))
  
  psi <- drop(crossprod((Z - beta0)/LambdaN, Z - beta0) - crossprod((Z0 - beta0) * mult/Lambda, Z0 - beta0) + psi_0)
  
  loglik <- -N/2 * log((nu - 2) * pi) + ldetKi/2 - 1/2 * sum((mult - 1) * log(Lambda) + log(mult)) + lgamma((nu + N)/2) - lgamma(nu/2) - (nu + N)/2 * log(1 + psi/(nu - 2))
  
  if(penalty){
    nu_hat_var <- drop(crossprod(Delta - nmean, Kgi) %*% (Delta - nmean))/length(Delta)
    
    ## To avoid 0 variance, e.g., when Delta = nmean
    if(nu_hat_var < eps) return(loglik)
    
    # if(hardpenalty)
    #   return(loglik + min(0, - n/2 * log(nu_hat_var) - sum(log(diag(Kg_c))) - n/2*log(2*pi) - n/2))
    
    pen <- - n/2 * log(nu_hat_var) - sum(log(diag(Kg_c))) - n/2*log(2*pi) - n/2
    
    if(loglik < hom_ll && pen > 0){
      if(trace> 0) warning("Penalty is desactivated when unpenalized likelihood is lower than its homTP equivalent")
      return(loglik)
    }
    
    return(loglik + pen)
  }
  return(loglik)
}

dlogLikHetTP <- function(X0, Z0, Z, mult, Delta, theta, nu, sigma2, g, k_theta_g = NULL, theta_g = NULL, beta0 = NULL, pX = NULL,
                         logN = TRUE, components = NULL, eps = sqrt(.Machine$double.eps), covtype = "Gaussian", 
                         penalty = T, hom_ll = NULL, env = NULL){
  
  ## Verifications
  if(is.null(k_theta_g) && is.null(theta_g))
    cat("Either k_theta_g or theta_g must be provided \n")
  
  ## Initialisations
  
  # Which terms need to be computed
  if(is.null(components)){
    components <- c(list("theta", "Delta", "g", "nu"))
    if(is.null(k_theta_g)){
      components <- c(components, list("theta_g"))
    }else{
      components <- c(components, list("k_theta_g"))
    }
    if(!is.null(pX))
      components <- c(components, list("pX"))
  }
  
  if(is.null(theta_g))
    theta_g <- k_theta_g * theta
  
  n <- nrow(X0)
  N <- length(Z)
  
  if(!is.null(env)){
    Cg <- env$Cg
    Kg_c <- env$Kg_c
    Kgi <- env$Kgi
  }else{
    Cg <- cov_gen(X1 = X0, theta = theta_g, type = covtype)
    Kg_c <- chol(Cg + diag(eps + g/mult))
    Kgi <- chol2inv(Kg_c)
  }

  # M <- Cg %*% Kgi
  M <- add_diag(Kgi * (-eps - g / mult), rep(1, n))
  
  ## Precomputations for reuse
  
  rSKgi <- rowSums(Kgi)
  sKgi <- sum(Kgi)
  
  nmean <- drop(rSKgi %*% Delta / sKgi) ## ordinary kriging mean
  
  ## Precomputations for reuse
  KgiD <- Kgi %*% (Delta - nmean)
  
  # if(penalty){
  #   nu_hat_var <- drop(crossprod(KgiD, (Delta - nmean)))/length(Delta) 
  #   # To prevent numerical issues when Delta = nmean, giving a positive penalty (or if the penalty is positive)
  #   if(nu_hat_var < eps || (hardpenalty && (- n/2 * log(nu_hat_var) - sum(log(diag(Kg_c))) - n/2*log(2*pi) - n/2) > 0)) penalty <- FALSE
  # }                                          
  
  Lambda <- drop(nmean + M %*% (Delta - nmean))
  if(logN){
    Lambda <- exp(Lambda)
  }
  else{
    Lambda[Lambda <= 0] <- eps
  }
  
  LambdaN <- rep(Lambda, times = mult)
  
  if(!is.null(env)){
    C <- env$C
    Ki <- env$Ki
    ldetKi <- env$ldetKi
  }else{
    C <- cov_gen(X1 = X0, theta = theta, type = covtype)
    Ki <- chol(sigma2 * C + diag(Lambda/mult + eps))
    ldetKi <- - 2 * sum(log(diag(Ki))) # log determinant from Cholesky
    Ki <- chol2inv(Ki)
  }

  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% Z0 / sum(Ki))
  
  ## Precomputations for reuse
  KiZ0 <- Ki %*% (Z0 - beta0)
  rsM <- rowSums(M)
  
  psi_0 <- drop(crossprod(KiZ0, Z0 - beta0))
  
  psi <- drop(crossprod((Z - beta0)/LambdaN, Z - beta0) - crossprod((Z0 - beta0) * mult/Lambda, Z0 - beta0) + psi_0)
  
  if(penalty){
    nu_hat_var <- drop(crossprod(KgiD, (Delta - nmean)))/length(Delta) 
    
    # To prevent numerical issues when Delta = nmean, resulting in divisions by zero
    if(nu_hat_var < eps){
      penalty <- FALSE
    }else{
      loglik <- -N/2 * log(2*pi) - N/2 * log(psi/N)  + 1/2 * ldetKi - 1/2 * sum((mult - 1) * log(Lambda) + log(mult)) - N/2
      pen <- - n/2 * log(nu_hat_var) - sum(log(diag(Kg_c))) - n/2*log(2*pi) - n/2
      if(loglik < hom_ll && pen > 0) penalty <- FALSE
    }
    # if(nu_hat_var < eps || (hardpenalty && (- n/2 * log(nu_hat_var) - sum(log(diag(Kg_c))) - n/2*log(2*pi) - n/2) > 0)) penalty <- FALSE
  }   
  
  dLogL_dtheta <- dLogL_dDelta <- dLogL_dkthetag <- dLogL_dthetag <- dLogL_dg <- dLogL_dnu <- dLogL_dsigma2 <- NULL
  
  # First component, derivative of logL with respect to theta
  if("theta" %in% components){
    dLogL_dtheta <- rep(NA, length(theta))
    for(i in 1:length(theta)){
      if(length(theta) == 1){
        dC_dthetak <- partial_cov_gen(X1 = X0, theta = theta, arg = "theta_k", type = covtype) * C # partial derivative of C with respect to theta
      }else{
        dC_dthetak <- partial_cov_gen(X1 = X0[, i, drop = FALSE], theta = theta[i], arg = "theta_k", type = covtype) * C # partial derivative of C with respect to theta
      }
      
      if("k_theta_g" %in% components){
        if(length(theta) == 1){
          dCg_dthetak <- partial_cov_gen(X1 = X0, theta = k_theta_g * theta, arg = "theta_k", type = covtype) * k_theta_g * Cg # partial derivative of Cg with respect to theta[i]
        }else{
          dCg_dthetak <- partial_cov_gen(X1 = X0[, i, drop = FALSE], theta = k_theta_g * theta[i], arg = "theta_k", type = covtype) * k_theta_g * Cg # partial derivative of Cg with respect to theta[i]
        }
        
        # Derivative Lambda / theta_k (first part)
        dLdtk <- dCg_dthetak %*% KgiD - M %*% (dCg_dthetak %*% KgiD)
        
        # (second part)
        dLdtk <- dLdtk - (1 - rsM) * drop(rSKgi %*% dCg_dthetak %*% (Kgi %*% Delta) * sKgi - rSKgi %*% Delta * (rSKgi %*% dCg_dthetak %*% rSKgi))/sKgi^2
        
        if(logN)
          dLdtk <- dLdtk * Lambda
        
        dK_dthetak <- add_diag(sigma2 * dC_dthetak, drop(dLdtk)/mult) # dK/dtheta[k]
        dLogL_dtheta[i] <- (nu + N)/(2 * (nu + psi - 2)) * (crossprod((Z - beta0)/LambdaN * rep(dLdtk, times = mult), (Z - beta0)/LambdaN) - crossprod((Z0 - beta0)/Lambda * mult * dLdtk, (Z0 - beta0)/Lambda) +
                                                              crossprod(KiZ0, dK_dthetak) %*% KiZ0) - 1/2 * trace_sym(Ki, dK_dthetak)
        dLogL_dtheta[i] <- dLogL_dtheta[i] - 1/2 * sum((mult - 1) * dLdtk/Lambda) # derivative of the sum(a_i - 1)log(lambda_i)
        
        if(penalty)
          dLogL_dtheta[i] <- dLogL_dtheta[i]  + 1/2 * crossprod(KgiD, dCg_dthetak) %*% KgiD / nu_hat_var  - 1/2 * trace_sym(Kgi, dCg_dthetak)
        
      }else{
        dLogL_dtheta[i] <- sigma2 * (nu + N)/(2 * (nu + psi - 2)) * crossprod(KiZ0, dC_dthetak) %*% KiZ0  - sigma2/2 * trace_sym(Ki, dC_dthetak)
      }
    }
  }
  
  # Derivative of logL with respect to Lambda
  if(any(c("Delta", "g", "k_theta_g", "theta_g", "pX") %in% components)){
    dLogLdLambda <- (nu + N)/(2 * (nu + psi - 2)) * ((fast_tUY2(mult, (Z - beta0)^2) - (Z0 - beta0)^2 * mult)/Lambda^2 + KiZ0^2/mult) - (mult - 1)/(2*Lambda) - 1/(2*mult) * diag(Ki)
    
    if(logN){
      dLogLdLambda <- Lambda * dLogLdLambda
    }
  }
  
  
  ## Derivative of Lambda with respect to Delta
  if("Delta" %in% components)
    # dLogL_dDelta <- crossprod(M - tcrossprod(M, matrix(rep(rSKgi / sKgi, ncol(Kgi)), ncol(Kgi))) + matrix(rep(rSKgi / sKgi, nrow(M)), nrow(M), byrow = T), dLogLdLambda) #chain rule
    dLogL_dDelta <- crossprod(M, dLogLdLambda) + rSKgi/sKgi*sum(dLogLdLambda) -  rSKgi / sKgi * sum(crossprod(M, dLogLdLambda)) #chain rule
  
  # Derivative Lambda / k_theta_g
  if("k_theta_g" %in% components){
    dCg_dk <- partial_cov_gen(X1 = X0, theta = theta, k_theta_g = k_theta_g, arg = "k_theta_g", type = covtype) * Cg
    
    dLogL_dkthetag <- dCg_dk %*% KgiD - M %*% (dCg_dk %*% KgiD) -
      (1 - rsM) * drop(rSKgi %*% dCg_dk %*% (Kgi %*% Delta) * sKgi - rSKgi %*% Delta * (rSKgi %*% dCg_dk %*% rSKgi))/sKgi^2
    
    dLogL_dkthetag <- crossprod(dLogL_dkthetag, dLogLdLambda) ## chain rule
  }
  
  # Derivative Lambda / theta_g
  if("theta_g" %in% components){
    dLogL_dthetag <- rep(NA, length(theta_g))
    
    for(i in 1:length(theta_g)){
      
      if(length(theta_g) == 1){
        dCg_dthetagk <- partial_cov_gen(X1 = X0, theta = theta_g, arg = "theta_k", type = covtype) * Cg # partial derivative of Cg with respect to theta
      }else{
        dCg_dthetagk <- partial_cov_gen(X1 = X0[, i, drop = FALSE], theta = theta_g[i], arg = "theta_k", type = covtype) * Cg # partial derivative of Cg with respect to theta
      }
      
      dLogL_dthetag[i] <- crossprod(dCg_dthetagk %*% KgiD - M %*% (dCg_dthetagk %*% KgiD) -
                                      (1 - rsM) * drop(rSKgi %*% dCg_dthetagk %*% (Kgi %*% Delta) * sKgi - rSKgi %*% Delta * (rSKgi %*% dCg_dthetagk %*% rSKgi))/sKgi^2, dLogLdLambda) #chain rule
      # Penalty term
      if(penalty) dLogL_dthetag[i] <- dLogL_dthetag[i] + 1/2 * crossprod(KgiD, dCg_dthetagk) %*% KgiD/nu_hat_var - trace_sym(Kgi, dCg_dthetagk)/2 
    }
  }
  
  if("sigma2" %in% components) 
    dLogL_dsigma2 <- (nu + N)/(2 * (nu + psi - 2)) * crossprod(KiZ0, C) %*% KiZ0 -1/2 * trace_sym(Ki, C)
  
  if("nu" %in% components) 
    dLogL_dnu <- -N / (2*(nu - 2)) + 1/2*digamma((nu + N)/2) - 1/2*digamma(nu/2) - 1/2 * log(1 + psi/(nu - 2)) + psi * (nu + N)/(2*(nu - 2)^2 + 2 * psi * (nu - 2))
  
  ## Derivative Lambda / g
  if("g" %in% components){
    dLogL_dg <- crossprod(-M %*% (KgiD/mult) - (1 - rsM) * drop(Delta %*% (Kgi %*% (rSKgi/mult)) * sKgi - rSKgi %*% Delta * sum(rSKgi^2/mult))/sKgi^2, dLogLdLambda) #chain rule
  }
  
  
  # Additional penalty terms on Delta
  if(penalty){
    if("Delta" %in% components){
      dLogL_dDelta <- dLogL_dDelta - KgiD / nu_hat_var
    }
    if("k_theta_g" %in% components){
      dLogL_dkthetag <- dLogL_dkthetag + 1/2 * crossprod(KgiD, dCg_dk) %*% KgiD / nu_hat_var - trace_sym(Kgi, dCg_dk)/2 
    }
    if("g" %in% components){
      dLogL_dg <- dLogL_dg + 1/2 * crossprod(KgiD/mult, KgiD) / nu_hat_var - sum(diag(Kgi)/mult)/2 
    }
    
  }
  
  return(c(dLogL_dtheta,
           dLogL_dDelta,
           dLogL_dkthetag,
           dLogL_dthetag,
           dLogL_dsigma2,
           dLogL_dnu,
           dLogL_dg))
  
}


##' @title Student-t process modeling with heteroskedastic noise
##' @description 
##' Student-t process regression under input dependent noise based on maximum likelihood estimation of the hyperparameters.
##' A GP is used to model latent (log-) variances.
##' This function is enhanced to deal with replicated observations.
##' @param X matrix of all designs, one per row, or list with elements:
##' \itemize{
##'   \item \code{X0} matrix of unique design locations, one point per row
##'   \item \code{Z0} vector of averaged observations, of length \code{nrow(X0)}
##'   \item \code{mult} number of replicates at designs in \code{X0}, of length \code{nrow(X0)}
##' } 
##' @param Z vector of all observations. If using a list with \code{X}, \code{Z} has to be ordered with respect to \code{X0}, and of length \code{sum(mult)}
##' @param lower,upper bounds for the \code{theta} parameter (see \code{\link[hetGP]{cov_gen}} for the exact parameterization).
##' In the multivariate case, it is possible to give vectors for bounds (resp. scalars) for anisotropy (resp. isotropy)
##' @param noiseControl list with elements related to optimization of the noise process parameters:
##' \itemize{
##' \item \code{g_min}, \code{g_max} minimal and maximal noise to signal ratio (of the mean process)
##' \item \code{lowerDelta}, \code{upperDelta} optional vectors (or scalars) of bounds on \code{Delta}, of length \code{nrow(X0)} (default to \code{rep(eps, nrow(X0))} and \code{rep(noiseControl$g_max, nrow(X0))} resp., or their \code{log}) 
## ' \item lowerpX, upperpX optional vectors of bounds of the input domain if pX is used.
##' \item \code{lowerTheta_g}, \code{upperTheta_g} optional vectors of bounds for the lengthscales of the noise process if \code{linkThetas == 'none'}.
##' Same as for \code{theta} if not provided.
##' \item \code{k_theta_g_bounds} if \code{linkThetas == 'joint'}, vector with minimal and maximal values for \code{k_theta_g} (default to \code{c(1, 100)}). See Details.
##' \item \code{g_bounds} vector for minimal and maximal noise to signal ratios for the noise of the noise process, i.e., the smoothing parameter for the noise process.
##' (default to \code{c(1e-6, 1)}).
##' \item \code{sigma2_bounds}, vector providing minimal and maximal signal variance.
##' \item \code{nu_bounds}, vector providing minimal and maximal values for the degrees of freedom. 
##'}
##' @param settings list for options about the general modeling procedure, with elements:
##' \itemize{
##'   \item \code{linkThetas} defines the relation between lengthscales of the mean and noise processes.
##'   Either \code{'none'}, \code{'joint'}(default) or \code{'constr'}, see Details.
##'   \item \code{logN}, when \code{TRUE} (default), the log-noise process is modeled.
##'   \item \code{initStrategy} one of \code{'simple'}, \code{'residuals'} (default) and \code{'smoothed'} to obtain starting values for \code{Delta}, see Details
##'   \item \code{penalty} when \code{TRUE}, the penalized version of the likelihood is used (i.e., the sum of the log-likelihoods of the mean and variance processes, see References).
## '   \item \code{hardpenalty} is \code{TRUE}, the log-likelihood from the noise GP is taken into account only if negative.
##'   \item \code{checkHom} when \code{TRUE}, if the log-likelihood with a homoskedastic model is better, then return it.
##'   \item \code{trace} optional scalar (default to \code{0}). If positive, tracing information on the fitting process.
##' If \code{1}, information is given about the result of the heterogeneous model optimization.
##' Level \code{2} gives more details. Level {3} additionaly displays all details about initialization of hyperparameters.
##' \item \code{return.matrices} boolean too include the inverse covariance matrix in the object for further use (e.g., prediction).
##' \item Arguments \code{factr} (default to 1e9) and \code{pgtol} are available to be passed to \code{control} for L-BFGS-B in \code{\link[stats]{optim}}.   
##' }
##' @param eps jitter used in the inversion of the covariance matrix for numerical stability
##' @param init,known optional lists of starting values for mle optimization or that should not be optimized over, respectively.
##' Values in \code{known} are not modified, while it can happen to those of \code{init}, see Details. 
##' One can set one or several of the following:
##' \itemize{
##' \item \code{theta} lengthscale parameter(s) for the mean process either one value (isotropic) or a vector (anistropic)
##' \item \code{Delta} vector of nuggets corresponding to each design in \code{X0}, that are smoothed to give \code{Lambda}
##' (as the global covariance matrix depend on \code{Delta} and \code{nu_hat}, it is recommended to also pass values for \code{theta})
##' \item \code{beta0} constant trend of the mean process
##' \item \code{k_theta_g} constant used for link mean and noise processes lengthscales, when \code{settings$linkThetas == 'joint'}
##' \item \code{theta_g} either one value (isotropic) or a vector (anistropic) for lengthscale parameter(s) of the noise process, when \code{settings$linkThetas != 'joint'}
##' \item \code{g} scalar nugget of the noise process
##' \item \code{nu} degree of freedom parameter
##' \item \code{sigma2} scale variance
##' \item \code{g_H} scalar homoskedastic nugget for the initialisation with a \code{\link[hetGP]{mleHomGP}}. See Details.
## '\item pX matrix of fixed pseudo inputs locations of the noise process corresponding to Delta
##' }
##' @param covtype covariance kernel type, either \code{'Gaussian'}, \code{'Matern5_2'} or \code{'Matern3_2'}, see \code{\link[hetGP]{cov_gen}}
##' @param maxit maximum number of iterations for \code{L-BFGS-B} of \code{\link[stats]{optim}} dedicated to maximum likelihood optimization
##' @details
##' 
##' The global covariance matrix of the model is parameterized as \code{K = sigma2 * C + Lambda * diag(1/mult)},
##' with \code{C} the correlation matrix between unique designs, depending on the family of kernel used (see \code{\link[hetGP]{cov_gen}} for available choices).
##' \code{Lambda} is the prediction on the noise level given by a (homoskedastic) GP: \cr
##' \deqn{\Lambda = C_g(C_g + \mathrm{diag}(g/\mathrm{mult}))^{-1} \Delta} \cr
##' with \code{C_g} the correlation matrix between unique designs for this second GP, with lengthscales hyperparameters \code{theta_g} and nugget \code{g}
##' and \code{Delta} the variance level at \code{X0} that are estimated.
##' 
##' It is generally recommended to use \code{\link[hetGP]{find_reps}} to pre-process the data, to rescale the inputs to the unit cube and to normalize the outputs.
##' 
##' The noise process lengthscales can be set in several ways:
##' \itemize{
##' \item using \code{k_theta_g} (\code{settings$linkThetas == 'joint'}), supposed to be greater than one by default. 
##' In this case lengthscales of the noise process are multiples of those of the mean process.
##' \item if \code{settings$linkThetas == 'constr'}, then the lower bound on \code{theta_g} correspond to estimated values of an homoskedastic GP fit.
##' \item else lengthscales between the mean and noise process are independent (both either anisotropic or not).
##' }
##'
##' When no starting nor fixed parameter values are provided with \code{init} or \code{known}, 
##' the initialization process consists of fitting first an homoskedastic model of the data, called \code{modHom}.
##' Unless provided with \code{init$theta}, initial lengthscales are taken at 10\% of the range determined with \code{lower} and \code{upper},
##' while \code{init$g_H} may be use to pass an initial nugget value.
##' The resulting lengthscales provide initial values for \code{theta} (or update them if given in \code{init}). \cr \cr
##' If necessary, a second homoskedastic model, \code{modNugs}, is fitted to the empirical residual variance between the prediction
##'  given by \code{modHom} at \code{X0} and \code{Z} (up to \code{modHom$nu_hat}).
##' Note that when specifying \code{settings$linkThetas == 'joint'}, then this second homoskedastic model has fixed lengthscale parameters.
##' Starting values for \code{theta_g} and \code{g} are extracted from \code{modNugs}.\cr \cr
##' Finally, three initialization schemes for \code{Delta} are available with \code{settings$initStrategy}: 
##' \itemize{
##' \item for \code{settings$initStrategy == 'simple'}, \code{Delta} is simply initialized to the estimated \code{g} value of \code{modHom}. 
##' Note that this procedure may fail when \code{settings$penalty == TRUE}.
##' \item for \code{settings$initStrategy == 'residuals'}, \code{Delta} is initialized to the estimated residual variance from the homoskedastic mean prediction.
##' \item for \code{settings$initStrategy == 'smoothed'}, \code{Delta} takes the values predicted by \code{modNugs} at \code{X0}.
##' }
##'
##' Notice that \code{lower} and \code{upper} bounds cannot be equal for \code{\link[stats]{optim}}.
##'
##' @return a list which is given the S3 class \code{"hetTP"}, with elements:
##' \itemize{
##' \item \code{theta}: unless given, maximum likelihood estimate (mle) of the lengthscale parameter(s),
##' \item \code{Delta}: unless given, mle of the nugget vector (non-smoothed),
##' \item \code{Lambda}: predicted input noise variance at \code{X0}, 
##' \item \code{sigma2}: plugin estimator of the variance,
##' \item \code{theta_g}: unless given, mle of the lengthscale(s) of the noise/log-noise process,
##' \item \code{k_theta_g}: if \code{settings$linkThetas == 'joint'}, mle for the constant by which lengthscale parameters of \code{theta} are multiplied to get \code{theta_g},
##' \item \code{g}: unless given, mle of the nugget of the noise/log-noise process,
##' \item \code{trendtype}: either "\code{SK}" if \code{beta0} is provided, else "\code{OK}",
##' \item \code{beta0} constant trend of the mean process, plugin-estimator unless given,
##' \item \code{nmean}: plugin estimator for the constant noise/log-noise process mean,
## ' \item \code{pX}: if used, matrix of pseudo-inputs locations for the noise/log-noise process,
##' \item \code{ll}: log-likelihood value, (\code{ll_non_pen}) is the value without the penalty,
##' \item \code{nit_opt}, \code{msg}: \code{counts} and \code{message} returned by \code{\link[stats]{optim}}
##' \item \code{modHom}: homoskedastic GP model of class \code{homGP} used for initialization of the mean process,
##' \item \code{modNugs}: homoskedastic GP model of class \code{homGP} used for initialization of the noise/log-noise process,
##' \item \code{nu_hat_var}: variance of the noise process,
##' \item \code{used_args}: list with arguments provided in the call to the function, which is saved in \code{call},
##' \item \code{X0}, \code{Z0}, \code{Z}, \code{eps}, \code{logN}, \code{covtype}: values given in input,
##' \item \code{time}: time to train the model, in seconds.
##'}
##' @seealso \code{\link[hetGP]{predict.hetTP}} for predictions. 
##' \code{summary} and \code{plot} functions are available as well. 
##' @references 
##' M. Binois, Robert B. Gramacy, M. Ludkovski (2017+), Practical heteroskedastic Gaussian process modeling for large simulation experiments, arXiv preprint arXiv:1611.05902. \cr \cr
##' 
##' A. Shah, A. Wilson, Z. Ghahramani (2014), Student-t processes as alternatives to Gaussian processes, Artificial Intelligence and Statistics, 877--885.
##' 
##' @export
##' @importFrom stats optim var
##' @import methods
##' @examples 
##' ##------------------------------------------------------------
##' ## Example 1: Heteroskedastic TP modeling on the motorcycle data
##' ##------------------------------------------------------------
##' set.seed(32)
##' 
##' ## motorcycle data
##' library(MASS)
##' X <- matrix(mcycle$times, ncol = 1)
##' Z <- mcycle$accel
##' nvar <- 1
##' plot(X, Z, ylim = c(-160, 90), ylab = 'acceleration', xlab = "time")
##' 
##' 
##' ## Model fitting
##' model <- mleHetTP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar),
##'                   covtype = "Matern5_2")
##'             
##' ## Display averaged observations
##' points(model$X0, model$Z0, pch = 20)
##' 
##' ## A quick view of the fit                  
##' summary(model)
##' 
##' ## Create a prediction grid and obtain predictions
##' xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 
##' preds <- predict(x = xgrid, object =  model)
##' 
##' ## Display mean predictive surface
##' lines(xgrid, preds$mean, col = 'red', lwd = 2)
##' ## Display 95% confidence intervals
##' lines(xgrid, preds$mean + sqrt(preds$sd2) * qt(0.05, df = model$nu + nrow(X)), col = 2, lty = 2)
##' lines(xgrid, preds$mean + sqrt(preds$sd2) * qt(0.95, df = model$nu + nrow(X)), col = 2, lty = 2)
##' ## Display 95% prediction intervals
##' lines(xgrid, preds$mean + sqrt(preds$sd2 + preds$nugs) * qt(0.05, df = model$nu + nrow(X)),
##'   col = 3, lty = 2)
##' lines(xgrid, preds$mean + sqrt(preds$sd2 + preds$nugs) * qt(0.95, df = model$nu + nrow(X)), 
##'   col = 3, lty = 2)
##' 
##' ##------------------------------------------------------------
##' ## Example 2: 2D Heteroskedastic TP modeling
##' ##------------------------------------------------------------
##' set.seed(1)
##' nvar <- 2
##'   
##' ## Branin redefined in [0,1]^2
##' branin <- function(x){
##'   if(is.null(nrow(x)))
##'     x <- matrix(x, nrow = 1)
##'     x1 <- x[,1] * 15 - 5
##'     x2 <- x[,2] * 15
##'     (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
##' }
##' 
##' ## Noise field via standard deviation
##' noiseFun <- function(x){
##'   if(is.null(nrow(x)))
##'     x <- matrix(x, nrow = 1)
##'   return(1/5*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
##' }
##' 
##' ## data generating function combining mean and noise fields
##' ftest <- function(x){
##'   return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
##' }
##' 
##' ## Grid of predictive locations
##' ngrid <- 51
##' xgrid <- matrix(seq(0, 1, length.out = ngrid), ncol = 1) 
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' 
##' ## Unique (randomly chosen) design locations
##' n <- 100
##' Xu <- matrix(runif(n * 2), n)
##' 
##' ## Select replication sites randomly
##' X <- Xu[sample(1:n, 20*n, replace = TRUE),]
##' 
##' ## obtain training data response at design locations X
##' Z <- ftest(X)
##' 
##' ## Formating of data for model creation (find replicated observations) 
##' prdata <- find_reps(X, Z, rescale = FALSE, normalize = FALSE)
##'
##' ## Model fitting
##' model <- mleHetTP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z, ,
##'                   lower = rep(0.01, nvar), upper = rep(10, nvar),
##'                   covtype = "Matern5_2")
##'
##' ## a quick view into the data stored in the "hetTP"-class object
##' summary(model)                  
##'              
##' ## prediction from the fit on the grid     
##' preds <- predict(x = Xgrid, object =  model)
##' 
##' ## Visualization of the predictive surface
##' par(mfrow = c(2, 2))
##' contour(x = xgrid,  y = xgrid, z = matrix(branin(Xgrid), ngrid), 
##'   main = "Branin function", nlevels = 20)
##' points(X, col = 'blue', pch = 20)
##' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid), 
##'   main = "Predicted mean", nlevels = 20)
##' points(X, col = 'blue', pch = 20)
##' contour(x = xgrid,  y = xgrid, z = matrix(noiseFun(Xgrid), ngrid), 
##'   main = "Noise standard deviation function", nlevels = 20)
##' points(X, col = 'blue', pch = 20)
##' contour(x = xgrid,  y= xgrid, z = matrix(sqrt(preds$nugs), ngrid), 
##'   main = "Predicted noise values", nlevels = 20)
##' points(X, col = 'blue', pch = 20)
##' par(mfrow = c(1, 1))
##
mleHetTP <- function(X, Z, lower = NULL, upper = NULL,
                     noiseControl = list(k_theta_g_bounds = c(1, 100), g_max = 1e4, g_bounds = c(1e-6, 0.1),
                                         nu_bounds = c(2 + 1e-3, 30), sigma2_bounds = c(sqrt(.Machine$double.eps), 1e4)),
                     settings = list(linkThetas = 'joint', logN = TRUE, initStrategy = 'residuals', checkHom = TRUE,
                                     penalty = TRUE, trace = 0, return.matrices = TRUE, return.hom = FALSE, factr = 1e9), 
                     covtype = c("Gaussian", "Matern5_2", "Matern3_2"), maxit = 100, known = list(beta0 = 0),
                     init = list(nu = 3), eps = sqrt(.Machine$double.eps)){
  
  if(typeof(X)=="list"){
    X0 <- X$X0
    Z0 <- X$Z0
    mult <- X$mult
    if(sum(mult) != length(Z)) stop("Length(Z) should be equal to sum(mult)")
    if(is.null(dim(X0))) X0 <- matrix(X0, ncol = 1)
    if(length(Z0) != nrow(X0)) stop("Dimension mismatch between Z0 and X0")
  }else{
    if(is.null(dim(X))) X <- matrix(X, ncol = 1)
    if(nrow(X) != length(Z)) stop("Dimension mismatch between Z and X")
    elem <- find_reps(X, Z, return.Zlist = F)
    X0 <- elem$X0
    Z0 <- elem$Z0
    Z <- elem$Z
    mult <- elem$mult
  }
  
  covtype <- match.arg(covtype)
  
  if(is.null(lower) || is.null(upper)){
    auto_thetas <- auto_bounds(X = X0, covtype = covtype)
    if(is.null(lower)) lower <- auto_thetas$lower
    if(is.null(upper)) upper <- auto_thetas$upper
  }
  
  if(length(lower) != length(upper)) stop("upper and lower should have the same size")
  
  ## Save time to train model
  tic <- proc.time()[3]
  
  ## Initial checks
  
  n <- nrow(X0)
  
  if(is.null(n))
    stop("X0 should be a matrix. \n")

  jointThetas <- constrThetas <- FALSE
  
  if(is.null(settings$linkThetas)){
    jointThetas <- TRUE
  }else{
    if(settings$linkThetas == 'joint'){
      jointThetas <- TRUE
    }else{
      if(settings$linkThetas == 'constr')
        constrThetas <- TRUE
    }
  }
  
  logN <- TRUE
  if(!is.null(settings$logN)){
    logN <- settings$logN
  }
  
  if(is.null(settings$return.matrices))
    settings$return.matrices <- TRUE
  
  if(is.null(settings$return.hom))
    settings$return.hom <- FALSE
  
  if(jointThetas && is.null(noiseControl$k_theta_g_bounds))
    noiseControl$k_theta_g_bounds <- c(1, 100)
  
  if(is.null(settings$initStrategy))
    settings$initStrategy <- 'residuals'
  
  if(is.null(settings$factr))
    settings$factr <- 1e9
  
  penalty <- TRUE
  if(!is.null(settings$penalty))
    penalty <- settings$penalty
  
  # hardpenalty <- FALSE
  # if(!is.null(settings$hardpenalty))
  #   hardpenalty <- settings$hardpenalty
  
  if(is.null(settings$checkHom))
    settings$checkHom <- TRUE
  
  trace <- 0
  if(!is.null(settings$trace))
    trace <- settings$trace
  
  components <- NULL
  if(is.null(known$theta)){
    components <- c(components, list("theta"))
  }else{
    init$theta <- known$theta
  }
  
  if(is.null(known$Delta)){
    components <- c(components, list("Delta"))
  }else{
    init$Delta <- known$Delta
  }
  
  if(jointThetas && is.null(known$k_theta_g)){
    components <- c(components, list("k_theta_g"))
    
  }else{
    init$k_theta_g <- known$k_theta_g
  }
  
  if(!jointThetas && is.null(known$theta_g)){
    components <- c(components, list("theta_g"))
  }else{
    init$theta_g <- known$theta_g
  }
  
  if(is.null(known$sigma2)){
    components <- c(components, list("sigma2"))
    if(is.null(init$sigma2)) init$sigma2 <- var(Z0)
  }else{
    init$sigma2 <- known$sigma2
  }
  
  if(is.null(known$nu)){
    components <- c(components, list("nu"))
    if(is.null(init$nu)) init$nu <- 3
  }else{
    init$nu <- known$nu
  }
  
  if(is.null(known$g)){
    components <- c(components, list("g"))
  }else{
    init$g <- known$g
  }
  
  if(!is.null(init$pX)){
    components <- c(components, list("pX"))
    idcs_pX <- which(duplicated(rbind(init$pX, X0))) - nrow(init$pX) ## Indices of starting points in pX
  }else{
    if(!is.null(known$pX)){
      init$pX <- known$pX
    }else{
      init$pX <- NULL
    }
  }
  
  if(penalty && "pX" %in% components){
    penalty <- FALSE
    warning("Penalty not available with pseudo-inputs for now")
  }
  
  trendtype <- 'SK'
  if(is.null(known$beta0))  known$beta0 <- 0
  beta0 <- known$beta0
  
  if(is.null(components) && is.null(known$theta_g)){
    known$theta_g <- known$k_theta_g * known$theta
  }
  
  if(is.null(noiseControl$g_bounds)){
    noiseControl$g_bounds <- c(1e-6, 0.1)
  }
  
  if(is.null(noiseControl$nu_bounds)){
    noiseControl$nu_bounds <- c(2 + 1e-3, 30)
  }
  
  if(is.null(noiseControl$sigma2_bounds)){
    noiseControl$sigma2_bounds <- c(sqrt(.Machine$double.eps), 1e4)
  }
  
  # ## Advanced option Single Nugget Kriging model for the noise process
  # if(!is.null(noiseControl$SiNK) && noiseControl$SiNK){
  #   SiNK <- TRUE
  #   if(is.null(noiseControl$SiNK_eps)){
  #     SiNK_eps <- 1e-4
  #   }else{
  #     SiNK_eps <- noiseControl$SiNK_eps
  #   }
  # }else{
  #   SiNK <- FALSE
  #   SiNK_eps <- 1e-4
  # } 
  
  
  ### Automatic Initialisation
  modHom <- modNugs <- NULL
  if(is.null(init$theta) || is.null(init$Delta) || is.null(init$sigma2)){
    ## A) homoskedastic mean process
    
    if(!is.null(known$g_H)){
      g_init <- NULL
      # noiseControl$g_max <- known$g_H + eps
      # noiseControl$g_min <- known$g_H
    }else{
      g_init <- init$g_H
      ## Initial value for g of the homoskedastic process: based on the mean variance at replicates compared to the variance of Z0
      if(any(mult > 5)){
        mean_var_replicates <- mean((fast_tUY2(mult, (Z - rep(Z0, times = mult))^2)/mult)[which(mult > 5)])
        
        if(is.null(g_init)) g_init <- mean_var_replicates 
        if(is.null(noiseControl$g_max))
          noiseControl$g_max <- max(1e2, 100 * g_init)
        if(is.null(noiseControl$g_min))
          noiseControl$g_min <- eps
        
      }else{
        if(is.null(g_init)) g_init <- 0.05
        
        if(is.null(noiseControl$g_max))
          noiseControl$g_max <- 1e2
        
        if(is.null(noiseControl$g_min))
          noiseControl$g_min <- eps
        
      }
      
    }
    
    if(is.null(init$theta))
      init$theta <- sqrt(upper * lower)
    
    if(settings$checkHom){
      rKI <- TRUE
    }else{
      rKI <- FALSE
    } 
    modHom <- mleHomTP(X = list(X0 = X0, Z0 = Z0, mult = mult), Z = Z, lower = lower,
                       known = list(theta = known$theta, g = known$g_H, nu = known$nu, sigma2 = known$sigma2, beta0 = beta0),
                       init = list(theta = init$theta, g = g_init, nu = init$nu, sigma2 = init$sigma2),
                       upper = upper, covtype = covtype, maxit = maxit,
                       noiseControl = list(g_bounds = c(noiseControl$g_min, noiseControl$g_max),
                                           nu_bounds = noiseControl$nu_bounds), eps = eps, 
                       settings = list(return.Ki = rKI, factr = settings$factr, pgtol = settings$pgtol))
    
    ## Update init values
    if(is.null(known$theta)) init$theta <- modHom$theta
    
    if(is.null(known$sigma2)) init$sigma2 <- modHom$sigma2
    
    if(is.null(known$nu) && is.null(init$nu)) init$nu <- modHom$nu # Estimating nu may be hard, and even more for homTP
    
    if(is.null(init$Delta)){
      predHom <- suppressWarnings(predict(x = X0, object = modHom)$mean)
      nugs_est <- (Z - rep(predHom, times = mult))^2 #squared deviation from the homoskedastic prediction mean to the actual observations

      if(logN){
        nugs_est <- log(nugs_est)
      }
      
      nugs_est0 <- drop(fast_tUY2(mult, nugs_est))/mult # average
      
    }else{
      nugs_est0 <- init$Delta
    }
    
    if(constrThetas){
      noiseControl$lowerTheta_g <- modHom$theta
    }
    
    if(settings$initStrategy == 'simple'){
      if(logN){
        init$Delta <- rep(log(modHom$g), nrow(X0))
      }else{
        init$Delta <- rep(modHom$g, nrow(X0))
      }
    }
    if(settings$initStrategy == 'residuals')
      init$Delta <- nugs_est0
    
  }
  
  
  if((is.null(init$theta_g) && is.null(init$k_theta_g)) || is.null(init$g)){
    
    ## B) Homegeneous noise process
    if(jointThetas){
      if(is.null(init$k_theta_g)){
        init$k_theta_g <- 1
        init$theta_g <- init$theta
      }else{
        init$theta_g <- init$k_theta_g * init$theta
      }
      
      if(is.null(noiseControl$lowerTheta_g)){
        noiseControl$lowerTheta_g <- init$theta_g - eps
      }
      
      if(is.null(noiseControl$upperTheta_g)){
        noiseControl$upperTheta_g <- init$theta_g + eps
      }
    }
    
    if(!jointThetas && is.null(init$theta_g)){
      init$theta_g <- init$theta
    }
    
    if(is.null(noiseControl$lowerTheta_g)){
      noiseControl$lowerTheta_g <- lower
    }
    
    if(is.null(noiseControl$upperTheta_g)){
      noiseControl$upperTheta_g <- upper
    }
    
    ## If an homegeneous process of the mean has already been computed, it is used for estimating the parameters of the noise process
    if(exists("nugs_est")){
      
      if(is.null(init$g)){
        mean_var_replicates_nugs <- mean((fast_tUY2(mult, (nugs_est - rep(nugs_est0, times = mult))^2)/mult))
        init$g <- mean_var_replicates_nugs
      }
      
      modNugs <- mleHomGP(X = list(X0 = X0, Z0 = nugs_est0, mult = mult), Z = nugs_est,
                          lower = noiseControl$lowerTheta_g, upper = noiseControl$upperTheta_g,
                          init = list(theta = init$theta_g, g =  init$g), covtype = covtype, noiseControl = noiseControl,
                          maxit = maxit, eps = eps, settings = list(return.Ki = F, factr = settings$factr, pgtol = settings$pgtol))
      prednugs <- suppressWarnings(predict(x = X0, object = modNugs))
      
    }else{
      if(!exists("nugs_est0")) nugs_est0 <- init$Delta
      
      if(is.null(init$g)){
        init$g <- 0.05 * var(nugs_est0)
      }
      
      modNugs <- mleHomGP(X = list(X0 = X0, Z0 = nugs_est0, mult = rep(1, nrow(X0))),
                          Z = nugs_est0,
                          lower = noiseControl$lowerTheta_g, upper = noiseControl$upperTheta_g,
                          init = list(theta = init$theta_g, g =  init$g), covtype = covtype, noiseControl = noiseControl,
                          maxit = maxit, eps = eps, settings = list(return.Ki = F, factr = settings$factr, pgtol = settings$pgtol))
      prednugs <- suppressWarnings(predict(x = X0, object = modNugs))
    }
    
    if(settings$initStrategy == 'simple'){
      if(logN){
        init$Delta <- rep(log(modHom$g), nrow(X0))
      }else{
        init$Delta <- rep(modHom$g, nrow(X0))
      }
    }
    if(settings$initStrategy == 'residuals')
      init$Delta <- nugs_est0
    if(settings$initStrategy == 'smoothed')
      init$Delta <- prednugs$mean
    
    if(is.null(known$g))
      init$g <- modNugs$g
    
    if(jointThetas && is.null(init$k_theta_g))
      init$k_theta_g <- 1
    
    if(!jointThetas && is.null(init$theta_g))
      init$theta_g <- modNugs$theta
    
  }
  
  if(is.null(noiseControl$lowerTheta_g)){
    noiseControl$lowerTheta_g <- lower
  }
  
  if(is.null(noiseControl$upperTheta_g)){
    noiseControl$upperTheta_g <- upper
  }
  
  ### Start of optimization of the log-likelihood
  fn <- function(par, X0, Z0, Z, mult, Delta = NULL, theta = NULL, g = NULL, k_theta_g = NULL, theta_g = NULL,
                 nu = NULL, sigma2 = NULL, logN = FALSE, beta0 = NULL, hom_ll, env = NULL){
    
    idx <- 1 # to store the first non used element of par
    
    if(is.null(theta)){
      idx <- length(init$theta)
      theta <- par[1:idx]
      idx <- idx + 1
    }
    if(is.null(Delta)){
      Delta <- par[idx:(idx - 1 + length(init$Delta))]
      idx <- idx + length(init$Delta)
    }
    if(jointThetas && is.null(k_theta_g)){
      k_theta_g <- par[idx]
      idx <- idx + 1
    }
    if(!jointThetas && is.null(theta_g)){
      theta_g <- par[idx:(idx - 1 + length(init$theta_g))]
      idx <- idx + length(init$theta_g)
    }
    
    if(is.null(sigma2)){
      sigma2 <- par[idx]
      idx <- idx + 1
    }
    
    if(is.null(nu)){
      nu <- par[idx]
      idx <- idx + 1
    }
    
    if(is.null(g)){
      g <- par[idx]
    }
    
    loglik <-logLik_HetTP(X0 = X0, Z0 = Z0, Z = Z, mult = mult, Delta = Delta, theta = theta, g = g, k_theta_g = k_theta_g, theta_g = theta_g,
                          nu = nu, sigma2 = sigma2, logN = logN, beta0 = beta0, covtype = covtype, eps = eps,
                          penalty = penalty, hom_ll = hom_ll, env = env, trace = trace)
    
    if(!is.null(env) && !is.na(loglik)){
      if(is.null(env$max_loglik) || loglik > env$max_loglik){
        env$max_loglik <- loglik
        env$arg_max <- par
      }
    }
    
    return(loglik)
  }
  
  
  gr <- function(par, X0, Z0, Z, mult, Delta = NULL, theta = NULL, g = NULL, k_theta_g = NULL, theta_g = NULL, nu = NULL, sigma2 = NULL, logN = FALSE,
                 beta0 = NULL, pX = NULL, hom_ll, env = NULL){
    
    idx <- 1 # to store the first non used element of par
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
    }
    if(is.null(Delta)){
      Delta <- par[idx:(idx - 1 + length(init$Delta))]
      idx <- idx + length(init$Delta)
    }
    if(jointThetas && is.null(k_theta_g)){
      k_theta_g <- par[idx]
      idx <- idx + 1
    }
    if(!jointThetas && is.null(theta_g)){
      theta_g <- par[idx:(idx - 1 + length(init$theta_g))]
      idx <- idx + length(init$theta_g)
    }
    if(is.null(sigma2)){
      sigma2 <- par[idx]
      idx <- idx + 1
    }
    if(is.null(nu)){
      nu <- par[idx]
      idx <- idx + 1
    }
    if(is.null(g)){
      g <- par[idx]
    }
    
    # c(grad(logLikHet_Wood_Ani_snugs_pX, X0 = X0, Z0 = Z0, Z = Z, mult = mult, Delta = Delta, x = theta, g = g, k_theta_g = k_theta_g, theta_g = theta_g,
    #        logN = logN, SiNK = SiNK, beta0 = beta0, pX = pX),
    #   grad(logLikHet_Wood_Ani_snugs_pX, X0 = X0, Z0 = Z0, Z = Z, mult = mult, x = Delta, theta = theta, g = g, k_theta_g = k_theta_g, theta_g = theta_g,
    #        logN = logN, SiNK = SiNK, beta0 = beta0, pX = pX),
    #   grad(logLikHet_Wood_Ani_snugs_pX, X0 = X0, Z0 = Z0, Z = Z, mult = mult, Delta = Delta, theta = theta, g = g, x = k_theta_g, theta_g = theta_g,
    #        logN = logN, SiNK = SiNK, beta0 = beta0, pX = pX),
    #   grad(logLikHet_Wood_Ani_snugs_pX, X0 = X0, Z0 = Z0, Z = Z, mult = mult, Delta = Delta, theta = theta, x = g, k_theta_g = k_theta_g, theta_g = theta_g,
    #        logN = logN, SiNK = SiNK, beta0 = beta0, pX = pX)
    # )
    # grad(fn, x = par, X0 = X0, Z0 = Z0, Z = Z, mult = mult, logN = logN, beta0 = beta0, method.args=list(r = 6))
    
    return(dlogLikHetTP(X0 = X0, Z0 = Z0, Z = Z, mult = mult, Delta = Delta, theta = theta, g = g, k_theta_g = k_theta_g, theta_g = theta_g,
                        nu = nu, sigma2 = sigma2, logN = logN, beta0 = beta0, components = components, covtype = covtype, eps = eps,
                        penalty = penalty, hom_ll = hom_ll, env = env))
    
  }
  
  ## Pre-processing for heterogeneous fit
  parinit <- lowerOpt <- upperOpt <- NULL
  
  if(trace > 2)
    cat("Initial value of the parameters:\n")
  
  if(is.null(known$theta)){
    parinit <- init$theta
    lowerOpt <- lower
    upperOpt <- upper
    if(trace > 2) cat("Theta: ", init$theta, "\n")
  }
  
  if(is.null(known$Delta)){
    
    if(is.null(noiseControl$lowerDelta)){
      if(logN){
        noiseControl$lowerDelta <- log(eps)
      }else{
        noiseControl$lowerDelta <- eps
      }
    }
    
    if(length(noiseControl$lowerDelta) == 1){
      noiseControl$lowerDelta <- rep(noiseControl$lowerDelta, n)
    }
    
    if(is.null(noiseControl$g_max)) noiseControl$g_max <- 1e2 
    
    if(is.null(noiseControl$upperDelta)){
      if(logN){
        noiseControl$upperDelta <- log(noiseControl$g_max)
      }else{
        noiseControl$upperDelta <- noiseControl$g_max
      }
    }
    
    if(length(noiseControl$upperDelta) == 1){
      noiseControl$upperDelta <- rep(noiseControl$upperDelta, n)
    }
    
    lowerOpt <- c(lowerOpt, noiseControl$lowerDelta)
    upperOpt <- c(upperOpt, noiseControl$upperDelta)
    parinit <- c(parinit, init$Delta)
    
    if(trace > 2) cat("Delta: ", init$Delta, "\n")
  }
  
  if(jointThetas && is.null(known$k_theta_g)){
    parinit <- c(parinit, init$k_theta_g)
    lowerOpt <- c(lowerOpt, noiseControl$k_theta_g_bounds[1])
    upperOpt <- c(upperOpt, noiseControl$k_theta_g_bounds[2])
    if(trace > 2) cat("k_theta_g: ", init$k_theta_g, "\n")
  }
  
  if(!jointThetas && is.null(known$theta_g)){
    parinit <- c(parinit, init$theta_g)
    lowerOpt <- c(lowerOpt, noiseControl$lowerTheta_g)
    upperOpt <- c(upperOpt, noiseControl$upperTheta_g)
    if(trace > 2) cat("theta_g: ", init$theta_g, "\n")
  }
  
  if(is.null(known$sigma2)){
    parinit <- c(parinit, init$sigma2)
    lowerOpt <- c(lowerOpt, noiseControl$sigma2_bounds[1])
    upperOpt <- c(upperOpt, noiseControl$sigma2_bounds[2])
  }
  
  if(is.null(known$nu)){
    parinit <- c(parinit, init$nu)
    lowerOpt <- c(lowerOpt, noiseControl$nu_bounds[1])
    upperOpt <- c(upperOpt, noiseControl$nu_bounds[2])
  }
  
  if(is.null(known$g)){
    parinit <- c(parinit, init$g)
    if(trace > 2) cat("g: ", init$g, "\n")
    
    lowerOpt <- c(lowerOpt, noiseControl$g_bounds[1])
    upperOpt <- c(upperOpt, noiseControl$g_bounds[2])
  }
  
  if("pX" %in% components){
    parinit <- c(parinit, as.vector(init$pX))
    lowerOpt <- c(lowerOpt, as.vector(rep(noiseControl$lowerpX, times = nrow(init$pX))))
    upperOpt <- c(upperOpt, as.vector(rep(noiseControl$upperpX, times = nrow(init$pX))))
    if(trace > 2) cat("pX: ", as.vector(init$pX), "\n")
  }
  
  
  ## Case when some parameters need to be estimated
  mle_par <- known # Store infered and known parameters
  if(!is.null(components)){
    
    if(!is.null(modHom)){
      hom_ll <- modHom$ll
    }else{
      ## Compute reference homoskedastic likelihood, with fixed theta for speed
      modHom_tmp <- mleHomTP(X = list(X0 = X0, Z0 = Z0, mult = mult), Z = Z, lower = lower, upper = upper, upper,
                             known = list(theta = known[["theta"]], g = known$g_H, nu = known$nu, sigma2 = known$sigma2, beta0 = 0),
                             covtype = covtype, init = init,
                             noiseControl = noiseControl, eps = eps, settings = list(return.Ki = F))
      
      hom_ll <- modHom_tmp$ll
    } 
    
    ## Maximization of the log-likelihood
    envtmp <- environment()
    out <- try(optim(par = parinit, fn = fn, gr = gr, method = "L-BFGS-B", lower = lowerOpt, upper = upperOpt,
                 X0 = X0, Z0 = Z0, Z = Z, mult = mult, logN = logN, Delta = known$Delta, theta = known$theta,
                 g = known$g, k_theta_g = known$k_theta_g, theta_g = known$theta_g, hom_ll = hom_ll, env = envtmp,
                 nu = known$nu, beta0 = known$beta0, sigma2 = known$sigma2, 
                 control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol)))
    
    ## Catch errors when at least one likelihood evaluation worked
    if(class(out) == "try-error")
      out <- list(par = envtmp$arg_max, value = envtmp$max_loglik, counts = NA,
                  message = "Optimization stopped due to NAs, use best value so far")
    
    ## Temporary
    if(trace > 0){
      print(out$message)
      cat("Number of variables at boundary values:" , length(which(out$par == upperOpt)) + length(which(out$par == lowerOpt)), "\n")
    }
    
    if(trace > 1)
      cat("Name | Value | Lower bound | Upper bound \n")
    
    ## Post-processing
    idx <- 1
    if(is.null(known[["theta"]])){
      mle_par$theta <- out$par[1:length(init$theta)]
      idx <- idx + length(init$theta)
      if(trace > 1) cat("Theta |", mle_par$theta, " | ", lower, " | ", upper, "\n")
    }

    if(is.null(known$Delta)){
      mle_par$Delta <- out$par[idx:(idx - 1 + length(init$Delta))]
      idx <- idx + length(init$Delta)
      if(trace > 1){
        for(ii in 1:ceiling(length(mle_par$Delta)/5)){
          i_tmp <- (1 + 5*(ii-1)):min(5*ii, length(mle_par$Delta))
          if(logN) cat("Delta |", mle_par$Delta[i_tmp], " | ", pmax(log(eps * mult[i_tmp]), init$Delta[i_tmp] - log(1000)), " | ", init$Delta[i_tmp] + log(100), "\n")
          if(!logN) cat("Delta |", mle_par$Delta[i_tmp], " | ", pmax(mult[i_tmp] * eps, init$Delta[i_tmp] / 1000), " | ", init$Delta[i_tmp] * 100, "\n")
        }
      }
    }
    
    if(jointThetas){
      if(is.null(known$k_theta_g)){
        mle_par$k_theta_g <- out$par[idx]
        idx <- idx + 1
      } 
      mle_par$theta_g <- mle_par$k_theta_g * mle_par$theta
      
      if(trace > 1) cat("k_theta_g |", mle_par$k_theta_g, " | ", noiseControl$k_theta_g_bounds[1], " | ", noiseControl$k_theta_g_bounds[2], "\n")
    }
    
    if(!jointThetas && is.null(known$theta_g)){
      mle_par$theta_g <- out$par[idx:(idx - 1 + length(init$theta_g))]
      idx <- idx + length(init$theta_g)
      if(trace > 1) cat("theta_g |", mle_par$theta_g, " | ", noiseControl$lowerTheta_g, " | ", noiseControl$upperTheta_g, "\n")
    }
    
    if(is.null(known$sigma2)){
      mle_par$sigma2 <- out$par[idx]
      idx <- idx + 1
    }
    
    if(is.null(known$nu)){
      mle_par$nu <- out$par[idx]
      idx <- idx + 1
    }
    
    if(is.null(known$g)){
      mle_par$g <- out$par[idx]
      idx <- idx + 1
      if(trace > 1) cat("g |", mle_par$g, " | ", noiseControl$g_bounds[1], " | ", noiseControl$g_bounds[2], "\n")
    }
    # if(idx != (length(out$par) + 1)){
    #   mle_par$pX <- matrix(out$par[idx:length(out$par)], ncol = ncol(X0))
    #   if(trace > 1) cat("pX |", as.vector(mle_par$pX), " | ", as.vector(rep(noiseControl$lowerpX, times = nrow(init$pX))), " | ", as.vector(rep(noiseControl$upperpX, times = nrow(init$pX))), "\n")
    # }
    
  }else{
    out <- list(message = "All hyperparameters given, no optimization \n", count = 0, value = NULL)
  }
  
  ## Computation of nu2
  
  if(penalty){
    ll_non_pen <- logLik_HetTP(X0 = X0, Z0 = Z0, Z = Z, mult = mult, Delta = mle_par$Delta, theta = mle_par$theta, g = mle_par$g,
                               k_theta_g = mle_par$k_theta_g, theta_g = mle_par$theta_g, sigma2 = mle_par$sigma2, nu = mle_par$nu,
                               logN = logN, beta0 = mle_par$beta0, covtype = covtype, eps = eps, penalty = FALSE, trace = trace)
  }else{
    ll_non_pen <- out$value
  }
  # if(!is.null(modHom)){
  #   if(modHom$ll >= ll_non_pen){
  #     if(trace >= 0) cat("homoskedastic model has higher log-likelihood: ", modHom$ll, " compared to ", ll_non_pen, "\n")
  #     if(settings$checkHom){
  #       if(trace >= 0) cat("Return homoskedastic model \n")
  #       return(modHom)
  #     }
  #   }
  # }
  
  # if(is.null(mle_par$pX)){
  if(jointThetas){
    Cg <- cov_gen(X1 = X0, theta = mle_par$k_theta_g * mle_par$theta, type = covtype)
  }else{
    Cg <- cov_gen(X1 = X0, theta = mle_par$theta_g, type = covtype)
  } 
  Kgi <- chol2inv(chol(Cg + diag(eps + mle_par$g/mult)))
  
  nmean <- drop(rowSums(Kgi) %*% mle_par$Delta / sum(Kgi)) ## ordinary kriging mean
  
  nu_hat_var <- max(eps, drop(crossprod(mle_par$Delta - nmean, Kgi) %*% (mle_par$Delta - nmean))/length(mle_par$Delta))
  
  # if(SiNK){
  #   rhox <- 1 / rho_AN(xx = X0, X0 = mle_par$pX, theta_g = mle_par$theta_g, g = mle_par$g, type = covtype, eps = eps, SiNK_eps = SiNK_eps, mult = mult)
  #   M <-  rhox * Cg %*% (Kgi %*% (mle_par$Delta - nmean))
  # }else{
  M <- Cg %*% (Kgi %*% (mle_par$Delta - nmean))
  # }
  # }else{
  #   Cg <- cov_gen(X1 = mle_par$pX, theta = mle_par$theta_g, type = covtype)
  #   Kgi <- chol2inv(chol(Cg + diag(eps + mle_par$g/mult)))
  #   
  #   kg <- cov_gen(X1 = X0, X2 = mle_par$pX, theta = mle_par$theta_g, type = covtype)
  #   
  #   nmean <- drop(rowSums(Kgi) %*% mle_par$Delta / sum(Kgi)) ## ordinary kriging mean
  #   
  #   nu_hat_var <- max(eps, drop(crossprod(mle_par$Delta - nmean, Kgi) %*% (mle_par$Delta - nmean))/length(mle_par$Delta))
  #   
  #   if(SiNK){
  #     rhox <- 1 / rho_AN(xx = X0, X0 = mle_par$pX, theta_g = mle_par$theta_g, g = mle_par$g, type = covtype, eps = eps, SiNK_eps = SiNK_eps, mult = mult)
  #     M <-  rhox * kg %*% (Kgi %*% (mle_par$Delta - nmean))
  #   }else{
  #     M <- kg %*% (Kgi %*% (mle_par$Delta - nmean))
  #   }
  # }
  
  
  Lambda <- drop(nmean + M)
  
  if(logN){
    Lambda <- exp(Lambda)
  }
  else{
    Lambda[Lambda <= 0] <- eps
  }
  # Lambda <- pmax(mult * eps, Lambda)
  
  LambdaN <- rep(Lambda, times = mult)
  
  Ki <- chol2inv(chol(add_diag(mle_par$sigma2 * cov_gen(X1 = X0, theta = mle_par$theta, type = covtype), Lambda/mult + eps))) 
  ldetKi <- determinant(Ki, logarithm = TRUE)[[1]][1]
  
  if(is.null(known$beta0))
    mle_par$beta0 <- drop(colSums(Ki) %*% Z0 / sum(Ki))
  
  psi_0 <- drop(crossprod(Z0 - mle_par$beta0, Ki) %*% (Z0 - mle_par$beta0))
  
  psi <- (crossprod((Z - mle_par$beta0)/LambdaN, Z - mle_par$beta0) - crossprod((Z0 - mle_par$beta0) * mult/Lambda, Z0 - mle_par$beta0) + psi_0)
  
  res <- list(theta = mle_par$theta, Delta = mle_par$Delta, psi = as.numeric(psi), beta0 = mle_par$beta0,
              k_theta_g = mle_par$k_theta_g, theta_g = mle_par$theta_g, g = mle_par$g, nmean = nmean, Lambda = Lambda,
              ll = out$value, ll_non_pen = ll_non_pen, nit_opt = out$counts, logN = logN, covtype = covtype, pX = mle_par$pX, msg = out$message,
              X0 = X0, Z0 = Z0, Z = Z, mult = mult, trendtype = trendtype, eps = eps,
              nu_hat_var = nu_hat_var, call = match.call(), nu = mle_par$nu, sigma2 = mle_par$sigma2,
              used_args = list(noiseControl = noiseControl, settings = settings, lower = lower, upper = upper, known = known),
              time = proc.time()[3] - tic)
  
  if(settings$return.matrices){
    res <- c(res, list(Ki = Ki, Kgi = Kgi))
  }
  
  if(settings$return.hom){
    res <- c(res, list(modHom = modHom, modNugs = modNugs))
  }
  
  class(res) <- "hetTP"
  
  return(res)
  
}

##'Student-t process predictions using a heterogeneous noise TP object (of class \code{hetTP}) 
##' @param x matrix of designs locations to predict at
##' @param object an object of class \code{hetTP}; e.g., as returned by \code{\link[hetGP]{mleHetTP}}
##' @param noise.var should the variance of the latent variance process be returned?
##' @param xprime optional second matrix of predictive locations to obtain the predictive covariance matrix between \code{x} and \code{xprime}
##' @param nugs.only if TRUE, only return noise variance prediction
##' @param ... no other argument for this method.
##' @return list with elements
##' \itemize{
##' \item \code{mean}: kriging mean;
##' \item \code{sd2}: kriging variance (filtered, e.g. without the nugget values)
##' \item \code{nugs}: noise variance
##' \item \code{sd2_var}: (optional) kriging variance of the noise process (i.e., on log-variances if \code{logN = TRUE})
##' \item \code{cov}: (optional) predictive covariance matrix between x and xprime
##' }
##' @details The full predictive variance corresponds to the sum of \code{sd2} and \code{nugs}.
##' @method predict hetTP 
##' @export
predict.hetTP <- function(object, x, noise.var = FALSE, xprime = NULL, nugs.only = FALSE, ...){
  
  if(is.null(dim(x))){
    x <- matrix(x, nrow = 1)
    if(ncol(x) != ncol(object$X0)) stop("x is not a matrix")
  }
  
  if(!is.null(xprime) && is.null(dim(xprime))){
    xprime <- matrix(xprime, nrow = 1)
    if(ncol(xprime) != ncol(object$X0)) stop("xprime is not a matrix")
  }
  
  n1 <- length(object$Z)
  if(is.null(object$Kgi)){
    if(is.null(object$pX)){
      Cg <- cov_gen(X1 = object$X0, theta = object$theta_g, type = object$covtype)
    }else{
      Cg <- cov_gen(X1 = object$pX, theta = object$theta_g, type = object$covtype)
    }
    object$Kgi <- chol2inv(chol(add_diag(Cg, object$eps + object$g/object$mult)))
  }
  
  # if(is.null(object$pX)){
  kg <- cov_gen(X1 = x, X2 = object$X0, theta = object$theta_g, type = object$covtype)
  # }else{
  #   kg <- cov_gen(X1 = x, X2 = object$pX, theta = object$theta_g, type = object$covtype)
  # }
  
  if(is.null(object$Ki))
    object$Ki <- chol2inv(chol(add_diag(object$sigma2 * cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype), object$Lambda/object$mult + object$eps)))
  
  # if(object$SiNK){
  #   M <-  1/rho_AN(xx = x, X0 = object$pX, theta_g = object$theta_g, g = object$g,
  #                  type = object$covtype, SiNK_eps = object$SiNK_eps, eps = object$eps, mult = object$mult) * kg %*% (object$Kgi %*% (object$Delta - object$nmean))
  # }else{
  M <- kg %*% (object$Kgi %*% (object$Delta - object$nmean))
  # }
  
  if(object$logN){
    nugs <- exp(drop(object$nmean + M))
  }else{
    nugs <- pmax(0, drop(object$nmean + M))
  }
  
  if(nugs.only){
    return(list(nugs = nugs))
  }
  
  # if(noise.var){
  #   if(is.null(object$nu_hat_var))
  #     object$nu_hat_var <- max(object$eps, drop(crossprod(object$Delta - object$nmean, object$Kgi) %*% (object$Delta - object$nmean))/length(object$Delta)) ## To avoid 0 variance
  #   sd2var = object$nu_hat * object$nu_hat_var* drop(1 - fast_diag(kg, tcrossprod(object$Kgi, kg)) + (1 - tcrossprod(rowSums(object$Kgi), kg))^2/sum(object$Kgi))
  # }else{
  #   sd2var = NULL
  # }
  
  kx <- object$sigma2 * cov_gen(X1 = x, X2 = object$X0, theta = object$theta, type = object$covtype)
  
  # if(object$trendtype == 'SK'){
  sd2 <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * drop(object$sigma2 - fast_diag(kx, tcrossprod(object$Ki, kx)))
  # }else{
  #   sd2 <- drop(object$sigma2 *  - fast_diag(kx, tcrossprod(object$Ki, kx)) + (1 - tcrossprod(rowSums(object$Ki), kx))^2/sum(object$Ki))
  # }
  
  ## In case of numerical errors, some sd2 values may become negative
  if(any(sd2 < 0)){
    # object$Ki <- ginv(add_diag(object$sigma2 * cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype), object$Lambda/object$mult + object$eps))
    # sd2 <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * drop(object$sigma2 - fast_diag(kx, tcrossprod(object$Ki, kx)))
    
    sd2 <- pmax(0, sd2)
    warning("Numerical errors caused some negative predictive variances to be thresholded to zero. Consider using ginv via rebuild.hetTP")
  }
  
  
  if(!is.null(xprime)){
    kxprime <- object$sigma2 * cov_gen(X1 = object$X0, X2 = xprime, theta = object$theta, type = object$covtype)
    # if(object$trendtype == 'SK'){
    if(nrow(x) > nrow(xprime)){
      cov <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * (object$sigma2 * cov_gen(X1 = x, X2 = xprime, theta = object$theta, type = object$covtype) - kx %*% object$Ki %*% kxprime)
    }else{
      cov <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * (object$sigma2 * cov_gen(X1 = x, X2 = xprime, theta = object$theta, type = object$covtype) - (kx %*% object$Ki) %*% kxprime)
    }
    
    # }else{
    #   cov <- (object$nu + object$psi - 2) / (object$nu + n1 - 2) * (object$sigma2 * cov_gen(X1 = x, X2 = xprime, theta = object$theta, type = object$covtype) - kx %*% tcrossprod(object$Ki, kxprime) + crossprod(1 - tcrossprod(rowSums(object$Ki), kx), 1 - tcrossprod(rowSums(object$Ki), kxprime))/sum(object$Ki))
    # }
  }else{
    cov = NULL
  }
  
  return(list(mean = drop(object$beta0 + kx %*% (object$Ki %*% (object$Z0 - object$beta0))),
              sd2 = sd2,
              nugs = nugs,
              # sd2var = sd2var,
              cov = cov))
  
}


##' @method summary hetTP
##' @export
summary.hetTP <- function(object,...){
  ans <- object
  class(ans) <- "summary.hetTP"
  ans
}

##' @export
print.summary.hetTP <- function(x, ...){
  cat("N = ", length(x$Z), " n = ", length(x$Z0), " d = ", ncol(x$X0), "\n")
  cat(x$covtype, " covariance lengthscale values of the main process: ", x$theta, "\n")
  cat("Variance/scale hyperparameter: ", x$sigma2, "\n")
  cat("Degree of freedom nu: ", x$nu, "\n")
  
  cat("Summary of Lambda values: \n")
  print(summary(x$Lambda))
  
  cat("Given constant trend value: ", x$beta0, "\n")
  
  
  if(x$logN){
    cat(x$covtype, " covariance lengthscale values of the log-noise process: ", x$theta_g, "\n")
    cat("Nugget of the log-noise process: ", x$g, "\n")
    cat("Estimated constant trend value of the log-noise process: ", x$nmean, "\n")
  }else{
    cat(x$covtype, " covariance lengthscale values of the log-noise process: ", x$theta_g, "\n")
    cat("Nugget of the noise process: ", x$g, "\n")
    cat("Estimated constant trend value of the noise process: ", x$nmean, "\n")
  }
  
  
  cat("MLE optimization: \n", "Log-likelihood = ", x$ll, "; Nb of evaluations (obj, gradient) by L-BFGS-B: ", x$nit_opt, "; message: ", x$msg, "\n")
  
}

##' @method print hetTP
##' @export
print.hetTP <- function(x, ...){
  print(summary(x))
}

##' @method plot hetTP
##' @export
plot.hetTP <- function(x, ...){
  LOOpreds <- LOO_preds(x)
  plot(x$Z, LOOpreds$mean[rep(1:nrow(x$X0), times = x$mult)], xlab = "Observed values", ylab = "Predicted values",
       main = "Leave-one-out predictions")
  arrows(x0 = LOOpreds$mean + sqrt(LOOpreds$sd2) * qt(0.05, df = x$nu + nrow(x$X0)),
         x1 = LOOpreds$mean + sqrt(LOOpreds$sd2) * qt(0.95, df = x$nu + nrow(x$X0)),
         y0 = LOOpreds$mean, length = 0, col = "blue")
  points(x$Z0[which(x$mult > 1)], LOOpreds$mean[which(x$mult > 1)], pch = 20, col = 2)
  abline(a = 0, b = 1, lty = 3)
  legend("topleft", pch = c(1, 20, NA), lty = c(NA, NA, 1), col = c(1, 2, 4),
         legend = c("observations", "averages (if > 1 observation)", "LOO prediction interval"))
}

## ' Rebuild inverse covariance matrices of \code{hetTP} (e.g., if exported without inverse matrices \code{Kgi} and/or \code{Ki})
## ' @param object \code{hetGP} model without slots \code{Ki} and/or \code{Kgi} (inverse covariance matrices)
## ' @param robust if \code{TRUE} \code{\link[MASS]{ginv}} is used for matrix inversion, otherwise it is done via Cholesky.
##' @rdname ExpImp
##' @method rebuild hetTP
##' @export
rebuild.hetTP <- function(object, robust = FALSE){
  if(robust){
    object$Kgi <- ginv(add_diag(cov_gen(X1 = object$X0, theta = object$theta_g, type = object$covtype), object$eps + object$g/object$mult))
    
    object$Ki <- ginv(add_diag(object$sigma2 * cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype), object$Lambda/object$mult + object$eps))
  }else{
    object$Kgi <- chol2inv(chol(add_diag(cov_gen(X1 = object$X0, theta = object$theta_g, type = object$covtype), object$eps + object$g/object$mult)))
    
    object$Ki <- chol2inv(chol(add_diag(object$sigma2 * cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype), object$Lambda/object$mult + object$eps)))
  }
  return(object)
}


