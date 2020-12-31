################################################################################
### Contains update of homGP/hetGP models
################################################################################

if(!isGeneric("update")) {
  setGeneric(name = "update",
             def = function(object, ...) standardGeneric("update")
  )
}

##' Fast update of existing \code{hetGP} model with new observations. 
##' @title Update \code{"hetGP"}-class model fit with new observations
##' @param object previously fit \code{"hetGP"}-class model
##' @param Xnew matrix of new design locations; \code{ncol(Xnew)} must match the input dimension encoded in \code{object}
##' @param Znew vector new observations at those design locations, of length \code{nrow(X)}. \code{NA}s can be passed, see Details
##' @param ginit minimal value of the smoothing parameter (i.e., nugget of the noise process) for optimization initialisation.
##' It is compared to the \code{g} hyperparameter in the object.   
##' @param lower,upper,noiseControl,settings,known optional bounds for mle optimization, see \code{\link[hetGP]{mleHetGP}}. 
##' If not provided, they are extracted from the existing model 
##' @param maxit maximum number of iterations for the internal L-BFGS-B optimization method; see \code{\link{optim}} for more details
##' @param method one of \code{"quick"}, \code{"mixed"} see Details.
##' @param ... no other argument for this method.
##' @details
##' 
##' The update can be performed with or without re-estimating hyperparameter.
##' In the first case, \code{\link[hetGP]{mleHetGP}} is called, based on previous values for initialization. 
##' The only missing values are the latent variables at the new points, that are initialized based on two possible update schemes in \code{method}:
##' \itemize{
##'   \item \code{"quick"} the new delta value is the predicted nugs value from the previous noise model;
##'   \item \code{"mixed"} new values are taken as the barycenter between prediction given by the noise process and empirical variance. 
##' }
##' The subsequent number of MLE computations can be controlled with \code{maxit}.
##' 
##' In case hyperparameters need not be updated, \code{maxit} can be set to \code{0}. 
##' In this case it is possible to pass \code{NA}s in \code{Znew}, then the model can still be used to provide updated variance predictions.
##' 
##' @export
##' @method update hetGP
##' @importFrom stats rnorm
##' @examples 
##' ##------------------------------------------------------------
##' ## Sequential update example
##' ##------------------------------------------------------------
##' set.seed(42)
##'
##' ## Spatially varying noise function
##' noisefun <- function(x, coef = 1){
##'   return(coef * (0.05 + sqrt(abs(x)*20/(2*pi))/10))
##' }
##' 
##' ## Initial data set
##' nvar <- 1
##' n <- 20
##' X <- matrix(seq(0, 2 * pi, length=n), ncol = 1)
##' mult <- sample(1:10, n, replace = TRUE)
##' X <- rep(X, mult)
##' Z <- sin(X) + rnorm(length(X), sd = noisefun(X))
##' 
##' ## Initial fit
##' testpts <- matrix(seq(0, 2*pi, length = 10*n), ncol = 1)
##' model <- model_init <- mleHetGP(X = X, Z = Z, lower = rep(0.1, nvar), 
##'   upper = rep(50, nvar), maxit = 1000)
##'
##' ## Visualizing initial predictive surface
##' preds <- predict(x = testpts, model_init) 
##' plot(X, Z)
##' lines(testpts, preds$mean, col = "red")
##' 
##' ## 10 fast update steps
##' nsteps <- 5
##' npersteps <- 10
##' for(i in 1:nsteps){
##'   newIds <- sort(sample(1:(10*n), npersteps))
##'   
##'   newX <- testpts[newIds, drop = FALSE] 
##'   newZ <- sin(newX) + rnorm(length(newX), sd = noisefun(newX))
##'   points(newX, newZ, col = "blue", pch = 20)
##'   model <- update(object = model, Xnew = newX, Znew = newZ)
##'   X <- c(X, newX)
##'   Z <- c(Z, newZ)
##'   plot(X, Z)
##'   print(model$nit_opt)
##' }
##'
##' ## Final predictions after 10 updates
##' p_fin <- predict(x=testpts, model) 
##'
##' ## Visualizing the result by augmenting earlier plot
##' lines(testpts, p_fin$mean, col = "blue")
##' lines(testpts, qnorm(0.05, p_fin$mean, sqrt(p_fin$sd2)), col = "blue", lty = 2)
##' lines(testpts, qnorm(0.95, p_fin$mean, sqrt(p_fin$sd2)), col = "blue", lty = 2)
##' lines(testpts, qnorm(0.05, p_fin$mean, sqrt(p_fin$sd2 + p_fin$nugs)), 
##'   col = "blue", lty = 3)
##' lines(testpts, qnorm(0.95, p_fin$mean, sqrt(p_fin$sd2 + p_fin$nugs)), 
##'   col = "blue", lty = 3)
##' 
##' ## Now compare to what you would get if you did a full batch fit instead
##' model_direct <-  mleHetGP(X = X, Z = Z, maxit = 1000,
##'                           lower = rep(0.1, nvar), upper = rep(50, nvar),
##'                           init = list(theta = model_init$theta, k_theta_g = model_init$k_theta_g))
##' p_dir <- predict(x = testpts, model_direct)
##' print(model_direct$nit_opt)
##' lines(testpts, p_dir$mean, col = "green")
##' lines(testpts, qnorm(0.05, p_dir$mean, sqrt(p_dir$sd2)), col = "green", 
##'   lty = 2)
##' lines(testpts, qnorm(0.95, p_dir$mean, sqrt(p_dir$sd2)), col = "green", 
##'   lty = 2)
##' lines(testpts, qnorm(0.05, p_dir$mean, sqrt(p_dir$sd2 + p_dir$nugs)), 
##'   col = "green", lty = 3)
##' lines(testpts, qnorm(0.95, p_dir$mean, sqrt(p_dir$sd2 + p_dir$nugs)), 
##'   col = "green", lty = 3)
##' lines(testpts, sin(testpts), col = "red", lty = 2)
##' 
##' ## Compare outputs
##' summary(model_init)
##' summary(model)
##' summary(model_direct)
##' 
## ' ##------------------------------------------------------------
## ' ## Example 2: update keeping hyperparameters
## ' ##------------------------------------------------------------
## ' 
## ' nvar <- 2
## ' n <- 10
## ' X <- matrix(runif(nvar*n), ncol = nvar)
## ' mult <- sample(1:10, n, replace = TRUE)
## ' X <- X[rep(1:n, times = mult),, drop = FALSE]
## ' Z <- sin(rowSums(X)) + rnorm(nrow(X), sd = noisefun(rowSums(X)))
## ' 
## ' model_init <- mleHetGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar), known = list(beta = 0))
## ' 
## ' n_new <- 10
## ' Xnew <- matrix(runif(nvar * n_new), ncol = nvar)
## ' Znew <- sin(rowSums(Xnew)) + rnorm(nrow(Xnew), sd = noisefun(rowSums(X)))
## ' 
## ' model_up <- update(model_init, Xnew = Xnew, Znew = Znew, maxit = 0)
## ' 
## ' model_ref <- mleHetGP(X = rbind(X, Xnew), Z = c(Z, Znew), upper = model_init$theta + 1e-8, lower = model_init$theta,
## ' noiseControl = list(g_bounds = c(model_init$g, model_init$g + 1e-8)),
## ' known = list(beta = 0, theta = model_init$theta, k_theta_g = model_init$k_theta_g, g = model_init$g, theta_g = model_init$theta_g,
## ' Delta = c(model_init$Delta, log(predict(model_init, Xnew, nugs.only = T)$nugs/model_init$nu_hat))))
## ' 
## ' range(model_ref$Ki - model_up$Ki)
## ' K_ref <- solve(cov_gen(model_up$X0, theta = model_up$theta, type = model_up$covtype) + diag(model_up$Lambda/model_up$mult + 1e-8))
update.hetGP <- function(object, Xnew, Znew, ginit = 1e-2, lower = NULL, upper = NULL, noiseControl = NULL, settings = NULL,
                         known = NULL, maxit = 100, method = 'quick', ...){
  
  # first reduce Xnew/Znew in case of potential replicates
  newdata <- find_reps(Xnew, Znew, normalize = FALSE, rescale = FALSE)
  
  # 'mixed' requires new data
  if(any(is.na(Znew))) method = 'quick'
  
  # copy object for update 
  m_new <- object
  
  # check if newdata contains designs already in X0
  if(any(duplicated(rbind(m_new$X0, newdata$X0)))){
    id_exists <- NULL
    for(i in 1:nrow(newdata$X0)){
      tmp <- duplicated(rbind(newdata$X0[i,], m_new$X0))
      if(any(tmp)){
        id_exists <- c(id_exists, i)
        id_X0 <- which(tmp) - 1
        m_new$Z0[id_X0] <- (m_new$mult[id_X0] * m_new$Z0[id_X0] + newdata$Z0[i] * newdata$mult[i])/(m_new$mult[id_X0] + newdata$mult[i])
        idZ <- cumsum(m_new$mult)
        m_new$Z <- append(m_new$Z, newdata$Zlist[[i]], after = idZ[id_X0])
        
        ## Inverse matrices are updated if MLE is not performed 
        if(maxit == 0){
          m_new$Ki <- update_Ki_rep(id_X0, m_new, nrep = newdata$mult[i])
          m_new$Kgi <- update_Kgi_rep(id_X0, m_new, nrep = newdata$mult[i])
        }
        
        m_new$mult[id_X0] <- m_new$mult[id_X0] + newdata$mult[i]
        
        ### Update Delta value depending on the selected scheme
        
        # if method == 'quick': nothing to be done for replicates
        
        # use object to use previous predictions of Lambda/mean
        if(method == 'mixed'){
          # model predictions
          delta_loo <- LOO_preds_nugs(object, id_X0) ## LOO mean and variance at X0[id_X0,] (only variance is used)
          
          # empirical estimates
          sd2h <- mean((m_new$Z[idZ[id_X0]:(idZ[id_X0] + m_new$mult[id_X0] - 1)] - predict(object, newdata$X0[i,,drop = FALSE])$mean)^2)/object$nu_hat
          sd2sd2h <- 2*sd2h^2/m_new$mult[id_X0] #variance of the estimator of the variance
          
          if(object$logN){
            ## Add correction terms for the log transform
            sd2h <- log(sd2h) - digamma((m_new$mult[id_X0] - 1)/2) - log(2) + log(m_new$mult[id_X0] - 1)
            sd2sd2h <- trigamma((m_new$mult[id_X0] - 1)/2)
          }
          
          delta_loo$mean <- m_new$Delta[id_X0]
          newdelta_sd2 <- 1/(1/delta_loo$sd2 + 1/sd2sd2h)
          new_delta_mean <- (delta_loo$mean/delta_loo$sd2 + sd2h/sd2sd2h) * newdelta_sd2 
          newdata$Delta[id_X0] <- new_delta_mean
        }
        
      } 
    }
    
    # remove duplicates now
    newdata$X0 <- newdata$X0[-id_exists,,drop = FALSE]
    newdata$Z0 <- newdata$Z0[-id_exists]
    newdata$mult <- newdata$mult[-id_exists]
    newdata$Zlist <- newdata$Zlist[-id_exists]
  }
  
  ## Now deal with new data
  if(nrow(newdata$X0) > 0){
    
    if(method == 'quick'){
      delta_new_init <- predict(x = newdata$X0, object, nugs_only = T)$nugs/object$nu_hat
      if(object$logN) delta_new_init <- log(delta_new_init)
    }
    
    if(method == 'mixed'){
      delta_new_init <- rep(NA, nrow(newdata$X0))
      pred_deltas <- predict(x = newdata$X0, object = object, noise.var = T)
      
      for(i in 1:nrow(newdata$X0)){
        sd2h <- mean((newdata$Zlist[[i]] - predict(object, newdata$X0[i,,drop = FALSE])$mean)^2)/object$nu_hat
        sd2sd2h <- 2*sd2h^2/newdata$mult[i] #variance of the estimator of the variance
        
        pdelta <- pred_deltas$nugs[i]/object$nu_hat
        if(object$logN){
          pdelta <- log(pdelta)
          
          ## Correction terms for the log transform
          sd2h <- log(sd2h) - digamma(newdata$mult[i]/2) - log(2) + log(newdata$mult[i])
          sd2sd2h <- trigamma(newdata$mult[i]/2)
        }
        
        newdelta_sd2 <- 1/(object$nu_hat/pred_deltas$sd2var[i] + 1/sd2sd2h)
        delta_new_init[i] <- (object$nu_hat*pdelta/pred_deltas$sd2var[i] + sd2h/sd2sd2h) * newdelta_sd2
        
      }
    }
    
    if(maxit == 0){
      for(i in 1:nrow(newdata$X0)){
        if(object$logN){ ##   Lambda = exp(Delta)
          m_new$Ki <- update_Ki(newdata$X0[i,,drop = F], m_new, nrep = newdata$mult[i], new_lambda = exp(delta_new_init[i]))
        }else{
          m_new$Ki <- update_Ki(newdata$X0[i,,drop = F], m_new, new_lambda = delta_new_init[i], nrep = newdata$mult[i])
        }
        m_new$Kgi <- update_Kgi(newdata$X0[i,,drop = F], m_new, nrep = newdata$mult[i])
        m_new$X0 <- rbind(m_new$X0, newdata$X0[i,, drop = F])
      }
      if(object$logN){ ##   Lambda = exp(Delta)
        m_new$Lambda <- c(m_new$Lambda, exp(delta_new_init))
      }else{
        m_new$Lambda <- c(m_new$Lambda, delta_new_init)
      }
      
    }else{
      m_new$X0 <- rbind(m_new$X0, newdata$X0)
    }
    
    m_new$Z0 <- c(m_new$Z0, newdata$Z0)
    m_new$mult <- c(m_new$mult, newdata$mult)
    m_new$Z <- c(m_new$Z, unlist(newdata$Zlist))
    m_new$Delta <- c(m_new$Delta, delta_new_init)
    
  }
  
  if(maxit == 0){
    m_new$nit_opt <- 0
    m_new$msg <- "Not optimized \n"
    
  }else{
    
    if(is.null(upper)) upper <- object$used_args$upper
    if(is.null(lower)) lower <- object$used_args$lower
    if(is.null(noiseControl)){
      noiseControl <- object$used_args$noiseControl
      noiseControl$lowerDelta <- noiseControl$upperDelta <- NULL ## must be given to noiseControl in update
    } 
    if(is.null(settings)) settings <- object$used_args$settings
    if(is.null(known)) known <- object$used_args$known
    
    m_new <- mleHetGP(X = list(X0 = m_new$X0, Z0 = m_new$Z0, mult = m_new$mult), Z = m_new$Z,
                      noiseControl = noiseControl, lower = lower, upper = upper, covtype = object$covtype, settings = settings,
                      init = list(theta = object$theta, theta_g = object$theta_g, k_theta_g = object$k_theta_g, Delta = m_new$Delta, g = max(object$g, ginit)),
                      known = known, eps = object$eps, maxit = maxit)
  }
  return(m_new)
}



##' Update existing \code{homGP} model with new observations
##' @title Fast \code{homGP}-update
##' @param object initial model of class \code{homGP} 
##' @param Xnew matrix of new design locations; \code{ncol(Xnew)} must match the input dimension encoded in object
##' @param Znew vector new observations at those new design locations, of length \code{nrow(X)}. \code{NA}s can be passed, see Details
##' @param lower,upper,noiseControl,known optional bounds for MLE optimization, see \code{\link[hetGP]{mleHomGP}}.
##' If not provided, they are extracted from the existing model 
##' @param maxit maximum number of iterations for the internal L-BFGS-B optimization method; see \code{\link{optim}} for more details
##' @param ... no other argument for this method.
##' @details 
##' In case hyperparameters need not be updated, \code{maxit} can be set to \code{0}. 
##' In this case it is possible to pass \code{NA}s in \code{Znew}, then the model can still be used to provide updated variance predictions.
##' @export
##' @method update homGP
##' @examples 
##' \dontrun{
##' ##------------------------------------------------------------
##' ## Example : Sequential Homoskedastic GP modeling 
##' ##------------------------------------------------------------
##' set.seed(42)
##'
##' ## Spatially varying noise function
##' noisefun <- function(x, coef = 1){
##'   return(coef * (0.05 + sqrt(abs(x)*20/(2*pi))/10))
##' }
##' 
##' nvar <- 1
##' n <- 10
##' X <- matrix(seq(0, 2 * pi, length=n), ncol = 1)
##' mult <- sample(1:10, n)
##' X <- rep(X, mult)
##' Z <- sin(X) + rnorm(length(X), sd = noisefun(X))
##' 
##' testpts <- matrix(seq(0, 2*pi, length = 10*n), ncol = 1)
##' model <- model_init <- mleHomGP(X = X, Z = Z,
##'                                 lower = rep(0.1, nvar), upper = rep(50, nvar))
##' preds <- predict(x = testpts, object = model_init) 
##' plot(X, Z)
##' lines(testpts, preds$mean, col = "red")
##' 
##' 
##' nsteps <- 10
##' for(i in 1:nsteps){
##'   newIds <- sort(sample(1:(10*n), 10))
##'   
##'   newX <- testpts[newIds, drop = FALSE] 
##'   newZ <- sin(newX) + rnorm(length(newX), sd = noisefun(newX))
##'   points(newX, newZ, col = "blue", pch = 20)
##'   model <- update(object = model, newX, newZ)
##'   X <- c(X, newX)
##'   Z <- c(Z, newZ)
##'   plot(X, Z)
##'   print(model$nit_opt)
##' }
##' p_fin <- predict(x = testpts, object = model) 
##' lines(testpts, p_fin$mean, col = "blue")
##' lines(testpts, qnorm(0.05, p_fin$mean, sqrt(p_fin$sd2)), col = "blue", lty = 2)
##' lines(testpts, qnorm(0.95, p_fin$mean, sqrt(p_fin$sd2)), col = "blue", lty = 2)
##' lines(testpts, qnorm(0.05, p_fin$mean, sqrt(p_fin$sd2 + p_fin$nugs)),
##'       col = "blue", lty = 3)
##' lines(testpts, qnorm(0.95, p_fin$mean, sqrt(p_fin$sd2 + p_fin$nugs)),
##'       col = "blue", lty = 3)
##' 
##' model_direct <-  mleHomGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar))
##' p_dir <- predict(x = testpts, object = model_direct)
##' print(model_direct$nit_opt)
##' lines(testpts, p_dir$mean, col = "green")
##' lines(testpts, qnorm(0.05, p_dir$mean, sqrt(p_dir$sd2)), col = "green", lty = 2)
##' lines(testpts, qnorm(0.95, p_dir$mean, sqrt(p_dir$sd2)), col = "green", lty = 2)
##' lines(testpts, qnorm(0.05, p_dir$mean, sqrt(p_dir$sd2 + p_dir$nugs)),
##'       col = "green", lty = 3)
##' lines(testpts, qnorm(0.95, p_dir$mean, sqrt(p_dir$sd2 + p_dir$nugs)),
##'       col = "green", lty = 3)
##' 
##' lines(testpts, sin(testpts), col = "red", lty = 2)
##' 
##' ## Compare outputs
##' summary(model_init)
##' summary(model)
##' summary(model_direct)
##' 
##' 
## ' ##------------------------------------------------------------
## ' ## Example 2: update keeping hyperparameters
## ' ##------------------------------------------------------------
## ' 
## ' nvar <- 2
## ' n <- 10
## ' X <- matrix(runif(nvar*n), ncol = nvar)
## ' mult <- sample(1:10, n, replace = TRUE)
## ' X <- X[rep(1:n, times = mult),]
## ' Z <- sin(rowSums(X)) + rnorm(nrow(X), sd = noisefun(rowSums(X)))
## ' 
## ' model_init <- mleHomGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar))
## ' 
## ' n_new <- 10
## ' Xnew <- matrix(runif(nvar * n_new), ncol = nvar)
## ' Znew <- sin(rowSums(Xnew)) + rnorm(nrow(Xnew), sd = noisefun(rowSums(X)))
## ' 
## ' model_up <- update(model_init, Xnew = Xnew, Znew = Znew, maxit = 0, beta = 0)
## ' 
## ' model_ref <- mleHomGP(X = rbind(X, Xnew), Z = c(Z, Znew), upper = model_init$theta + 1e-8, lower = model_init$theta,
## '                       noiseControl = list(g_bounds = c(model_init$g, model_init$g + 1e-8)), beta = 0)
## ' range(model_ref$Ki - model_up$Ki)
##' }
update.homGP <- function(object, Xnew, Znew = NULL, lower = NULL, upper = NULL, noiseControl = NULL, known = NULL, maxit = 100, ...){
  # first reduce Xnew/Znew in case of potential replicates
  newdata <- find_reps(Xnew, Znew, rescale = FALSE, normalize = FALSE)
  
  # check if newdata contains designs already in X0
  if(any(duplicated(rbind(object$X0, newdata$X0)))){
    id_exists <- NULL
    for(i in 1:nrow(newdata$X0)){
      tmp <- duplicated(rbind(newdata$X0[i,], object$X0))
      if(any(tmp)){
        id_exists <- c(id_exists, i)
        id_X0 <- which(tmp) - 1
        object$Z0[id_X0] <- (object$mult[id_X0] * object$Z0[id_X0] + newdata$Z0[i] * newdata$mult[i])/(object$mult[id_X0] + newdata$mult[i])
        idZ <- cumsum(object$mult)
        object$Z <- append(object$Z, newdata$Zlist[[i]], after = idZ[id_X0])
        
        if(maxit == 0){
          object$Ki <- update_Ki_rep(id_X0, object, nrep = newdata$mult[i])
          object$nit_opt <- 0
          object$msg <- "Not optimized \n"
        } 
        object$mult[id_X0] <- object$mult[id_X0] + newdata$mult[i]
      } 
    }
    # remove duplicates now
    newdata$X0 <- newdata$X0[-id_exists,,drop = FALSE]
    newdata$Z0 <- newdata$Z0[-id_exists]
    newdata$mult <- newdata$mult[-id_exists]
    newdata$Zlist <- newdata$Zlist[-id_exists]
  }
  
  if(nrow(newdata$X0) > 0 && maxit == 0){
    for(i in 1:nrow(newdata$X0)){
      object$Ki <- update_Ki(newdata$X0[i,,drop = F], object, nrep = newdata$mult[i])
      object$X0 <- rbind(object$X0, newdata$X0[i,])
      object$Z0 <- c(object$Z0, newdata$Z0[i])
      object$mult <- c(object$mult, newdata$mult[i])
      object$Z <- c(object$Z, newdata$Zlist[[i]])
    }
    
    object$nit_opt <- 0
    object$msg <- "Not optimized \n"
  }
  
  if(maxit> 0){
    
    
    if(is.null(upper)) upper <- object$used_args$upper
    if(is.null(lower)) lower <- object$used_args$lower
    if(is.null(noiseControl)) noiseControl <- object$used_args$noiseControl
    if(is.null(known)) known <- object$used_args$known
    
    init <- NULL
    if(is.null(known$theta)) init <- list(theta = object$theta)
    if(is.null(known$g)) init <- c(init, list(g = object$g))
    
    object <- mleHomGP(X = list(X0 = rbind(object$X0, newdata$X0), Z0 = c(object$Z0, newdata$Z0), mult = c(object$mult, newdata$mult)), Z = c(object$Z, unlist(newdata$Zlist)), 
                       lower = lower, upper = upper, noiseControl = noiseControl, covtype = object$covtype, 
                       init = init, 
                       known = known, eps = object$eps, maxit = maxit)
  }
  
  return(object)
  
}

##' Update existing \code{homTP} model with new observations
##' @title Fast \code{homTP}-update
##' @param object initial model of class \code{homTP} 
##' @param Xnew matrix of new design locations; \code{ncol(Xnew)} must match the input dimension encoded in object
##' @param Znew vector new observations at those new design locations, of length \code{nrow(X)}. \code{NA}s can be passed, see Details
##' @param lower,upper,noiseControl,known optional bounds for MLE optimization, see \code{\link[hetGP]{mleHomTP}}.
##' If not provided, they are extracted from the existing model 
##' @param maxit maximum number of iterations for the internal L-BFGS-B optimization method; see \code{\link{optim}} for more details
##' @param ... no other argument for this method.
##' @details 
##' In case hyperparameters need not be updated, \code{maxit} can be set to \code{0}. 
##' In this case it is possible to pass \code{NA}s in \code{Znew}, then the model can still be used to provide updated variance predictions.
##' @export
##' @method update homTP
##' @examples 
##' \dontrun{
##' ##------------------------------------------------------------
##' ## Example : Sequential Homoskedastic TP moding 
##' ##------------------------------------------------------------
##' set.seed(42)
##'
##' ## Spatially varying noise function
##' noisefun <- function(x, coef = 1){
##'   return(coef * (0.05 + sqrt(abs(x)*20/(2*pi))/10))
##' }
##' 
##' df_noise <- 3
##' nvar <- 1
##' n <- 10
##' X <- matrix(seq(0, 2 * pi, length=n), ncol = 1)
##' mult <- sample(1:50, n, replace = TRUE)
##' X <- rep(X, mult)
##' Z <- sin(X) + noisefun(X) * rt(length(X), df = df_noise)
##' 
##' testpts <- matrix(seq(0, 2*pi, length = 10*n), ncol = 1)
##' mod <- mod_init <- mleHomTP(X = X, Z = Z, covtype = "Matern5_2",
##'                                 lower = rep(0.1, nvar), upper = rep(50, nvar))
##' preds <- predict(x = testpts, object = mod_init) 
##' plot(X, Z)
##' lines(testpts, preds$mean, col = "red")
##' 
##' 
##' nsteps <- 10
##' for(i in 1:nsteps){
##'   newIds <- sort(sample(1:(10*n), 5))
##'   
##'   newX <- testpts[rep(newIds, times = sample(1:50, length(newIds), replace = TRUE)), drop = FALSE] 
##'   newZ <- sin(newX) + noisefun(newX) * rt(length(newX), df = df_noise)
##'   points(newX, newZ, col = "blue", pch = 20)
##'   mod <- update(object = mod, newX, newZ)
##'   X <- c(X, newX)
##'   Z <- c(Z, newZ)
##'   plot(X, Z)
##'   print(mod$nit_opt)
##' }
##' p_fin <- predict(x = testpts, object = mod) 
##' lines(testpts, p_fin$mean, col = "blue")
##' lines(testpts, p_fin$mean + sqrt(p_fin$sd2) * qt(0.05, df = mod$nu + length(Z)),
##'       col = "blue", lty = 2)
##' lines(testpts, p_fin$mean + sqrt(p_fin$sd2) * qt(0.95, df = mod$nu + length(Z)),
##'       col = "blue", lty = 2)
##' lines(testpts, p_fin$mean + sqrt(p_fin$sd2 + p_fin$nugs) * qt(0.05, df = mod$nu + length(Z)),
##'       col = "blue", lty = 3)
##' lines(testpts, p_fin$mean + sqrt(p_fin$sd2 + p_fin$nugs) * qt(0.95, df = mod$nu + length(Z)),
##'       col = "blue", lty = 3)
##' 
##' mod_dir <-  mleHomTP(X = X, Z = Z, covtype = "Matern5_2",
##'                           lower = rep(0.1, nvar), upper = rep(50, nvar))
##' p_dir <- predict(x = testpts, object = mod_dir)
##' print(mod_dir$nit_opt)
##' lines(testpts, p_dir$mean, col = "green")
##' lines(testpts, p_dir$mean + sqrt(p_dir$sd2) * qt(0.05, df = mod_dir$nu + length(Z)),
##'       col = "green", lty = 2)
##' lines(testpts, p_dir$mean + sqrt(p_dir$sd2) * qt(0.95, df = mod_dir$nu + length(Z)),
##'       col = "green", lty = 2)
##' lines(testpts, p_dir$mean + sqrt(p_dir$sd2 + p_dir$nugs) * qt(0.05, df = mod_dir$nu + length(Z)),
##'       col = "green", lty = 3)
##' lines(testpts, p_dir$mean + sqrt(p_dir$sd2 + p_dir$nugs) * qt(0.95, df = mod_dir$nu + length(Z)),
##'       col = "green", lty = 3)
##' 
##' lines(testpts, sin(testpts), col = "red", lty = 2)
##' 
##' ## Compare outputs
##' summary(mod_init)
##' summary(mod)
##' summary(mod_dir)
##' 
##' 
## ' ##------------------------------------------------------------
## ' ## Example 2: update keeping hyperparameters
## ' ##------------------------------------------------------------
## ' 
## ' nvar <- 2
## ' n <- 10
## ' X <- matrix(runif(nvar*n), ncol = nvar)
## ' mult <- sample(1:10, n, replace = TRUE)
## ' X <- X[rep(1:n, times = mult),]
## ' Z <- sin(rowSums(X)) + rnorm(nrow(X), sd = noisefun(rowSums(X)))
## ' 
## ' model_init <- mleHomTP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar))
## ' 
## ' n_new <- 10
## ' Xnew <- matrix(runif(nvar * n_new), ncol = nvar)
## ' Znew <- sin(rowSums(Xnew)) + rnorm(nrow(Xnew), sd = noisefun(rowSums(X)))
## ' 
## ' model_up <- update(model_init, Xnew = Xnew, Znew = Znew, maxit = 0, beta = 0)
## ' 
## ' model_ref <- mleHomTP(X = rbind(X, Xnew), Z = c(Z, Znew), upper = model_init$theta + 1e-8, lower = model_init$theta,
## '                       noiseControl = list(g_bounds = c(model_init$g, model_init$g + 1e-8)), beta = 0)
## ' range(model_ref$Ki - model_up$Ki)
##' }
update.homTP <- function(object, Xnew, Znew = NULL, lower = NULL, upper = NULL, noiseControl = NULL, known = NULL, maxit = 100, ...){
  # first reduce Xnew/Znew in case of potential replicates
  newdata <- find_reps(Xnew, Znew, rescale = FALSE, normalize = FALSE)
  
  # check if newdata contains designs already in X0
  if(any(duplicated(rbind(object$X0, newdata$X0)))){
    id_exists <- NULL
    for(i in 1:nrow(newdata$X0)){
      tmp <- duplicated(rbind(newdata$X0[i,], object$X0))
      if(any(tmp)){
        id_exists <- c(id_exists, i)
        id_X0 <- which(tmp) - 1
        object$Z0[id_X0] <- (object$mult[id_X0] * object$Z0[id_X0] + newdata$Z0[i] * newdata$mult[i])/(object$mult[id_X0] + newdata$mult[i])
        idZ <- cumsum(object$mult)
        object$Z <- append(object$Z, newdata$Zlist[[i]], after = idZ[id_X0])
        
        if(maxit == 0){
          object$Ki <- update_Ki_rep(id_X0, object, nrep = newdata$mult[i])
          object$nit_opt <- 0
          object$msg <- "Not optimized \n"
          if(!any(is.na(Znew))){
            object$psi <- (crossprod(object$Z - object$beta0) - crossprod((object$Z0 - object$beta0) * object$mult, object$Z0 - object$beta0))/object$g
            object$psi <- drop(object$psi) + drop(crossprod(object$Z0 - object$beta0, object$Ki) %*% (object$Z0 - object$beta0))
          } 
          
        } 
        object$mult[id_X0] <- object$mult[id_X0] + newdata$mult[i]
      } 
    }
    # remove duplicates now
    newdata$X0 <- newdata$X0[-id_exists,,drop = FALSE]
    newdata$Z0 <- newdata$Z0[-id_exists]
    newdata$mult <- newdata$mult[-id_exists]
    newdata$Zlist <- newdata$Zlist[-id_exists]
  }
  
  if(nrow(newdata$X0) > 0 && maxit == 0){
    for(i in 1:nrow(newdata$X0)){
      object$Ki <- update_Ki(newdata$X0[i,,drop = F], object, nrep = newdata$mult[i])
      object$X0 <- rbind(object$X0, newdata$X0[i,])
      object$Z0 <- c(object$Z0, newdata$Z0[i])
      object$mult <- c(object$mult, newdata$mult[i])
      object$Z <- c(object$Z, newdata$Zlist[[i]])
      if(!any(is.na(Znew))){
        object$psi <- (crossprod(object$Z - object$beta0) - crossprod((object$Z0 - object$beta0) * object$mult, object$Z0 - object$beta0))/object$g
        object$psi <- drop(object$psi) + drop(crossprod(object$Z0 - object$beta0, object$Ki) %*% (object$Z0 - object$beta0))
      } 
    }
    
    object$nit_opt <- 0
    object$msg <- "Not optimized \n"
  }
  
  if(maxit> 0){
    
    
    if(is.null(upper)) upper <- object$used_args$upper
    if(is.null(lower)) lower <- object$used_args$lower
    if(is.null(noiseControl)) noiseControl <- object$used_args$noiseControl
    if(is.null(known)) known <- object$used_args$known
    
    init <- NULL
    if(is.null(known$theta)) init <- list(theta = object$theta)
    if(is.null(known$g)) init <- c(init, list(g = object$g))
    if(is.null(known$nu)) init <- c(init, list(nu = object$nu))
    if(is.null(known$sigma2)) init <- c(init, list(sigma2 = object$sigma2))
    
    object <- mleHomTP(X = list(X0 = rbind(object$X0, newdata$X0), Z0 = c(object$Z0, newdata$Z0), mult = c(object$mult, newdata$mult)), Z = c(object$Z, unlist(newdata$Zlist)), 
                       lower = lower, upper = upper, noiseControl = noiseControl, covtype = object$covtype, 
                       init = init, 
                       known = known, eps = object$eps, maxit = maxit)
  }
  
  return(object)
  
}

##' Fast update of existing \code{hetTP} model with new observations. 
##' @title Update \code{"hetTP"}-class model fit with new observations
##' @param object previously fit \code{"hetTP"}-class model
##' @param Xnew matrix of new design locations; \code{ncol(Xnew)} must match the input dimension encoded in \code{object}
##' @param Znew vector new observations at those design locations, of length \code{nrow(X)}. \code{NA}s can be passed, see Details
##' @param ginit minimal value of the smoothing parameter (i.e., nugget of the noise process) for optimization initialisation.
##' It is compared to the \code{g} hyperparameter in the object.   
##' @param lower,upper,noiseControl,settings,known optional bounds for mle optimization, see \code{\link[hetGP]{mleHetTP}}. 
##' If not provided, they are extracted from the existing model 
##' @param maxit maximum number of iterations for the internal L-BFGS-B optimization method; see \code{\link{optim}} for more details
##' @param method one of \code{"quick"}, \code{"mixed"} see Details.
##' @param ... no other argument for this method.
##' @details
##' 
##' The update can be performed with or without re-estimating hyperparameter.
##' In the first case, \code{\link[hetGP]{mleHetTP}} is called, based on previous values for initialization. 
##' The only missing values are the latent variables at the new points, that are initialized based on two possible update schemes in \code{method}:
##' \itemize{
##'   \item \code{"quick"} the new delta value is the predicted nugs value from the previous noise model;
##'   \item \code{"mixed"} new values are taken as the barycenter between prediction given by the noise process and empirical variance. 
##' }
##' The subsequent number of MLE computations can be controlled with \code{maxit}.
##' 
##' In case hyperparameters need not be updated, \code{maxit} can be set to \code{0}. 
##' In this case it is possible to pass \code{NA}s in \code{Znew}, then the model can still be used to provide updated variance predictions.
##' 
##' @export
##' @method update hetTP
##' @importFrom stats rnorm
##' @examples 
##' ##------------------------------------------------------------
##' ## Sequential update example
##' ##------------------------------------------------------------
##' set.seed(42)
##'
##' ## Spatially varying noise function
##' noisefun <- function(x, coef = 1){
##'   return(coef * (0.05 + sqrt(abs(x)*20/(2*pi))/10))
##' }
##' 
##' ## Initial data set
##' nvar <- 1
##' n <- 20
##' X <- matrix(seq(0, 2 * pi, length=n), ncol = 1)
##' mult <- sample(1:10, n, replace = TRUE)
##' X <- rep(X, mult)
##' Z <- sin(X) + noisefun(X) * rt(length(X), df = 10)
##' 
##' ## Initial fit
##' testpts <- matrix(seq(0, 2*pi, length = 10*n), ncol = 1)
##' model <- model_init <- mleHetTP(X = X, Z = Z, lower = rep(0.1, nvar), 
##'   upper = rep(50, nvar), maxit = 1000)
##'
##' ## Visualizing initial predictive surface
##' preds <- predict(x = testpts, model_init) 
##' plot(X, Z)
##' lines(testpts, preds$mean, col = "red")
##' 
##' ## 10 fast update steps
##' nsteps <- 5
##' npersteps <- 10
##' for(i in 1:nsteps){
##'   newIds <- sort(sample(1:(10*n), npersteps))
##'   
##'   newX <- testpts[newIds, drop = FALSE] 
##'   newZ <- sin(newX) + noisefun(newX) * rt(length(newX), df = 10)
##'   points(newX, newZ, col = "blue", pch = 20)
##'   model <- update(object = model, Xnew = newX, Znew = newZ)
##'   X <- c(X, newX)
##'   Z <- c(Z, newZ)
##'   plot(X, Z)
##'   print(model$nit_opt)
##' }
##'
##' ## Final predictions after 10 updates
##' p_fin <- predict(x=testpts, model) 
##'
##' ## Visualizing the result by augmenting earlier plot
##' lines(testpts, p_fin$mean, col = "blue")
##' lines(testpts, qnorm(0.05, p_fin$mean, sqrt(p_fin$sd2)), col = "blue", lty = 2)
##' lines(testpts, qnorm(0.95, p_fin$mean, sqrt(p_fin$sd2)), col = "blue", lty = 2)
##' lines(testpts, qnorm(0.05, p_fin$mean, sqrt(p_fin$sd2 + p_fin$nugs)), 
##'   col = "blue", lty = 3)
##' lines(testpts, qnorm(0.95, p_fin$mean, sqrt(p_fin$sd2 + p_fin$nugs)), 
##'   col = "blue", lty = 3)
##' 
##' ## Now compare to what you would get if you did a full batch fit instead
##' model_direct <-  mleHetTP(X = X, Z = Z, maxit = 1000,
##'                           lower = rep(0.1, nvar), upper = rep(50, nvar),
##'                           init = list(theta = model_init$theta, k_theta_g = model_init$k_theta_g))
##' p_dir <- predict(x = testpts, model_direct)
##' print(model_direct$nit_opt)
##' lines(testpts, p_dir$mean, col = "green")
##' lines(testpts, qnorm(0.05, p_dir$mean, sqrt(p_dir$sd2)), col = "green", 
##'   lty = 2)
##' lines(testpts, qnorm(0.95, p_dir$mean, sqrt(p_dir$sd2)), col = "green", 
##'   lty = 2)
##' lines(testpts, qnorm(0.05, p_dir$mean, sqrt(p_dir$sd2 + p_dir$nugs)), 
##'   col = "green", lty = 3)
##' lines(testpts, qnorm(0.95, p_dir$mean, sqrt(p_dir$sd2 + p_dir$nugs)), 
##'   col = "green", lty = 3)
##' lines(testpts, sin(testpts), col = "red", lty = 2)
##' 
##' ## Compare outputs
##' summary(model_init)
##' summary(model)
##' summary(model_direct)
##' 
## ' ##------------------------------------------------------------
## ' ## Example 2: update keeping hyperparameters
## ' ##------------------------------------------------------------
## ' 
## ' nvar <- 2
## ' n <- 10
## ' X <- matrix(runif(nvar*n), ncol = nvar)
## ' mult <- sample(1:10, n, replace = TRUE)
## ' X <- X[rep(1:n, times = mult),, drop = FALSE]
## ' Z <- sin(rowSums(X)) + rnorm(nrow(X), sd = noisefun(rowSums(X)))
## ' 
## ' model_init <- mleHetTP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar), known = list(beta0 = 0))
## ' 
## ' n_new <- 10
## ' Xnew <- matrix(runif(nvar * n_new), ncol = nvar)
## ' Znew <- sin(rowSums(Xnew)) + rnorm(nrow(Xnew), sd = noisefun(rowSums(X)))
## ' 
## ' model_up <- update(model_init, Xnew = Xnew, Znew = Znew, maxit = 0)
## ' 
## ' model_ref <- mleHetTP(X = rbind(X, Xnew), Z = c(Z, Znew), upper = model_init$theta + 1e-8, lower = model_init$theta,
## ' noiseControl = list(g_bounds = c(model_init$g, model_init$g + 1e-8)),
## ' known = list(beta0 = 0, theta = model_init$theta, k_theta_g = model_init$k_theta_g, g = model_init$g, theta_g = model_init$theta_g,
## '              Delta = c(model_init$Delta, log(predict(model_init, Xnew, nugs.only = T)$nugs)), nu = model_up$nu, sigma2 = model_init$sigma2))
## ' 
## ' range(model_ref$Ki - model_up$Ki)
update.hetTP <- function(object, Xnew, Znew, ginit = 1e-2, lower = NULL, upper = NULL, noiseControl = NULL, settings = NULL,
                         known = NULL, maxit = 100, method = 'quick', ...){
  
  # first reduce Xnew/Znew in case of potential replicates
  newdata <- find_reps(Xnew, Znew, normalize = FALSE, rescale = FALSE)
  
  # 'mixed' requires new data
  if(any(is.na(Znew))) method = 'quick'
  
  # copy object for update 
  m_new <- object
  
  # check if newdata contains designs already in X0
  if(any(duplicated(rbind(m_new$X0, newdata$X0)))){
    id_exists <- NULL
    for(i in 1:nrow(newdata$X0)){
      tmp <- duplicated(rbind(newdata$X0[i,], m_new$X0))
      if(any(tmp)){
        id_exists <- c(id_exists, i)
        id_X0 <- which(tmp) - 1
        m_new$Z0[id_X0] <- (m_new$mult[id_X0] * m_new$Z0[id_X0] + newdata$Z0[i] * newdata$mult[i])/(m_new$mult[id_X0] + newdata$mult[i])
        idZ <- cumsum(m_new$mult)
        m_new$Z <- append(m_new$Z, newdata$Zlist[[i]], after = idZ[id_X0])
        
        m_new$mult[id_X0] <- m_new$mult[id_X0] + newdata$mult[i]
        
        ## Inverse matrices are updated if MLE is not performed 
        if(maxit == 0){
          m_new$Ki <- update_Ki_rep(id_X0, m_new, nrep = newdata$mult[i])
          m_new$Kgi <- update_Kgi_rep(id_X0, m_new, nrep = newdata$mult[i])
          if(!any(is.na(Znew))){
            m_new$psi <- drop(crossprod((m_new$Z - m_new$beta0)/rep(m_new$Lambda, times = m_new$mult), m_new$Z - m_new$beta0) - 
                                crossprod((m_new$Z0 - m_new$beta0) * m_new$mult/m_new$Lambda, m_new$Z0 - m_new$beta0))
            m_new$psi <- m_new$psi + drop(crossprod(m_new$Z0 - m_new$beta0, m_new$Ki) %*% (m_new$Z0 - m_new$beta0))
          } 
        }
        
        ### Update Delta value depending on the selected scheme
        
        # if method == 'quick': nothing to be done for replicates
        
        # use object to use previous predictions of Lambda/mean
        if(method == 'mixed'){
          # model predictions
          delta_loo <- LOO_preds_nugs(object, id_X0) ## LOO mean and variance at X0[id_X0,] (only variance is used)
          
          # empirical estimates
          sd2h <- mean((m_new$Z[idZ[id_X0]:(idZ[id_X0] + m_new$mult[id_X0] - 1)] - predict(object, newdata$X0[i,,drop = FALSE])$mean)^2)
          sd2sd2h <- 2*sd2h^2/m_new$mult[id_X0] #variance of the estimator of the variance
          
          if(object$logN){
            ## Add correction terms for the log transform
            sd2h <- log(sd2h) - digamma((m_new$mult[id_X0] - 1)/2) - log(2) + log(m_new$mult[id_X0] - 1)
            sd2sd2h <- trigamma((m_new$mult[id_X0] - 1)/2)
          }
          
          delta_loo$mean <- m_new$Delta[id_X0]
          newdelta_sd2 <- 1/(1/delta_loo$sd2 + 1/sd2sd2h)
          new_delta_mean <- (delta_loo$mean/delta_loo$sd2 + sd2h/sd2sd2h) * newdelta_sd2 
          newdata$Delta[id_X0] <- new_delta_mean
        }
        
      } 
    }
    
    # remove duplicates now
    newdata$X0 <- newdata$X0[-id_exists,,drop = FALSE]
    newdata$Z0 <- newdata$Z0[-id_exists]
    newdata$mult <- newdata$mult[-id_exists]
    newdata$Zlist <- newdata$Zlist[-id_exists]
  }
  
  ## Now deal with new data
  if(nrow(newdata$X0) > 0){
    
    if(method == 'quick'){
      delta_new_init <- predict(x = newdata$X0, object, nugs_only = T)$nugs
      if(object$logN) delta_new_init <- log(delta_new_init)
    }
    
    if(method == 'mixed'){
      delta_new_init <- rep(NA, nrow(newdata$X0))
      pred_deltas <- predict(x = newdata$X0, object = object, noise.var = T)
      
      for(i in 1:nrow(newdata$X0)){
        sd2h <- mean((newdata$Zlist[[i]] - predict(object, newdata$X0[i,,drop = FALSE])$mean)^2)
        sd2sd2h <- 2*sd2h^2/newdata$mult[i] #variance of the estimator of the variance
        
        pdelta <- pred_deltas$nugs[i]
        if(object$logN){
          pdelta <- log(pdelta)
          
          ## Correction terms for the log transform
          sd2h <- log(sd2h) - digamma(newdata$mult[i]/2) - log(2) + log(newdata$mult[i])
          sd2sd2h <- trigamma(newdata$mult[i]/2)
        }
        
        newdelta_sd2 <- 1/(object$nu_hat/pred_deltas$sd2var[i] + 1/sd2sd2h)
        delta_new_init[i] <- (object$nu_hat*pdelta/pred_deltas$sd2var[i] + sd2h/sd2sd2h) * newdelta_sd2
        
      }
    }
    
    if(maxit == 0){
      for(i in 1:nrow(newdata$X0)){
        if(object$logN){ ##   Lambda = exp(Delta)
          m_new$Ki <- update_Ki(newdata$X0[i,,drop = F], m_new, nrep = newdata$mult[i], new_lambda = exp(delta_new_init[i]))
        }else{
          m_new$Ki <- update_Ki(newdata$X0[i,,drop = F], m_new, new_lambda = delta_new_init[i], nrep = newdata$mult[i])
        }
        m_new$Kgi <- update_Kgi(newdata$X0[i,,drop = F], m_new, nrep = newdata$mult[i])
        m_new$X0 <- rbind(m_new$X0, newdata$X0[i,, drop = F])
      }
      if(object$logN){ ##   Lambda = exp(Delta)
        m_new$Lambda <- c(m_new$Lambda, exp(delta_new_init))
      }else{
        m_new$Lambda <- c(m_new$Lambda, delta_new_init)
      }
      
    }else{
      m_new$X0 <- rbind(m_new$X0, newdata$X0)
    }
    
    m_new$Z0 <- c(m_new$Z0, newdata$Z0)
    m_new$mult <- c(m_new$mult, newdata$mult)
    m_new$Z <- c(m_new$Z, unlist(newdata$Zlist))
    m_new$Delta <- c(m_new$Delta, delta_new_init)
    
    if(maxit == 0 && !any(is.na(Znew))){
      m_new$psi <- drop(crossprod((m_new$Z - m_new$beta0)/rep(m_new$Lambda, times = m_new$mult), m_new$Z - m_new$beta0) - 
                          crossprod((m_new$Z0 - m_new$beta0) * m_new$mult/m_new$Lambda, m_new$Z0 - m_new$beta0))
      m_new$psi <- m_new$psi + drop(crossprod(m_new$Z0 - m_new$beta0, m_new$Ki) %*% (m_new$Z0 - m_new$beta0))
    } 
    
  }
  
  if(maxit == 0){
    m_new$nit_opt <- 0
    m_new$msg <- "Not optimized \n"
    
  }else{
    
    if(is.null(upper)) upper <- object$used_args$upper
    if(is.null(lower)) lower <- object$used_args$lower
    if(is.null(noiseControl)){
      noiseControl <- object$used_args$noiseControl
      noiseControl$lowerDelta <- noiseControl$upperDelta <- NULL ## must be given to noiseControl in update
    } 
    if(is.null(settings)) settings <- object$used_args$settings
    if(is.null(known)) known <- object$used_args$known
    
    m_new <- mleHetTP(X = list(X0 = m_new$X0, Z0 = m_new$Z0, mult = m_new$mult), Z = m_new$Z,
                      noiseControl = noiseControl, lower = lower, upper = upper, covtype = object$covtype, settings = settings,
                      init = list(theta = object$theta, theta_g = object$theta_g, k_theta_g = object$k_theta_g,
                                  Delta = m_new$Delta, g = max(object$g, ginit), nu = object$new, sigma2 = object$sigma2),
                      known = known, eps = object$eps, maxit = maxit)
  }
  return(m_new)
}




## ' Compute the inverse covariance matrix when adding a new design
## ' @param x matrix for the new design
## ' @param model \code{homGP} or \code{hetGP} model
## ' @param new_delta optional vector. In case of a \code{hetGP} model, value of Delta at \code{x}. 
## ' If not provided, it is taken as the prediction of the latent GP. 
## ' For \code{homGP} models, it corresponds to \code{g}. 
## ' @details see e.g., Bobby's lecture on design, for the partition inverse equations
## ' @export
## ' @examples 
## ' \dontrun{
## ' ## Validation
## ' ## 1) run example("mleHetGP", ask = FALSE)
## ' 
## ' nvar <- 2
## ' ntest <- 10
## ' design <- matrix(runif(ntest * nvar), ntest, nvar)
## ' response <- sum(sin(2*pi*design)) + rnorm(ntest)
## ' 
## ' model <- mleHetGP(X = design, Z = response, upper = rep(1, nvar), lower = rep(0.01,nvar))
## ' 
## ' nrep <- 4
## ' xnew <- matrix(runif(nvar), 1)
## ' Kni <- update_Ki(xnew, model, nrep = nrep)
## ' 
## ' Lambdan <- c(model$Lambda, predict(model, x = xnew)$nugs/model$nu_hat)
## ' multn <- c(model$mult, nrep)
## ' Kn <- cov_gen(rbind(model$X0, xnew), theta = model$theta, type = model$covtype) + diag(model$eps + Lambdan/multn)
## ' Kni_ref <- chol2inv(chol(Kn))
## ' print(max(abs(Kni %*% Kn - diag(nrow(model$X0) + 1))))                                 
## ' }
update_Ki <- function(x, model, new_lambda = NULL, nrep = 1){
  if(class(model) %in% c("homTP", "hetTP")){
    kn1 <- model$sigma2 * cov_gen(x, model$X0, theta = model$theta, type = model$covtype)
    if(is.null(new_lambda))
      new_lambda <- predict(object = model, x = x, nugs.only = TRUE)$nugs
    vn <- drop(model$sigma2 - kn1 %*% tcrossprod(model$Ki, kn1)) + new_lambda/nrep + model$eps
  }else{
    kn1 <- cov_gen(x, model$X0, theta = model$theta, type = model$covtype)
    if(is.null(new_lambda))
      new_lambda <- predict(object = model, x = x, nugs.only = TRUE)$nugs/model$nu_hat
    vn <- drop(1 - kn1 %*% tcrossprod(model$Ki, kn1)) + new_lambda/nrep + model$eps
  }
  
  gn <- - tcrossprod(model$Ki, kn1) / vn
  Ki <- model$Ki + tcrossprod(gn) * vn
  
  return(rbind(cbind(Ki, gn), c(gn, 1/vn)))
}

## Same as update_Ki but for Kgi
update_Kgi <- function(x, model, nrep = 1){
  kn1 <- cov_gen(x, model$X0, theta = model$theta_g, type = model$covtype)
  nugsn <- model$g/nrep
  vn <- drop(1 - kn1 %*% tcrossprod(model$Kgi, kn1)) + nugsn/nrep + model$eps
  gn <- - tcrossprod(model$Kgi, kn1) / vn
  Kgi <- model$Kgi + tcrossprod(gn) * vn
  
  return(rbind(cbind(Kgi, gn), c(gn, 1/vn)))
}

## ## Verification
## ## 1) run example("mleHetGP", ask = FALSE)
## 
## Kni <- hetGP:::update_Ki_rep(16, model)
## multn <- model$mult
## multn[16] <- multn[16] + 1
## Kn <- cov_gen(model$X0, theta = model$theta, type = model$covtype) + diag(model$eps + model$Lambda/multn)
## Kni_ref <- chol2inv(chol(Kn))
## print(max(abs(Kni %*% Kn - diag(nrow(model$X0))))) 
##
## update of Ki in case model$X0[id,] is replicated nrep times
update_Ki_rep <- function(id, model, nrep = 1){
  if(is.null(model$Lambda)){
    tmp <- model$g
  }else{
    tmp <- model$Lambda[id]
  }
  B <- tcrossprod(model$Ki[id,]) / ((model$mult[id]*(model$mult[id] + nrep)) / (nrep*tmp) - model$Ki[id, id])
  Ki <- model$Ki + B
  return(Ki)
}

## update of Kgi in case model$X0[id,] is replicated nrep times
update_Kgi_rep <- function(id, model, nrep = 1){
  tmp <- model$g
  B <- tcrossprod(model$Kgi[id,]) / ((model$mult[id]*(model$mult[id] + nrep)) / (nrep*tmp) - model$Kgi[id, id])
  Kgi <- model$Kgi + B
  return(Kgi)
}


## @param model a hetGP model
## @param i index of the point to remove
LOO_preds_nugs <- function(model, i){
  
  model$Kgi <- model$Kgi - matrix(rowSums(model$Kgi), ncol = 1) %*% matrix(rowSums(model$Kgi), nrow = 1) / sum(model$Kgi)
  yih <- model$Delta[i] - (model$Kgi %*% (model$Delta - model$nmean))[i]/model$Kgi[i,i]
  sih <- 1/model$Kgi[i,i] - model$g  
  
  return(list(mean = yih, sd2 = sih))
}

##' Provide leave one out predictions, e.g., for model testing and diagnostics. 
##' This is used in the method plot available on GP and TP models.
##' @title Leave one out predictions
##' @param model \code{homGP} or \code{hetGP} model, TP version is not considered at this point
##' @param ids vector of indices of the unique design point considered (default to all)
##' @return list with mean and variance predictions at x_i assuming this point has not been evaluated
##' @export
##' @note For TP models, \code{psi} is considered fixed.
##' @references
##' O. Dubrule (1983), Cross validation of Kriging in a unique neighborhood, Mathematical Geology 15, 687--699. \cr \cr
##' 
##' F. Bachoc (2013), Cross Validation and Maximum Likelihood estimations of hyper-parameters of Gaussian processes 
##' with model misspecification, Computational Statistics & Data Analysis, 55--69.
##' 
##' @examples
##' set.seed(32)
##' ## motorcycle data
##' library(MASS)
##' X <- matrix(mcycle$times, ncol = 1)
##' Z <- mcycle$accel
##' nvar <- 1
## ' plot(X, Z, ylim = c(-160, 90), ylab = 'acceleration', xlab = "time")
##'
##' ## Model fitting
##' model <- mleHomGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(10, nvar),
##'                   covtype = "Matern5_2", known = list(beta0 = 0))
##' LOO_p <- LOO_preds(model)
##'  
##' # model minus observation(s) at x_i
##' d_mot <- find_reps(X, Z)
##' 
##' LOO_ref <- matrix(NA, nrow(d_mot$X0), 2)
##' for(i in 1:nrow(d_mot$X0)){
##'  model_i <- mleHomGP(X = list(X0 = d_mot$X0[-i,, drop = FALSE], Z0 = d_mot$Z0[-i],
##'                      mult = d_mot$mult[-i]), Z = unlist(d_mot$Zlist[-i]),
##'                      lower = rep(0.1, nvar), upper = rep(50, nvar), covtype = "Matern5_2",
##'                      known = list(theta = model$theta, k_theta_g = model$k_theta_g, g = model$g,
##'                                   beta0 = 0))
##'  model_i$nu_hat <- model$nu_hat
## ' # For hetGP, need to use the same Lambdas to get the same results  
## '  model_i$Lambda <- model$Lambda[-i] 
## '  model_i <- strip(model_i)
## '  model_i <- rebuild(model_i)
##'  p_i <- predict(model_i, d_mot$X0[i,,drop = FALSE])
##'  LOO_ref[i,] <- c(p_i$mean, p_i$sd2)
##' }
##' 
##' # Compare results
##' 
##' range(LOO_ref[,1] - LOO_p$mean)
##' range(LOO_ref[,2] - LOO_p$sd2)
##' 
##' # Use of LOO for diagnostics
##' plot(model)
LOO_preds <- function(model, ids = NULL){
  if(is.null(ids)) ids <- 1:nrow(model$X0)
  
  if(!is.null(model$trendtype) && model$trendtype == "OK") 
    model$Ki <- model$Ki - matrix(rowSums(model$Ki), ncol = 1) %*% matrix(rowSums(model$Ki), nrow = 1) / sum(model$Ki)
  
  if(class(model) == "homGP"){
    sds <- model$nu_hat * (1/diag(model$Ki)[ids] - model$g/model$mult[ids])
  } 
  
  if(class(model) == "hetGP"){
    sds <- model$nu_hat * (1/diag(model$Ki)[ids] - model$Lambda[ids]/model$mult[ids])
  }
  
  if(class(model) == "homTP"){
    sds <- (1/diag(model$Ki)[ids] - model$g/model$mult[ids])
    # TP correction
    sds <- (model$nu + model$psi - 2) / (model$nu + length(model$Z) - model$mult[ids] - 2) * sds
  }
  
  if(class(model) == "hetTP"){
    sds <- (1/diag(model$Ki)[ids] - model$Lambda[ids]/model$mult[ids])
    # TP correction
    sds <- (model$nu + model$psi - 2) / (model$nu + length(model$Z) - model$mult[ids] - 2) * sds
  }
  
  ys <- model$Z0[ids] - (model$Ki %*% (model$Z0 - model$beta0))[ids]/diag(model$Ki)[ids]
  return(list(mean = ys, sd2 = sds))
}



