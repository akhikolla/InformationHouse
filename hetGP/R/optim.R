#' Computes EI for minimization
#' @title Expected Improvement criterion
#' @param x matrix of new designs, one point per row (size n x d)
#' @param model \code{homGP} or \code{hetGP} model, or their TP equivalents, including inverse matrices
#' @param cst optional plugin value used in the EI, see details
#' @param preds optional predictions at \code{x} to avoid recomputing if already done
#' @note This is a beta version at this point.
#' @details 
#' \code{cst} is classically the observed minimum in the determininistic case. 
#' In the noisy case, the min of the predictive mean works fine.
#'  
#' @references
#' Mockus, J.; Tiesis, V. & Zilinskas, A. (1978).
#' The application of Bayesian methods for seeking the extremum Towards Global Optimization, Amsterdam: Elsevier, 2, 2.\cr \cr
#' 
#' Vazquez E, Villemonteix J, Sidorkiewicz M, Walter E (2008). 
#' Global Optimization Based on Noisy Evaluations: An Empirical Study of Two Statistical Approaches, 
#' Journal of Physics: Conference Series, 135, IOP Publishing.\cr \cr
#' 
#' A. Shah, A. Wilson, Z. Ghahramani (2014), Student-t processes as alternatives to Gaussian processes, Artificial Intelligence and Statistics, 877--885.
#' @importFrom stats dt qt
#' @export
#' @examples 
#' ## Optimization example
#' set.seed(42)
#' 
#' 
#' ## Noise field via standard deviation
#' noiseFun <- function(x, coef = 1.1, scale = 1){
#' if(is.null(nrow(x)))
#'  x <- matrix(x, nrow = 1)
#'    return(scale*(coef + cos(x * 2 * pi)))
#' }
#' 
#' ## Test function defined in [0,1]
#' ftest <- function(x){
#' if(is.null(nrow(x)))
#' x <- matrix(x, ncol = 1)
#' return(f1d(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
#' }
#' 
#' n_init <- 10 # number of unique designs
#' N_init <- 100 # total number of points
#' X <- seq(0, 1, length.out = 10)
#' X <- matrix(X[sample(1:n_init, N_init, replace = TRUE)], ncol = 1)
#' Z <- ftest(X)
#' 
#' ## Predictive grid
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- matrix(xgrid, ncol = 1)
#' 
#' model <- mleHetGP(X = X, Z = Z, lower = 0.001, upper = 1)
#' 
#' EIgrid <- crit_EI(Xgrid, model)
#' preds <- predict(x = Xgrid, model)
#' 
#' par(mar = c(3,3,2,3)+0.1)
#' plot(xgrid, f1d(xgrid), type = 'l', lwd = 1, col = "blue", lty = 3,
#' xlab = '', ylab = '', ylim = c(-8,16))
#' points(X, Z)
#' lines(Xgrid, preds$mean, col = 'red', lwd = 2)
#' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
#' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
#' par(new = TRUE)
#' plot(NA, NA, xlim = c(0, 1), ylim = c(0,max(EIgrid)), axes = FALSE, ylab = "", xlab = "")
#' lines(xgrid, EIgrid, lwd = 2, col = 'cyan')
#' axis(side = 4)
#' mtext(side = 4, line = 2, expression(EI(x)), cex = 0.8)
#' mtext(side = 2, line = 2, expression(f(x)), cex = 0.8)
crit_EI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X0)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  if(class(model) %in% c("homTP", "hetTP")){
    gamma <- (cst - preds$mean)/sqrt(preds$sd2)
    res <- (cst - preds$mean) * pt(gamma, df = model$nu + length(model$Z))
    res <- res + sqrt(preds$sd2) * (1 + (gamma^2 - 1)/(model$nu + length(model$Z) - 1)) * dt(x = gamma, df = model$nu + length(model$Z))
    res[which(res < 1e-12)] <- 0 # for stability
    return(res)
  }
  
  xcr <- (cst - preds$mean)/sqrt(preds$sd2)
  res <- (cst - preds$mean) * pnorm(xcr)
  res <- res + sqrt(preds$sd2) * dnorm(xcr)
  res[which(preds$sd2 < sqrt(.Machine$double.eps) | res < 1e-12)] <- 0 # for stability
  return(res) 
}

#' Derivative of EI criterion for GP models
#' @param x matrix for the new design (size 1 x d)
#' @param model \code{homGP} or \code{hetGP} model
#' @param cst threshold for contour criteria
#' @param preds pre-computed preds for contour criteria
#' @export
#' @seealso \code{\link[hetGP]{crit_EI}} for the criterion
#' @references 
#' Ginsbourger, D. Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables Ecole Nationale Superieure des Mines de Saint-Etienne, Ecole Nationale Superieure des Mines de Saint-Etienne, 2009 \cr \cr
#' Roustant, O., Ginsbourger, D., DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization, Journal of Statistical Software, 2012
deriv_crit_EI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  pred_gr <- predict_gr(model, x)
  
  z <- (cst - preds$mean)/sqrt(preds$sd2)
  
  if(class(model) %in% c("homGP", "hetGP")){
    res <- pred_gr$sd2 / (2 * sqrt(preds$sd2)) * dnorm(z) - pred_gr$mean * pnorm(z)
  }else{
    # dz = - dm/s - z ds/s = -dm/s - z * ds2/(2s2) 
    dz <- -pred_gr$mean/sqrt(preds$sd2) - z * pred_gr$sd2/(2 * matrix(preds$sd2, nrow = nrow(x), ncol(x)))
    
    # d( (cst - m(x)).pt(z(x))) = 
    p1 <- -pred_gr$mean * pt(z, df = model$nu + length(model$Z)) + (cst - preds$mean) * dz * dt(z, df = model$nu + length(model$Z))
    
    a <- model$nu + length(model$Z) - 1
    # d( s(x) (1 + (z^2-1)/(nu + N -1)) dt(z(x)) (in 2 lines)
    p2 <- (pred_gr$sd2/(2*sqrt(preds$sd2)) * (1 + (z^2 - 1)/a) + 2*sqrt(preds$sd2) * z * dz/a) * dt(z, df = model$nu + length(model$Z))
    p2 <- p2 + sqrt(preds$sd2) * (1 + (z^2 - 1)/a) * dz * dlambda(z, model$nu + length(model$Z))
    res <- p1 + p2
  }
  
  res[which(abs(res) < 1e-12)] <- 0 # for stability with optim
  
  return(res)
}

#' Prediction of the gradient
#' @param object GP/TP
#' @param x design location
#' @noRd
predict_gr <- function(object, x){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  kvec <-  cov_gen(X1 = object$X0, X2 = x, theta = object$theta, type = object$covtype)
  
  dm <- ds2 <- matrix(NA, nrow(x), ncol(x))
  
  for(i in 1:nrow(x)){
    dkvec <- matrix(NA, nrow(object$X0), ncol(x))
    for(j in 1:ncol(x)) dkvec[, j] <- drop(partial_cov_gen(X1 = x[i,,drop = F], X2 = object$X0, theta = object$theta, i1 = 1, i2 = j, arg = "X_i_j", type = object$covtype)) * kvec[,i]
    
    dm[i,] <- crossprod(object$Z0 - object$beta0, object$Ki) %*% dkvec
    if(class(object) %in% c("hetGP", "homGP") && object$trendtype == "OK") tmp <- drop(1 - colSums(object$Ki) %*% kvec[,i])/(sum(object$Ki)) * colSums(object$Ki) %*% dkvec else tmp <- 0
    ds2[i,] <- -2 * (crossprod(kvec[,i], object$Ki) %*% dkvec + tmp)
  }
  if(class(object) %in% c("hetGP", "homGP")){
    return(list(mean = dm, sd2 = object$nu_hat * ds2))
  }else{
    return(list(mean = object$sigma2 * dm, sd2 =  (object$nu + object$psi - 2) / (object$nu + length(object$Z) - 2) * object$sigma2^2 * ds2)) 
  }
}

#' Derivative of the student-t pdf
#' @param z input location
#' @param a degree of freedom parameter
#' @noRd
#' @example 
#' # grad(dt, x = 0.55, df = 3.6)
#' # dlambda(0.55, 3.6)
dlambda <- function(z, a){
  return(-(a + 1) * gamma((a + 1)/2)/(sqrt(pi * a) * a * gamma(a/2)) * z * ((a + z^2)/a)^(-(a +3)/2))
}

#' Search for best reduction in a criterion 
#' @title Criterion minimization
#' @param model \code{homGP} or \code{hetGP} model
#' @param crit considered criterion, one of \code{"crit_cSUR"}, \code{"crit_EI"}, \code{"crit_ICU"},
#'  \code{"crit_MCU"} and \code{"crit_tMSE"}. Note that \code{crit_IMSPE} has its dedicated method, see \code{\link[hetGP]{IMSPE_optim}}.
#' @param ... additional parameters of the criterion
#' @param replicate if \code{TRUE}, search only on existing designs
#' @param Xcand optional set of of candidates for discrete search
#' @param control list in case \code{Xcand == NULL}, with elements \code{multi.start},
#' to perform a multi-start optimization based on \code{\link[stats]{optim}}, with \code{maxit} iterations each.
#' Also, \code{tol_dist} defines the minimum distance to an existing design for a new point to be added, otherwise the closest existing design is chosen.
#' In a similar fashion, \code{tol_dist} is the minimum relative change of crit for adding a new design.
#' @param seed optional seed for the generation of designs with \code{\link[DiceDesign]{maximinSA_LHS}}
#' @param ncores number of CPU available (> 1 mean parallel TRUE), see \code{\link[parallel]{mclapply}}
#' @importFrom DiceDesign lhsDesign maximinSA_LHS
#' @return list with \code{par}, \code{value} elements, and additional slot \code{new} (boolean if it is or not a new design) and \code{id} giving the index of the duplicated design. 
#' @noRd
#' @examples 
#' ###############################################################################
#' ## Bi-variate example
#' ###############################################################################
#' 
#' nvar <- 2 
#' 
#' set.seed(42)
#' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
#' 
#' n <- 25 # must be a square
#' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
#' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
#' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
#' Z <- apply(X, 1, ftest)
#' 
#' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
#' 
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' 
#' preds <- predict(x = Xgrid, object =  model)
#' 
#' ## Initial plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit <- "crit_EI"
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Initial criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
#' ## Sequential crit search
#' nsteps <- 1 # Increase for better results
#' 
#' for(i in 1:nsteps){
#'   res <- crit.search(model, crit = crit, control = list(multi.start = 100, maxit = 50))
#'   newX <- res$par
#'   newZ <- ftest(newX)
#'   model <- update(object = model, Xnew = newX, Znew = newZ)
#' }
#' 
#' ## Final plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Final criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
crit.search <- function(model, crit, ..., replicate = FALSE, Xcand = NULL, 
                        control = list(tol_dist = 1e-6, tol_diff = 1e-6, multi.start = 20,
                                       maxit = 100, maximin = TRUE, Xstart = NULL), seed = NULL,
                        ncores = 1){
  # Only search on existing designs
  if(replicate){
    ## Discrete optimization
    res <- mclapply(1:nrow(model$X0), 
                    function(i) match.fun(crit)(x = model$X0[i,,drop = FALSE], model = model, ... = ...), mc.cores = ncores)
    res <- unlist(res)
    
    return(list(par = model$X0[which.max(res),,drop = FALSE], value = max(res), new = FALSE, id = which.max(res)))
  }
  
  if(is.null(control))
    control <- list(multi.start = 20, maxit = 100)
  
  if(is.null(control$multi.start))
    control$multi.start <- 20
  
  if(is.null(control$maxit))
    control$maxit <- 100
  
  if(is.null(control$maximin))
    control$maximin <- TRUE
  
  if(is.null(control$tol_dist)) control$tol_dist <- 1e-6
  if(is.null(control$tol_diff)) control$tol_diff <- 1e-6
  
  d <- ncol(model$X0)
  
  if(crit == "crit_EI") gr <- deriv_crit_EI else gr <- NULL
  crit <- match.fun(crit)
  
  ## Optimization
  if(is.null(Xcand)){
    ## Continuous optimization
    if(!is.null(control$Xstart)){
      Xstart <- control$Xstart
    }else{
      if(is.null(seed)) seed <- sample(1:2^15, 1) ## To be changed?
      if(control$maximin){
        if(d == 1){
          # perturbed 1d equidistant points
          Xstart <- matrix(seq(1/2, control$multi.start -1/2, length.out = control$multi.start) + runif(control$multi.start, min = -1/2, max = 1/2), ncol = 1)/control$multi.start
        }else{
          Xstart <- maximinSA_LHS(lhsDesign(control$multi.start, d, seed = seed)$design)$design
        }
      }else{
        Xstart <- lhsDesign(control$multi.start, d, seed = seed)$design
      }
    }
    
    res <- list(par = NA, value = -Inf, new = NA)
    
    local_opt_fun <- function(i){
      out <- try(optim(Xstart[i,, drop = FALSE], crit, ... = ..., method = "L-BFGS-B", gr = gr, 
                       lower = rep(0, d), upper = rep(1, d),
                       model = model, control = list(maxit = control$maxit, fnscale = -1)))
      if(class(out) == "try-error") return(NULL)
      return(out)
    }
    
    all_res <- mclapply(1:nrow(Xstart), local_opt_fun, mc.cores = ncores)
    res_max <- which.max(Reduce(c, lapply(all_res, function(x) x$value)))
    res <- list(par = apply(all_res[[res_max]]$par, c(1, 2), function(x) max(min(x , 1), 0)),
                value = all_res[[res_max]]$value, new = TRUE, id = NULL)
    
    # for(i in 1:nrow(Xstart)){
    #   out <- try(optim(Xstart[i,, drop = FALSE], crit, ... = ..., method = "L-BFGS-B", gr = gr, 
    #                    lower = rep(0, d), upper = rep(1, d),
    #                    model = model, control = list(maxit = control$maxit, fnscale = -1)))
    #   if(class(out) != "try-error"){
    #     if(out$value > res$value)
    #       res <- list(par = out$par, value = out$value, new = TRUE, id = NULL)
    #   }
    # }
    
    if(control$tol_dist > 0 || control$tol_diff > 0){
      ## Check if new design is not to close to existing design
      dists <- sqrt(distance_cpp(res$par, model$X0))
      if(min(dists) < control$tol_dist){
        res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                    value = crit(x = model$X0[which.min(dists),, drop = F], model = model, ... = ...),
                    new = FALSE, id = which.min(dists))
      }else{
        ## Check if crit difference between replication and new design is significative
        id_closest <- which.min(dists) # closest point to new design
        crit_rep <- crit(model$X0[which.min(dists),,drop = FALSE], model = model, ...=...)
        
        ## EI can be 0, in which case it is better to replicate even if crit_rep$value is also 0
        if(res$value == 0 || (res$value - crit_rep)/res$value < control$tol_diff){
          res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                      value = crit_rep,
                      new = FALSE, id = which.min(dists))
        }
      }
    }
    
    
    return(res)
    
    
  }else{
    ## Discrete optimization
    res <- mclapply(1:nrow(model$X0), 
                    function(i) match.fun(crit)(x = Xcand[i,,drop = FALSE], model = model, ... = ...), mc.cores = ncores)
    res <- unlist(res)
    
    tmp <- which(duplicated(rbind(model$X0, Xcand[which.max(res),,drop = FALSE]), fromLast = TRUE))
    if(length(tmp) > 0) return(list(par = Xcand[which.max(res),,drop = FALSE], value = max(res), new = FALSE, id = tmp))
    return(list(par = Xcand[which.max(res),,drop = FALSE], value = max(res), new = TRUE, id = NULL))
  }
}


#' Search for the best value of available criterion, possibly using a h-steps lookahead strategy to favor designs with replication
#' @title Criterion optimization
#' @param model \code{homGP} or \code{hetGP} model
#' @param crit considered criterion, one of \code{"crit_cSUR"}, \code{"crit_EI"}, \code{"crit_ICU"},
#'  \code{"crit_MCU"} and \code{"crit_tMSE"}. Note that \code{crit_IMSPE} has its dedicated method, see \code{\link[hetGP]{IMSPE_optim}}.
#' @param ... additional parameters of the criterion
##' @param Xcand optional discrete set of candidates (otherwise a maximin LHS is used to initialise continuous search)
#' @param control list in case \code{Xcand == NULL}, with elements \code{multi.start},
#' to perform a multi-start optimization based on \code{\link[stats]{optim}}, with \code{maxit} iterations each.
#' Also, \code{tol_dist} defines the minimum distance to an existing design for a new point to be added, otherwise the closest existing design is chosen.
#' In a similar fashion, \code{tol_dist} is the minimum relative change of crit for adding a new design.
#' @param seed optional seed for the generation of LHS designs with \code{\link[DiceDesign]{maximinSA_LHS}}
#' @param h horizon for multi-step ahead framework.
#' The decision is made between:
#' \itemize{
#'  \item sequential crit search starting by a new design (optimized first) then adding \code{h} replicates
#'  \item sequential crit searches starting by \code{1} to \code{h} replicates before adding a new point
#' }
#' Use \code{h = 0} for the myopic criterion, i.e., not looking ahead.
#' @param ncores number of CPU available (> 1 mean parallel TRUE), see \code{\link[parallel]{mclapply}}
#' @details 
#' When looking ahead, the kriging believer heuristic is used,
#'  meaning that the non-observed value is replaced by the mean prediction in the update.  
#' @return list with elements:
#' \itemize{
#' \item \code{par}: best first design,
#' \item \code{value}: criterion h-steps ahead starting from adding \code{par},
#' \item \code{path}: list of elements list(\code{par}, \code{value}, \code{new}) at each step \code{h}
#' }
#' @references
#' M. Binois, J. Huang, R. B. Gramacy, M. Ludkovski (2018+), Replication or exploration? Sequential design for stochastic simulation experiments,
#' Technometrics (to appear).\cr 
#' Preprint available on arXiv:1710.03206.
#' @export 
#' @examples
#' ###############################################################################
#' ## Bi-variate example (myopic version)
#' ###############################################################################
#' 
#' nvar <- 2 
#' 
#' set.seed(42)
#' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
#' 
#' n <- 25 # must be a square
#' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
#' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
#' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
#' Z <- apply(X, 1, ftest)
#' 
#' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
#' 
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' 
#' preds <- predict(x = Xgrid, object =  model)
#' 
#' ## Initial plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit <- "crit_EI"
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Initial criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
#' ## Sequential crit search
#' nsteps <- 1 # Increase for better results
#' 
#' for(i in 1:nsteps){
#'   res <- crit_optim(model, crit = crit, h = 0, control = list(multi.start = 50, maxit = 30))
#'   newX <- res$par
#'   newZ <- ftest(newX)
#'   model <- update(object = model, Xnew = newX, Znew = newZ)
#' }
#' 
#' ## Final plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Final criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
#' ###############################################################################
#' ## Bi-variate example (look-ahead version)
#' ###############################################################################
#' \dontrun{  
#' nvar <- 2 
#' 
#' set.seed(42)
#' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
#' 
#' n <- 25 # must be a square
#' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
#' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
#' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
#' Z <- apply(X, 1, ftest)
#' 
#' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
#' 
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' 
#' nsteps <- 5 # Increase for more steps
#' crit <- "crit_EI"
#' 
#' # To use parallel computation (turn off on Windows)
#' library(parallel)
#' parallel <- FALSE #TRUE #
#' if(parallel) ncores <- detectCores() else ncores <- 1
#' 
#' for(i in 1:nsteps){
#'   res <- crit_optim(model, h = 3, crit = crit, ncores = ncores,
#'                     control = list(multi.start = 100, maxit = 50))
#'   
#'   # If a replicate is selected
#'   if(!res$path[[1]]$new) print("Add replicate")
#'   
#'   newX <- res$par
#'   newZ <- ftest(newX)
#'   model <- update(object = model, Xnew = newX, Znew = newZ)
#'   
#'   ## Plots 
#'   preds <- predict(x = Xgrid, object =  model)
#'   contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'           main = "Predicted mean", nlevels = 20)
#'   points(model$X0, col = 'blue', pch = 20)
#'   points(newX, col = "red", pch = 20)
#'   
#'   crit_grid <- apply(Xgrid, 1, crit, model = model)
#'   filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                  nlevels = 20, color.palette = terrain.colors,
#'   plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' }
#' }
crit_optim <- function(model, crit, ..., h = 2, Xcand = NULL, control = list(multi.start = 10, maxit = 100), seed = NULL, ncores = 1){
  d <- ncol(model$X0)
  
  if(crit == "crit_IMSPE") stop("crit_IMSPE is intended to be optimized by IMSPE_optim")
  
  ## A) Setting to beat: first new point then replicate h times
  crit_A <- crit.search(model = model, crit = crit, ... = ..., control = control, Xcand = Xcand, seed = seed, ncores = ncores)
  new_designA <- crit_A$par ## store first considered design to be added
  path_A <- list(crit_A)
  
  if(h > 0){
    newmodelA <- model
    for(i in 1:h){
      ZnewA <- predict(newmodelA, crit_A$par)$mean
      newmodelA <- update(object = newmodelA, Xnew = crit_A$par, Znew = ZnewA, maxit = 0)
      crit_A <- crit.search(model = newmodelA, crit = crit, ... = ..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
      path_A <- c(path_A, list(crit_A))
    }
  }
  
  if(h == -1) return(list(par = new_designA, value = crit_A$value, path = path_A)) 
  
  ## B) Now compare with waiting to add new point
  newmodelB <- model
  
  if(h == 0){
    crit_B <- crit.search(model = newmodelB, crit = crit, ... =..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
    new_designB <- crit_B$par ## store considered design to be added
    
    # search from best replicate
    if(is.null(Xcand)){
      crit_C <- crit.search(model = newmodelB, crit = crit,
                            control = list(Xstart = crit_B$par, maxit = control$maxit,
                                           tol_dist = control$tol_dist, tol_diff = control$tol_diff), seed = seed, ncores = ncores)
    }else{
      crit_C <- crit_B
    }
    
    if(crit_C$value > max(crit_A$value, crit_B$value)) return(list(par = crit_C$par, value = crit_C$value, path = list(crit_C)))
    
    if(crit_B$value > crit_A$value){
      return(list(par = crit_B$par, value = crit_B$value, path = list(crit_B)))
    } 
  }else{
    for(i in 1:h){
      ## Add new replicate
      crit_B <- crit.search(model = newmodelB, crit = crit, ... = ..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
      
      if(i == 1){
        new_designB <- matrix(crit_B$par, nrow = 1) ##store first considered design to add
        path_B <- list()
      } 
      
      path_B <- c(path_B, list(crit_B))
      ZnewB <- predict(newmodelB, crit_B$par)$mean
      newmodelB <- update(object = newmodelB, Xnew = crit_B$par, Znew = ZnewB, maxit = 0)
      
      ## Add new design
      crit_C <- crit.search(model = newmodelB, crit = crit, ... = ..., control = control, Xcand = Xcand, seed = seed, ncores = ncores)
      path_C <- list(crit_C)
      
      if(i < h){
        newmodelC <- newmodelB
        
        for(j in i:(h-1)){
          ## Add remaining replicates
          ZnewC <- predict(newmodelC, crit_C$par)$mean
          newmodelC <- update(object = newmodelC, Xnew = crit_C$par, Znew = ZnewC, maxit = 0)
          crit_C <- crit.search(model = newmodelC, crit = crit, ... = ..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
          path_C <- c(path_C, list(crit_C))
        }
      }
      
      if(crit_C$value < crit_A$value) return(list(par = new_designB, value = crit_C$value, path = c(path_B, path_C)))
    }
  }
  
  return(list(par = new_designA, value = crit_A$value, path = path_A))
  
}

