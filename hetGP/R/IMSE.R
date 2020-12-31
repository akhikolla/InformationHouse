##' IMSPE of a given design
##' @title Integrated Mean Square Prediction Error
##' @param X \code{hetGP} or \code{homGP} model. Alternatively, one can provide a matrix of unique designs considered
##' @param theta lengthscales
##' @param Lambda diagonal matrix for the noise
##' @param mult number of replicates at each design
##' @param covtype either "Gaussian", "Matern3_2" or "Matern5_2"
##' @param nu variance parameter
##' @param eps numerical nugget
##' @details
##' One can provide directly a model of class \code{hetGP} or \code{homGP}, or provide \code{X} and all other arguments
##' @export
IMSPE <- function(X, theta = NULL, Lambda = NULL, mult = NULL, covtype = NULL, nu= NULL, eps = sqrt(.Machine$double.eps)){
  if(class(X) %in% c("homGP", "hetGP")){
    Wij <- Wij(mu1 = X$X0, theta = X$theta, type = X$covtype)
    if(X$trendtype == "OK"){
      tmp <- drop(1 - 2 * colSums(X$Ki) %*% mi(mu1 = X$X0, theta = X$theta, type = X$covtype) + colSums(X$Ki) %*% Wij %*% rowSums(X$Ki))/sum(X$Ki)
    }else 
      tmp <- 0
    return(X$nu_hat * (1 - sum(X$Ki * Wij) + tmp))
  }else{
    return(nu*(1 - sum(chol2inv(chol(add_diag(cov_gen(X, theta = theta, type = covtype), Lambda/mult + eps))) * Wij(mu1 = X, theta = theta, type = covtype))))
  }
}

##' Compute the integrated mean square prediction error after adding a new design
##' @title Sequential IMSPE criterion
##' @param x matrix for the new design (size 1 x d)
##' @param Wijs optional previously computed matrix of Wijs, to avoid recomputing it; see \code{\link[hetGP]{Wij}}
##' @param model \code{homGP} or \code{hetGP} model, including inverse matrices
##' @param id instead of providing \code{x}, one can provide the index of a considered existing design 
##' @export
##' @details The computations are scale free, i.e., values are not multiplied by \code{nu_hat} from \code{homGP} or \code{hetGP}.
##' Currently this function ignores the extra terms related to the estimation of the mean.
##' @seealso \code{\link[hetGP]{deriv_crit_IMSPE}} for the derivative
##' @examples
##' ## One-d toy example
##' 
##' set.seed(42)
##' ftest <- function(x, coef = 0.1) return(sin(2*pi*x) + rnorm(1, sd = coef))
##' 
##' n <- 9
##' designs <- matrix(seq(0.1, 0.9, length.out = n), ncol = 1)
##' X <- matrix(designs[rep(1:n, sample(1:10, n, replace = TRUE)),])
##' Z <- apply(X, 1, ftest)
##' 
##' prdata <- find_reps(X, Z, inputBounds = matrix(c(0,1), nrow = 2, ncol = 1))
##' Z <- prdata$Z
##' plot(prdata$X0[rep(1:n, times = prdata$mult),], prdata$Z, xlab = "x", ylab = "Y")
##' 
##' model <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult),
##'                   Z = Z, lower = 0.1, upper = 5)
##' 
##' ngrid <- 501
##' xgrid <- matrix(seq(0,1, length.out = ngrid), ncol = 1)
##' 
##' ## Precalculations
##' Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
##' 
## ' nref <- 2000
## ' # library(DiceDesign)
## ' # Xref <- maximinSA_LHS(lhsDesign(nref, nvar)$design)$design
##' 
##' t0 <- Sys.time()
##' 
##' IMSPE_grid <- apply(xgrid, 1, crit_IMSPE, Wijs = Wijs, model = model)
##' 
##' t1 <- Sys.time()
##' print(t1 - t0)
##' 
## ' ## Older version with candidate points (Warning: compare for fixed zero mean)
## ' Xref <- xgrid
## ' # deprecated, for testing only
## ' IMSPE_disc_grid <- apply(xgrid, 1, hetGP:::IMSPE_grid, model = model, Xref = Xref) 
## ' 
## ' t2 <- Sys.time()
## ' 
## ' print(t2-t1)
## ' 
## ' par(mfrow = c(1,2))
##' plot(xgrid, IMSPE_grid * model$nu_hat, xlab = "x", ylab = "crit_IMSPE values")
##' abline(v = designs)
## ' plot(xgrid, IMSPE_disc_grid / ngrid)
## ' abline(v = designs)
## ' par(mfrow = c(1,1))
##' 
##' ###############################################################################
##' ## Bi-variate case
##' 
##' nvar <- 2 
##' 
##' set.seed(2)
##' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
##' 
##' n <- 16 # must be a square
##' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
##' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
##' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
##' Z <- apply(X, 1, ftest)
##' 
##' prdata <- find_reps(X, Z, inputBounds = matrix(c(0,1), nrow = 2, ncol = 1))
##' Z <- prdata$Z
##' 
##' model <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = Z, 
##'  lower = rep(0.1, nvar), upper = rep(1, nvar))
##' ngrid <- 51
##' xgrid <- seq(0,1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' ## Precalculations
##' Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
##' t0 <- Sys.time()
##' 
##' IMSPE_grid <- apply(Xgrid, 1, crit_IMSPE, Wijs = Wijs, model = model)
##' filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid),
##'                nlevels = 20, color.palette = terrain.colors,
##'                main = "Sequential IMSPE values")
## ' filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid), nlevels = 20, color.palette = terrain.colors)
crit_IMSPE <- function(x, model, id = NULL, Wijs = NULL){

  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)

  if(!is.null(id)){
    if(is.null(model$Lambda)){
      tmp <- model$g
    }else{
      tmp <- model$Lambda[id]
    }
    return(1 - (sum(model$Ki*Wijs) + crossprod(model$Ki[,id], Wijs %*% model$Ki[,id]) / ((model$mult[id]*(model$mult[id] + 1)) / (1*tmp) - model$Ki[id, id])))
  }

  if(is.null(dim(x)))
    x <- matrix(x, 1)

  newWijs <- Wij(mu2 = x, mu1 = model$X0, theta = model$theta, type = model$covtype)
  W11 <- Wij(mu2 = x, mu1 = x, theta = model$theta, type = model$covtype)

  kn1 <- cov_gen(x, model$X0, theta = model$theta, type = model$covtype)
  
  new_lambda <- predict(object = model, x = x, nugs.only = TRUE)$nugs/model$nu_hat
  
  vn <- drop(1 - kn1 %*% tcrossprod(model$Ki, kn1)) + new_lambda + model$eps
  gn <- - tcrossprod(model$Ki, kn1) / vn
  
  return(1 - (sum(model$Ki*Wijs) + crossprod(gn, Wijs %*% gn)*vn + 2 * crossprod(newWijs, gn) + W11 / vn))
}


## Wrapper functions

d1 <- function(X, x, sigma, type){
  if(type == "Gaussian"){
    return(d_gauss_cpp(X = X, x = x, sigma = sigma))
  }
  if(type == "Matern5_2"){
    return(d_mat52_cpp(X = X, x = x, sigma = sigma))
  }
  if(type == "Matern3_2"){
    return(d_mat32_cpp(X = X, x = x, sigma = sigma))
  }
}

c1 <- function(X, x, sigma, W, type){
  if(type == "Gaussian"){
    return(c1_gauss_cpp(X = X, x = x, sigma = sqrt(sigma), W = W))
  }
  if(type == "Matern5_2"){
    return(c1_mat52_cpp(X = X, x = x, sigma = sigma, W = W))
  }
  if(type == "Matern3_2"){
    return(c1_mat32_cpp(X = X, x = x, sigma = sigma, W = W))
  }
}
  
c2 <- function(x, sigma, w, type){
  if(type == "Gaussian"){
    return(c2_gauss_cpp(x = x, t = sqrt(sigma), w = w))
  }
  if(type == "Matern5_2"){
    return(c2_mat52_cpp(x = x, t = sigma, w = w))
  }
  if(type == "Matern3_2"){
    return(c2_mat32_cpp(x = x, t = sigma, w = w))
  }
}

##' Derivative of crit_IMSPE
##' @param x matrix for the new design (size 1 x d)
##' @param model \code{homGP} or \code{hetGP} model
##' @param Wijs optional previously computed matrix of Wijs, see \code{\link[hetGP]{Wij}}
##' @return Derivative of the sequential IMSPE with respect to \code{x}
##' @seealso \code{\link[hetGP]{crit_IMSPE}} for the criterion
##' @export
## ' @examples
## ' \dontrun{
## ' nvar <- 2
## ' set.seed(42)
## ' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
## ' n <- nvar * 10
## ' designs <- matrix(runif(n*nvar), n)
## ' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
## ' Z <- apply(X, 1, ftest)
## ' 
## ' covtype <- "Gaussian"
## ' model <- mleHetGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(1, nvar), covtype = covtype)
## ' 
## ' library(numDeriv)
## ' xnew <- matrix(runif(nvar), nrow = 1)
## ' 
## ' print(grad(crit_IMSPE, x = xnew, model = model))
## ' deriv_crit_IMSPE(xnew, model) 
## ' 
## ' }
deriv_crit_IMSPE <- function(x, model, Wijs = NULL){
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
  
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  
  kn1 <- cov_gen(model$X0, x, theta = model$theta, type = model$covtype)
  if(class(model) == 'hetGP') kng1 <- cov_gen(model$X0, x, theta = model$theta_g, type = model$covtype)
  new_lambda <- predict(object = model, x = x, nugs.only = TRUE)$nugs/model$nu_hat
  k11 <- 1 + new_lambda
  
  W1 <- Wij(mu1 = model$X, mu2 = x, theta = model$theta, type = model$covtype)
  W11 <- Wij(mu1 = x, theta = model$theta, type = model$covtype)
  
  # compute derivative vectors and scalars
  Kikn1 <-  model$Ki %*% kn1
  v <- drop(k11 - crossprod(kn1, Kikn1))
  g <-  - Kikn1 / v
  
  tmp <- rep(NA, ncol(x))
  dlambda <- 0
  if(class(model) == 'hetGP') KgiD <- model$Kgi %*% (model$Delta - model$nmean)
  if(length(model$theta) < length(x)) model$theta <- rep(model$theta, length(x))
  if(class(model) == 'hetGP' && length(model$theta_g) < length(x)) model$theta_g <- rep(model$theta, length(x))
  
  Wig <- Wijs %*% g
  
  for (m in 1:length(x)){
    
    c1_v <- c1(X = model$X0[,m], x = x[m], sigma = model$theta[m], W = W1, type = model$covtype)
    
    dis <- drop(d1(X = model$X0[,m], x = x[m], sigma = model$theta[m], type = model$covtype) * kn1)
    
    if(class(model)=='hetGP'){
      dlambda <- crossprod(d1(X = model$X0[,m], x = x[m], sigma = model$theta_g[m], type = model$covtype) * kng1, KgiD) 
      if(model$logN)
        dlambda <- new_lambda *dlambda
    }
    
    v2 <- drop(2 * (- dis %*% Kikn1) + dlambda) 
    v1 <- -1/v^2 * v2
    h <- drop(-model$Ki %*% (v1 * kn1 + 1/v * dis))
    
    tmp[m] <- 2 * crossprod(c1_v, g) + c2(x = x[m], sigma = model$theta[m], w = W11, type = model$covtype) / v + crossprod(v2 * g + 2 * v * h, Wig) + 2 * crossprod(h, W1) + v1 * W11
  }
  return(-tmp)
}

##' Search for best reduction in IMSPE 
##' @title IMSPE minimization
##' @param model \code{homGP} or \code{hetGP} model
##' @param replicate if \code{TRUE}, search only on existing designs
##' @param Xcand optional set of of candidates for discrete search
##' @param control list in case \code{Xcand == NULL}, with elements \code{multi.start},
##' to perform a multi-start optimization based on \code{\link[stats]{optim}}, with \code{maxit} iterations each.
##' Also, \code{tol_dist} defines the minimum distance to an existing design for a new point to be added, otherwise the closest existing design is chosen.
##' In a similar fashion, \code{tol_dist} is the minimum relative change of IMSPE for adding a new design.
##' @param Wijs optional previously computed matrix of Wijs, see \code{\link[hetGP]{Wij}}
##' @param seed optional seed for the generation of designs with \code{\link[DiceDesign]{maximinSA_LHS}}
##' @param ncores number of CPU available (> 1 mean parallel TRUE), see \code{\link[parallel]{mclapply}}
##' @importFrom DiceDesign lhsDesign maximinSA_LHS
##' @importFrom parallel mclapply
##' @return list with \code{par}, \code{value} elements, and additional slot \code{new} (boolean if it is or not a new design) and \code{id} giving the index of the duplicated design. 
##' @noRd
##' @examples 
##' ###############################################################################
##' ## Bi-variate example
##' ###############################################################################
##' 
##' nvar <- 2 
##' 
##' set.seed(42)
##' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
##' 
##' n <- 25 # must be a square
##' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
##' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
##' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
##' Z <- apply(X, 1, ftest)
##' 
##' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
##' 
##' ngrid <- 51
##' xgrid <- seq(0,1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' 
##' preds <- predict(x = Xgrid, object =  model)
##' 
##' ## Initial plots
##' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
##'         main = "Predicted mean", nlevels = 20)
##' points(model$X0, col = 'blue', pch = 20)
##' 
##' IMSPE_grid <- apply(Xgrid, 1, crit_IMSPE, model = model)
##' filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid),
##'                nlevels = 20, color.palette = terrain.colors, 
##'                main = "Initial IMSPE criterion landscape",
##' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
##' 
##' ## Sequential IMSPE search
##' nsteps <- 1 # Increase for better results
##' 
##' for(i in 1:nsteps){
##'   res <- IMSPE.search(model, control = list(multi.start = 100, maxit = 50))
##'   newX <- res$par
##'   newZ <- ftest(newX)
##'   model <- update(object = model, Xnew = newX, Znew = newZ)
##' }
##' 
##' ## Final plots
##' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
##'         main = "Predicted mean", nlevels = 20)
##' points(model$X0, col = 'blue', pch = 20)
##' 
##' IMSPE_grid <- apply(Xgrid, 1, crit_IMSPE, model = model)
##' filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid),
##'                nlevels = 20, color.palette = terrain.colors, 
##'                main = "Final IMSPE criterion landscape",
##' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
##' 
IMSPE.search <- function(model, replicate = FALSE, Xcand = NULL, 
                        control = list(tol_dist = 1e-6, tol_diff = 1e-6, multi.start = 20,
                                       maxit = 100, maximin = TRUE, Xstart = NULL), Wijs = NULL, seed = NULL,
                        ncores = 1){
  # Only search on existing designs
  if(replicate){
    ## Discrete optimization
    res <- unlist(mclapply(1:nrow(model$X0), crit_IMSPE, Wijs = Wijs, model = model, x = NULL, mc.cores = ncores))
    
    return(list(par = model$X0[which.min(res),,drop = FALSE], value = min(res), new = FALSE, id = which.min(res)))
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
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
  

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
    
    res <- list(par = NA, value = Inf, new = NA)
    
    local_opt_fun <- function(i){
      out <- try(optim(Xstart[i,, drop = FALSE], crit_IMSPE, method = "L-BFGS-B", lower = rep(0, d), upper = rep(1, d),
                   Wijs = Wijs, model = model, control = list(maxit = control$maxit), gr = deriv_crit_IMSPE))
      if(class(out) == "try-error") return(NULL)
      return(out)
    }
    
    all_res <- mclapply(1:nrow(Xstart), local_opt_fun, mc.cores = ncores)
    res_min <- which.min(Reduce(c, lapply(all_res, function(x) x$value)))
    res <- list(par = apply(all_res[[res_min]]$par, c(1, 2), function(x) max(min(x , 1), 0)),
                value = all_res[[res_min]]$value, new = TRUE, id = NULL)
    
    # for(i in 1:nrow(Xstart)){
    #   out <- optim(Xstart[i,, drop = FALSE], crit_IMSPE, method = "L-BFGS-B", lower = rep(0, d), upper = rep(1, d),
    #                Wijs = Wijs, model = model, control = list(maxit = control$maxit), gr = deriv_crit_IMSPE)
    #   if(out$value < res$value){
    #     out$par <- apply(out$par, c(1,2), function(x) max(min(x , 1), 0))
    #     res <- list(par = out$par, value = out$value, new = TRUE, id = NULL)
    #   }
    # }
    
    if(control$tol_dist > 0 && control$tol_diff > 0){
      ## Check if new design is not to close to existing design
      dists <- sqrt(distance_cpp(res$par, model$X0))
      if(min(dists) < control$tol_dist){
        res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                    value = crit_IMSPE(x = model$X0[which.min(dists),, drop = F], model = model, id = which.min(dists), Wijs = Wijs),
                    new = FALSE, id = which.min(dists))
      }else{
        ## Check if IMSPE difference between replication and new design is significative
        id_closest <- which.min(dists) # closest point to new design
        imspe_rep <- crit_IMSPE(model = model, id = id_closest, Wijs = Wijs)
        if((imspe_rep - res$value)/res$value < control$tol_diff){ #(1 - sum(model$Ki * Wijs)) < ){
          res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                      value = imspe_rep,
                      new = FALSE, id = which.min(dists))
        }
      }
    }

    
    return(res)
    
    
  }else{
    ## Discrete optimization
    res <- unlist(mclapply(Xcand, 1, crit_IMSPE, Wijs = Wijs, model = model, ncores = ncores))
    
    tmp <- which(duplicated(rbind(model$X0, Xcand[which.min(res),,drop = FALSE]), fromLast = TRUE))
    if(length(tmp) > 0) return(list(par = Xcand[which.min(res),,drop = FALSE], value = min(res), new = FALSE, id = tmp))
    return(list(par = Xcand[which.min(res),,drop = FALSE], value = min(res), new = TRUE, id = NULL))
  }
    
  
}


#' Search for the best value of the IMSPE criterion, possibly using a h-steps lookahead strategy to favor designs with replication
##' @title IMSPE optimization
##' @param model \code{homGP} or \code{hetGP} model
##' @param Xcand optional discrete set of candidates (otherwise a maximin LHS is used to initialise continuous search)
##' @param control list in case \code{Xcand == NULL}, with elements \code{multi.start},
##' to perform a multi-start optimization based on \code{\link[stats]{optim}}, with \code{maxit} iterations each.
##' Also, \code{tol_dist} defines the minimum distance to an existing design for a new point to be added, otherwise the closest existing design is chosen.
##' In a similar fashion, \code{tol_dist} is the minimum relative change of IMSPE for adding a new design.
#' @param h horizon for multi-step ahead framework.
#' The decision is made between:
#' \itemize{
#'  \item sequential crit search starting by a new design (optimized first) then adding \code{h} replicates
#'  \item sequential crit searches starting by \code{1} to \code{h} replicates before adding a new point
#' }
#' Use \code{h = 0} for the myopic criterion, i.e., not looking ahead.
##' @param Wijs optional previously computed matrix of Wijs, see \code{\link[hetGP]{Wij}}
##' @param seed optional seed for the generation of designs with \code{\link[DiceDesign]{maximinSA_LHS}}
##' @param ncores number of CPU available (> 1 mean parallel TRUE), see \code{\link[parallel]{mclapply}}
##' @details The domain needs to be [0, 1]^d for now.
##' @return list with elements:
##' \itemize{
##' \item \code{par}: best first design,
##' \item \code{value}: IMSPE h-steps ahead starting from adding \code{par},
##' \item \code{path}: list of elements list(\code{par}, \code{value}, \code{new}) at each step \code{h}
##' }
##' @references M. Binois, J. Huang, R. Gramacy, M. Ludkovski (2017+), Replication or exploration? Sequential design for stochastic simulation experiments. arXiv preprint arXiv:1710.03206.
##' @export 
##' @examples
##' ###############################################################################
##' ## Bi-variate example (myopic version)
##' ###############################################################################
##' 
##' nvar <- 2 
##' 
##' set.seed(42)
##' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
##' 
##' n <- 25 # must be a square
##' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
##' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
##' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
##' Z <- apply(X, 1, ftest)
##' 
##' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
##' 
##' ngrid <- 51
##' xgrid <- seq(0,1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' 
##' preds <- predict(x = Xgrid, object =  model)
##' 
##' ## Initial plots
##' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
##'         main = "Predicted mean", nlevels = 20)
##' points(model$X0, col = 'blue', pch = 20)
##' 
##' IMSPE_grid <- apply(Xgrid, 1, crit_IMSPE, model = model)
##' filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid),
##'                nlevels = 20, color.palette = terrain.colors, 
##'                main = "Initial IMSPE criterion landscape",
##' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
##' 
##' ## Sequential IMSPE search
##' nsteps <- 1 # Increase for better results
##' 
##' for(i in 1:nsteps){
##'   res <- IMSPE_optim(model, control = list(multi.start = 30, maxit = 30))
##'   newX <- res$par
##'   newZ <- ftest(newX)
##'   model <- update(object = model, Xnew = newX, Znew = newZ)
##' }
##' 
##' ## Final plots
##' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
##'         main = "Predicted mean", nlevels = 20)
##' points(model$X0, col = 'blue', pch = 20)
##' 
##' IMSPE_grid <- apply(Xgrid, 1, crit_IMSPE, model = model)
##' filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid),
##'                nlevels = 20, color.palette = terrain.colors, 
##'                main = "Final IMSPE criterion landscape",
##' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
##' 
##' ###############################################################################
##' ## Bi-variate example (look-ahead version)
##' ###############################################################################
##' \dontrun{ 
##' nvar <- 2 
##' 
##' set.seed(42)
##' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
##' 
##' n <- 25 # must be a square
##' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
##' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
##' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
##' Z <- apply(X, 1, ftest)
##' 
##' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
##' 
##' ngrid <- 51
##' xgrid <- seq(0,1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' 
##' nsteps <- 5 # Increase for more steps
##' 
##' # To use parallel computation (turn off on Windows)
##' library(parallel)
##' parallel <- FALSE #TRUE #
##' if(parallel) ncores <- detectCores() else ncores <- 1
##' 
##' for(i in 1:nsteps){
##'   res <- IMSPE_optim(model, h = 3, control = list(multi.start = 100, maxit = 50),
##'    ncores = ncores)
##'   
##'   # If a replicate is selected
##'   if(!res$path[[1]]$new) print("Add replicate")
##'   
##'   newX <- res$par
##'   newZ <- ftest(newX)
##'   model <- update(object = model, Xnew = newX, Znew = newZ)
##'   
##'   ## Plots 
##'   preds <- predict(x = Xgrid, object =  model)
##'   contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
##'           main = "Predicted mean", nlevels = 20)
##'   points(model$X0, col = 'blue', pch = 20)
##'   points(newX, col = "red", pch = 20)
##'   
##'   ## Precalculations
##'   Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
##'   
##'   IMSPE_grid <- apply(Xgrid, 1, crit_IMSPE, Wijs = Wijs, model = model)
##'   filled.contour(x = xgrid, y = xgrid, matrix(IMSPE_grid, ngrid),
##'                  nlevels = 20, color.palette = terrain.colors,
##'   plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
##' }
##' }
IMSPE_optim <- function(model, h = 2, Xcand = NULL, control = list(tol_dist = 1e-6, tol_diff = 1e-6, multi.start = 20, maxit = 100),
                        Wijs = NULL, seed = NULL, ncores = 1){
  d <- ncol(model$X0)
  
  if(max(rbind(model$X0, Xcand)) > 1) stop("IMSPE works only with [0,1]^d domain for now.")
  if(max(rbind(model$X0, Xcand)) < 0) stop("IMSPE works only with [0,1]^d domain for now.")
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
  
  
  ## A) Setting to beat: first new point then replicate h times
  IMSPE_A <- IMSPE.search(model = model, control = control, Xcand = Xcand, Wijs = Wijs, ncores = ncores)
  new_designA <- IMSPE_A$par ## store first considered design to be added
  path_A <- list(IMSPE_A)
  
  if(h > 0){
    newmodelA <- model
    
    if(IMSPE_A$new){
      newWijs <- Wij(mu1 = model$X0, mu2 = new_designA, theta = model$theta, type = model$covtype)
      WijsA <- cbind(Wijs, newWijs)
      WijsA <- rbind(WijsA, c(newWijs, Wij(IMSPE_A$par, IMSPE_A$par, theta = model$theta, type = model$covtype)))
    }else{
      WijsA <- Wijs
    }
    
    for(i in 1:h){
      newmodelA <- update(object = newmodelA, Xnew = IMSPE_A$par, Znew = NA, maxit = 0)
      IMSPE_A <- IMSPE.search(model = newmodelA, replicate = TRUE, control = control, Wijs = WijsA, ncores = ncores)
      path_A <- c(path_A, list(IMSPE_A))
    }

  }
  
  if(h == -1) return(list(par = new_designA, value = IMSPE_A$value, path = path_A)) 
  
  ## B) Now compare with waiting to add new point
  newmodelB <- model
  
  if(h == 0){
    IMSPE_B <- IMSPE.search(model = newmodelB, replicate = TRUE, control = control, Wijs = Wijs, ncores = ncores)
    new_designB <- IMSPE_B$par ## store considered design to be added
    
    # search from best replicate
    if(is.null(Xcand)){
      IMSPE_C <- IMSPE.search(model = newmodelB, Wijs = Wijs,
                            control = list(Xstart = IMSPE_B$par, maxit = control$maxit,
                                           tol_dist = control$tol_dist, tol_diff = control$tol_diff), ncores = ncores)
    }else{
      IMSPE_C <- IMSPE_B
    }
    
    if(IMSPE_C$value < min(IMSPE_A$value, IMSPE_B$value)) return(list(par = IMSPE_C$par, value = IMSPE_C$value, path = list(IMSPE_C)))
    
    if(IMSPE_B$value < IMSPE_A$value){
      return(list(par = IMSPE_B$par, value = IMSPE_B$value, path = list(IMSPE_B)))
    } 
  }else{
    for(i in 1:h){
      ## Add new replicate
      IMSPE_B <- IMSPE.search(model = newmodelB, replicate = TRUE, control = control, Wijs = Wijs, ncores = ncores)
      
      if(i == 1){
        new_designB <- matrix(IMSPE_B$par, nrow = 1) ##store first considered design to add
        path_B <- list()
      } 
      
      path_B <- c(path_B, list(IMSPE_B))
      newmodelB <- update(object = newmodelB, Xnew = IMSPE_B$par, Znew = NA, maxit = 0)
      
      ## Add new design
      IMSPE_C <- IMSPE.search(model = newmodelB, control = control, Xcand = Xcand, Wijs = Wijs, ncores = ncores)
      path_C <- list(IMSPE_C)
      
      if(i < h){
        newmodelC <- newmodelB
        
        if(!any(duplicated(rbind(model$X0, IMSPE_C$par)))){ 
          newWijs <- Wij(mu1 = model$X0, mu2 = IMSPE_C$par, theta = model$theta, type = model$covtype)
          WijsC <- cbind(Wijs, newWijs)
          WijsC <- rbind(WijsC, c(newWijs, Wij(mu1 = IMSPE_C$par, theta = model$theta, type = model$covtype)))
        }else{
          WijsC <- Wijs
        }
        
        for(j in i:(h-1)){
          ## Add remaining replicates
          newmodelC <- update(object = newmodelC, Xnew = IMSPE_C$par, Znew = NA, maxit = 0)
          IMSPE_C <- IMSPE.search(model = newmodelC, replicate = TRUE, control = control, Wijs = WijsC, ncores = ncores)
          path_C <- c(path_C, list(IMSPE_C))
        }
      }
      
      if(IMSPE_C$value < IMSPE_A$value) return(list(par = new_designB, value = IMSPE_C$value, path = c(path_B, path_C)))
    }
  }
  
  return(list(par = new_designA, value = IMSPE_A$value, path = path_A))
  
}


##' Allocation of replicates on existing design locations, based on (29) from (Ankenman et al, 2010)
##' @title Allocation of replicates on existing designs
##' @param model \code{hetGP} model
##' @param N total budget of replication to allocate
##' @param Wijs optional previously computed matrix of \code{Wijs}, see \code{\link[hetGP]{Wij}}
##' @param use.Ki should \code{Ki} from \code{model} be used? 
##' Using the inverse of C (covariance matrix only, without noise, using \code{\link[MASS]{ginv}}) is also possible
##' @return vector with approximated best number of replicates per design
##' @references B. Ankenman, B. Nelson, J. Staum (2010), Stochastic kriging for simulation metamodeling, Operations research, pp. 371--382, 58
##' @export
##' @examples
##' ##------------------------------------------------------------
##' ## Example: Heteroskedastic GP modeling on the motorcycle data
##' ##------------------------------------------------------------
##' set.seed(32)
##' 
##' ## motorcycle data
##' library(MASS)
##' X <- matrix(mcycle$times, ncol = 1)
##' Z <- mcycle$accel
##' nvar <- 1
##' 
##' data_m <- find_reps(X, Z, rescale = TRUE)
##' 
##' plot(rep(data_m$X0, data_m$mult), data_m$Z, ylim = c(-160, 90),
##'      ylab = 'acceleration', xlab = "time")
##' 
##' 
##' ## Model fitting
##' model <- mleHetGP(X = list(X0 = data_m$X0, Z0 = data_m$Z0, mult = data_m$mult),
##'                   Z = Z, lower = rep(0.1, nvar), upper = rep(5, nvar),
##'                   covtype = "Matern5_2")
##' ## Compute best allocation                  
##' A <- allocate_mult(model, N = 1000)
##' 
##' ## Create a prediction grid and obtain predictions
##' xgrid <- matrix(seq(0, 1, length.out = 301), ncol = 1) 
##' predictions <- predict(x = xgrid, object =  model)
##' 
##' ## Display mean predictive surface
##' lines(xgrid, predictions$mean, col = 'red', lwd = 2)
##' ## Display 95% confidence intervals
##' lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
##' lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
##' ## Display 95% prediction intervals
##' lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
##' col = 3, lty = 2)
##' lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
##' col = 3, lty = 2)
##' 
##' par(new = TRUE)
##' plot(NA,NA, xlim = c(0,1), ylim = c(0,max(A)), axes = FALSE, ylab = "", xlab = "")
##' segments(x0 = model$X0, x1 = model$X0, 
##' y0 = rep(0, nrow(model$X)), y1 = A, col = 'grey')
##' axis(side = 4)
##' mtext(side = 4, line = 2, expression(a[i]), cex = 0.8)       
allocate_mult <- function(model, N, Wijs = NULL, use.Ki = FALSE){
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
  
  if(use.Ki){
    Ci <- pmax(0, diag(model$Ki %*% Wijs %*% model$Ki)) # pmax for numerical imprecision
  }else{
    Ci <- ginv(cov_gen(model$X0, theta = model$theta, type = model$covtype))
    Ci <- pmax(0, diag(Ci %*% Wijs %*% Ci))
  }
  
  if(class(model) == "hetGP") V <- model$Lambda
  else V <- rep(model$g, length(model$Z0))
  
  res <- N * sqrt(Ci * V) / sum(sqrt(Ci * V))
  
  # Now get integer results summing to N
  idxs <- sort(res %% 1, index.return = T, decreasing = T)
  bdg <- N - sum(floor(res))  # remaining number of points to add after truncating
  
  res <- floor(res)
  if(bdg > 0)
    res[idxs$ix[1:bdg]] <- res[idxs$ix[1:bdg]] + 1
  
  return(res)
}

##' Adapt the look-ahead horizon depending on the replicate allocation or a target ratio
##' @title Adapt horizon
##' @param model \code{hetGP} or \code{homGP} model
##' @param current_horizon horizon used for the previous iteration, see details
##' @param previous_ratio ratio before adding the previous new design
##' @param target scalar in ]0,1] for desired n/N
##' @param Wijs optional previously computed matrix of Wijs, see \code{\link[hetGP]{Wij}}
##' @return randomly selected horizon for next iteration (adpative) if no \code{target} is provided, 
##' otherwise returns the update horizon value.
##' @details 
##' If \code{target} is provided, along with \code{previous_ratio} and \code{current_horizon}:
##' \itemize{
##' \item the horizon is increased by one if more replicates are needed but a new ppint has been added at the previous iteration,
##' \item the horizon is decreased by one if new points are needed but a replicate has been added at the previous iteration,
##' \item otherwise it is unchanged.
##' }
##' 
##' If no \code{target} is provided, \code{\link[hetGP]{allocate_mult}} is used to obtain the best allocation of the existing replicates,
##' then the new horizon is sampled from the difference between the actual allocation and the best one, bounded below by 0.
##' See (Binois et al. 2017).
##' 
##' @references
##' M. Binois, J. Huang, R. B. Gramacy, M. Ludkovski (2018+), Replication or exploration? Sequential design for stochastic simulation experiments,
##' Technometrics (to appear).\cr 
##' Preprint available on arXiv:1710.03206.
##' @export
horizon <- function(model, current_horizon = NULL, previous_ratio = NULL, target = NULL, Wijs = NULL){
  if(is.null(target)){
    mult_star <- allocate_mult(model = model, N = sum(model$mult), Wijs = Wijs)
    tab <- table(pmax(mult_star - model$mult, 0))
    return(sample(as.numeric(names(tab)), 1, prob=tab))
  }
  
  if(is.null(current_horizon) || is.null(previous_ratio)) cat("Missing arguments to use target \n")
  
  ratio <- length(model$Z0)/length(model$Z)
  
  # Ratio increased while too small
  if(ratio < target && ratio < previous_ratio){
    return(max(-1, current_horizon - 1))
  }
  # Ratio decreased while too high
  if(ratio > target && ratio > previous_ratio)
    return(current_horizon + 1)
    
  return(current_horizon)
}

##' Compute double integral of the covariance kernel over a [0,1]^d domain
##' @param mu1,mu2 input locations considered
##' @param theta lengthscale hyperparameter of the kernel
##' @param type kernel type, one of "\code{Gaussian}", "\code{Matern5_2}" or "\code{Matern3_2}", see \code{\link[hetGP]{cov_gen}}
##' @export
##' @references M. Binois, J. Huang, R. Gramacy, M. Ludkovski (2017+), Replication or exploration? Sequential design for stochastic simulation experiments.
Wij <- function(mu1, mu2 = NULL, theta, type){
  if(ncol(mu1) > 1 && length(theta) == 1) theta <- rep(theta, ncol(mu1))
  
  if(type == "Gaussian"){
    if(is.null(mu2))
      return(Wijs_gauss_sym_cpp(mu1, sqrt(theta)))
    return(Wijs_gauss_cpp(mu1, mu2, sqrt(theta)))
  }
  if(type == "Matern5_2"){
    if(is.null(mu2))
      return(Wijs_mat52_sym_cpp(mu1, theta))
    return(Wijs_mat52_cpp(mu1, mu2, theta))
  }
  if(type == "Matern3_2"){
    if(is.null(mu2))
      return(Wijs_mat32_sym_cpp(mu1, theta))
    return(Wijs_mat32_cpp(mu1, mu2, theta))
  }
}

##' Compute integral of the covariance kernel over a [0,1]^d domain
##' @param mu1 input locations considered
##' @param theta lengthscale hyperparameter of the kernel
##' @param type kernel type, one of "\code{Gaussian}", "\code{Matern5_2}" or "\code{Matern3_2}", see \code{\link[hetGP]{cov_gen}}
##' @noRd
##' @references M. Binois, J. Huang, R. Gramacy, M. Ludkovski (2017+), Replication or exploration? Sequential design for stochastic simulation experiments.
mi <- function(mu1, theta, type){
  if(ncol(mu1) > 1 && length(theta) == 1) theta <- rep(theta, ncol(mu1))
  
  if(type == "Gaussian"){
    return(mi_gauss_cpp(mu1, sqrt(theta)))
  }
  if(type == "Matern5_2"){
    return(mi_mat52_cpp(mu1, theta))
  }
  if(type == "Matern3_2"){
    return(mi_mat32_cpp(mu1, theta))
  }
}
