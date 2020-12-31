##' Computes MEE infill criterion
##' @title Maximum Empirical Error criterion
##' @param x matrix of new designs, one point per row (size n x d)
##' @param model \code{homGP} or \code{hetGP} model, including inverse matrices
##' @param thres for contour finding
##' @param preds optional predictions at \code{x} to avoid recomputing if already done
##' @export
##' @importFrom stats pt pnorm dnorm
##' @references
##' Ranjan, P., Bingham, D. & Michailidis, G (2008). 
##' Sequential experiment design for contour estimation from complex computer codes, 
##' Technometrics, 50, pp. 527-541. \cr \cr
##' 
##' Bichon, B., Eldred, M., Swiler, L., Mahadevan, S. & McFarland, J. (2008).
##' Efficient global  reliability  analysis  for  nonlinear  implicit  performance  functions, 
##' AIAA Journal, 46, pp. 2459-2468. \cr \cr
##' 
##' Lyu, X., Binois, M. & Ludkovski, M. (2018). Evaluating Gaussian Process Metamodels and Sequential Designs for Noisy Level Set Estimation. arXiv:1807.06712. \cr
##' 
##' @examples 
##' ## Infill criterion example
##' set.seed(42)
##' branin <- function(x){
##'   m <- 54.8104; s <- 51.9496
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   xx <- 15 * x[,1] - 5; y <- 15 * x[,2]
##'   f <- (y - 5.1 * xx^2/(4 * pi^2) + 5 * xx/pi - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(xx) + 10
##'   f <- (f - m)/s
##'   return(f)
##' }
##' 
##' ftest <- function(x, sd = 0.1){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   return(apply(x, 1, branin) + rnorm(nrow(x), sd = sd))
##' }
##' 
##' ngrid <- 101; xgrid <- seq(0, 1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' Zgrid <- ftest(Xgrid)
##' 
##' n <- 20
##' N <- 500
##' X <- Xgrid[sample(1:nrow(Xgrid), n),]
##' X <- X[sample(1:n, N, replace = TRUE),]
##' Z <- ftest(X)
##' model <- mleHetGP(X, Z, lower = rep(0.001,2), upper = rep(1,2))
##' 
##' critgrid <- apply(Xgrid, 1, crit_MEE, model = model)
##' 
##' filled.contour(matrix(critgrid, ngrid), color.palette = terrain.colors, main = "MEE criterion")
##' 
crit_MEE <- function(x, model, thres = 0, preds = NULL){
  
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  ## TP case
  if(class(model) %in% c("homTP", "hetTP")){
    return(pt(-abs(preds$mean - thres)/sqrt(preds$sd2), df = model$nu + length(model$Z)))
  }
  
  ## GP case
  return(pnorm(-abs(preds$mean - thres)/sqrt(preds$sd2)))
}


##' Computes cSUR infill criterion
##' @title Contour Stepwise Uncertainty Reduction criterion
##' @param x matrix of new designs, one point per row (size n x d)
##' @param model \code{homGP} or \code{hetGP} model, including inverse matrices
##' @param thres for contour finding
##' @param preds optional predictions at \code{x} to avoid recomputing if already done (must contain \code{cov})
##' @references
##' Lyu, X., Binois, M. & Ludkovski, M. (2018). Evaluating Gaussian Process Metamodels and Sequential Designs for Noisy Level Set Estimation. arXiv:1807.06712. \cr
##' @export
## ' @details TODO: deal with replication
##' @examples 
##' ## Infill criterion example
##' set.seed(42)
##' branin <- function(x){
##'   m <- 54.8104; s <- 51.9496
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   xx <- 15 * x[,1] - 5; y <- 15 * x[,2]
##'   f <- (y - 5.1 * xx^2/(4 * pi^2) + 5 * xx/pi - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(xx) + 10
##'   f <- (f - m)/s
##'   return(f)
##' }
##' 
##' ftest <- function(x, sd = 0.1){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   return(apply(x, 1, branin) + rnorm(nrow(x), sd = sd))
##' }
##' 
##' ngrid <- 101; xgrid <- seq(0, 1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' Zgrid <- ftest(Xgrid)
##' 
##' n <- 20
##' N <- 500
##' X <- Xgrid[sample(1:nrow(Xgrid), n),]
##' X <- X[sample(1:n, N, replace = TRUE),]
##' Z <- ftest(X)
##' model <- mleHetGP(X, Z, lower = rep(0.001,2), upper = rep(1,2))
##' 
##' critgrid <- apply(Xgrid, 1, crit_cSUR, model = model)
##' 
##' filled.contour(matrix(critgrid, ngrid), color.palette = terrain.colors, main = "cSUR criterion")
##' 
crit_cSUR <- function(x, model, thres = 0, preds = NULL){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds) || is.null(preds$cov)) preds <- predict(model, x = x, xprime = x)
  
  if(class(model) %in% c("homTP", "hetTP")){
    
    # unscale the predictive variance and covariance (e.g., go back to the GP case)
    # (since psi is updated separately)
    pcov <- preds$cov * (model$nu + length(model$Z) - 2) / (model$nu + model$psi - 2)
    psd2 <- preds$sd2 * (model$nu + length(model$Z) - 2) / (model$nu + model$psi - 2)
    
    Ki_new <- chol2inv(chol(pcov + diag(model$eps + preds$nugs, nrow = length(preds$nugs))))
    sd2_new <- pmax(0, psd2 - fast_diag(pcov, tcrossprod(Ki_new, pcov)))
    
    # now update psi
    # kn1 <- model$sigma2 * cov_gen(model$X0, x, theta = model$theta, type = model$covtype)
    # hn <- - crossprod(model$Z0, model$Ki) %*% kn1
    # psi_n1 <- model$psi + (hn^2 + 2 * preds$mean * hn + preds$mean^2 + model$nu/(model$nu - 2)*preds$sd2)/preds$sd2
    psi_n1 <- model$psi + model$nu/(model$nu - 2)
    
    sd2_new <- (model$nu + psi_n1 - 2) / (model$nu + length(model$Z) - 1) * sd2_new # unscaled variance
    
    return(pt(-abs(preds$mean - thres)/sqrt(preds$sd2), df = model$nu + length(model$Z)) - 
             pt(-abs(preds$mean - thres)/sqrt(sd2_new), df = model$nu + length(model$Z) + 1))
  }else{
    
    Ki_new <- chol2inv(chol(preds$cov + diag(model$eps + preds$nugs, nrow = length(preds$nugs))))
    sd2_new <- pmax(0, preds$sd2 - fast_diag(preds$cov, tcrossprod(Ki_new, preds$cov)))
    
    return(pnorm(-abs(preds$mean - thres)/sqrt(preds$sd2)) - pnorm(-abs(preds$mean - thres)/sqrt(sd2_new)))
  }
  
}

##' Computes ICU infill criterion
##' @title Integrated Contour Uncertainty criterion
##' @param x matrix of new designs, one point per row (size n x d)
##' @param model \code{homGP} or \code{hetGP} model, including inverse matrices
##' @param Xref matrix of input locations to approximate the integral by a sum
##' @param w  optional weights vector of weights for \code{Xref} locations
##' @param thres for contour finding
##' @param preds optional predictions at \code{Xref} to avoid recomputing if already done
##' @param kxprime optional covariance matrix between \code{model$X0} and \code{Xref} to avoid its recomputation
##' @references
##' Lyu, X., Binois, M. & Ludkovski, M. (2018). Evaluating Gaussian Process Metamodels and Sequential Designs for Noisy Level Set Estimation. arXiv:1807.06712. \cr
##' @export
## ' @details TODO: deal with replication
##' @examples 
##' ## Infill criterion example
##' set.seed(42)
##' branin <- function(x){
##'   m <- 54.8104; s <- 51.9496
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   xx <- 15 * x[,1] - 5; y <- 15 * x[,2]
##'   f <- (y - 5.1 * xx^2/(4 * pi^2) + 5 * xx/pi - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(xx) + 10
##'   f <- (f - m)/s
##'   return(f)
##' }
##' 
##' ftest <- function(x, sd = 0.1){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   return(apply(x, 1, branin) + rnorm(nrow(x), sd = sd))
##' }
##' 
##' ngrid <- 51; xgrid <- seq(0, 1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' Zgrid <- ftest(Xgrid)
##' 
##' n <- 20
##' N <- 500
##' X <- Xgrid[sample(1:nrow(Xgrid), n),]
##' X <- X[sample(1:n, N, replace = TRUE),]
##' Z <- ftest(X)
##' model <- mleHetGP(X, Z, lower = rep(0.001,2), upper = rep(1,2))
##' 
##' # Precalculations for speedup
##' preds <- predict(model, x = Xgrid)
##' kxprime <- cov_gen(X1 = model$X0, X2 = Xgrid, theta = model$theta, type = model$covtype)
##'  
##' critgrid <- apply(Xgrid, 1, crit_ICU, model = model, Xref = Xgrid,
##'                   preds = preds, kxprime = kxprime)
##' 
##' filled.contour(matrix(critgrid, ngrid), color.palette = terrain.colors, main = "ICU criterion")
##' 
crit_ICU <- function(x, model, thres = 0, Xref, w = NULL, preds = NULL, kxprime = NULL){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = Xref)
  if(is.null(w)) w <- rep(1, nrow(Xref))
  
  predx <- predict(model, x = x)
  if(is.null(kxprime)){
    covnew <- predict(model, x = x, xprime = Xref)$cov
  }else{
    if(class(model) %in% c("homTP", "hetTP")){
      kxprime <- kxprime * model$sigma2
      kx <- model$sigma2 * cov_gen(X1 = x, X2 = model$X0, theta = model$theta, type = model$covtype)
      covnew <- (model$nu + model$psi - 2) / (model$nu + length(model$Z) - 2) * (model$sigma2 * cov_gen(X1 = x, X2 = Xref, theta = model$theta, type = model$covtype) - (kx %*% model$Ki) %*% kxprime)
    }else{
      kx <- model$nu_hat * cov_gen(X1 = x, X2 = model$X0, theta = model$theta, type = model$covtype)
      model$Ki <- model$Ki / model$nu_hat
      kxprime <- kxprime * model$nu_hat
      if(model$trendtype == 'SK'){
        covnew <- model$nu_hat * cov_gen(X1 = x, X2 = Xref, theta = model$theta, type = model$covtype) - (kx %*% model$Ki) %*% kxprime
      }else{
        covnew <- model$nu_hat * cov_gen(X1 = x, X2 = Xref, theta = model$theta, type = model$covtype) - (kx %*% model$Ki) %*% kxprime + crossprod(1 - tcrossprod(rowSums(model$Ki), kx), 1 - rowSums(model$Ki) %*% kxprime)/sum(model$Ki)
      }
    }
  }
  
  if(class(model) %in% c("homTP", "hetTP")){
    # unscale the predictive variance and covariances (e.g., go back to the GP case)
    # (since psi is updated separately)
    predx$sd2 <- (model$nu + length(model$Z) - 2) / (model$nu + model$psi - 2) * predx$sd2
    covnew <- (model$nu + length(model$Z) - 2) / (model$nu + model$psi - 2) * covnew
    preds$sd2 <- (model$nu + length(model$Z) - 2) / (model$nu + model$psi - 2) * preds$sd2
    
    sd2_new <- pmax(0, preds$sd2 - drop(covnew^2)/(predx$sd2 + predx$nugs + model$eps))
    
    # now update psi
    psi_n1 <- model$psi + model$nu/(model$nu - 2)
    
    # now rescale with updated psi
    sd2_new <- (model$nu + psi_n1 - 2) / (model$nu + length(model$Z) - 1) * sd2_new
    
    return(- sum(w * pt(-abs(preds$mean - thres)/sqrt(sd2_new), df = model$nu + length(model$Z) + 1)))
  }else{
    sd2_new <- pmax(0, preds$sd2 - drop(covnew^2)/(predx$sd2 + predx$nugs + model$eps))
    return(- sum(w * pnorm(-abs(preds$mean - thres)/sqrt(sd2_new))))
  }
  
}


##' Computes targeted mean squared error infill criterion
##' @title t-MSE criterion
##' @param x matrix of new designs, one point per row (size n x d)
##' @param model \code{homGP} or \code{hetGP} model, including inverse matrices
##' @param thres for contour finding
##' @param preds optional predictions at \code{x} to avoid recomputing if already done (must contain \code{cov})
##' @param seps parameter for the target window
##' @references
##' Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R., Kim, N. (2010).
##' Adaptive designs of experiments for accurate approximation of a target region,
##' Journal of Mechanical Design (132), p. 071008.\cr \cr
##' 
##' Lyu, X., Binois, M. & Ludkovski, M. (2018). 
##' Evaluating Gaussian Process Metamodels and Sequential Designs for Noisy Level Set Estimation. arXiv:1807.06712. \cr
##' @export
##' @examples 
##' ## Infill criterion example
##' set.seed(42)
##' branin <- function(x){
##'   m <- 54.8104; s <- 51.9496
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   xx <- 15 * x[,1] - 5; y <- 15 * x[,2]
##'   f <- (y - 5.1 * xx^2/(4 * pi^2) + 5 * xx/pi - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(xx) + 10
##'   f <- (f - m)/s
##'   return(f)
##' }
##' 
##' ftest <- function(x, sd = 0.1){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   return(apply(x, 1, branin) + rnorm(nrow(x), sd = sd))
##' }
##' 
##' ngrid <- 101; xgrid <- seq(0, 1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' Zgrid <- ftest(Xgrid)
##' 
##' n <- 20
##' N <- 500
##' X <- Xgrid[sample(1:nrow(Xgrid), n),]
##' X <- X[sample(1:n, N, replace = TRUE),]
##' Z <- ftest(X)
##' model <- mleHetGP(X, Z, lower = rep(0.001,2), upper = rep(1,2))
##' 
##' critgrid <- apply(Xgrid, 1, crit_tMSE, model = model)
##' 
##' filled.contour(matrix(critgrid, ngrid), color.palette = terrain.colors, main = "tMSE criterion")
##' 
crit_tMSE <- function(x, model, thres = 0, preds = NULL, seps = 0.05){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds) || is.null(preds$cov)) preds <- predict(model, x = x, xprime = x)
  
  w <- 1/sqrt(2 * pi * (preds$sd2 + seps)) * exp(-0.5 * (preds$mean - thres)^2 / (preds$sd2 + seps))
  
  return(w * preds$sd2)
}


##' Computes MCU infill criterion
##' @title Maximum Contour Uncertainty criterion
##' @param x matrix of new designs, one point per row (size n x d)
##' @param model \code{homGP} or \code{hetGP} model, including inverse matrices
##' @param thres for contour finding
##' @param gamma optional weight in -|f(x) - thres| + gamma * s(x). Default to 2.
##' @param preds optional predictions at \code{x} to avoid recomputing if already done
##' @export
##' @importFrom stats pt pnorm dnorm
##' @references
##' Srinivas, N., Krause, A., Kakade, S, & Seeger, M. (2012). 
##' Information-theoretic regret bounds for Gaussian process optimization 
##' in the bandit setting, IEEE Transactions on Information Theory, 58, pp. 3250-3265.\cr \cr
##' 
##' Bogunovic, J., Scarlett, J., Krause, A. & Cevher, V. (2016). 
##' Truncated variance reduction: A unified approach to Bayesian optimization and level-set estimation,
##' in Advances in neural information processing systems, pp. 1507-1515. \cr \cr
##' 
##' Lyu, X., Binois, M. & Ludkovski, M. (2018). 
##' Evaluating Gaussian Process Metamodels and Sequential Designs for Noisy Level Set Estimation. arXiv:1807.06712. \cr
##' 
##' @examples 
##' ## Infill criterion example
##' set.seed(42)
##' branin <- function(x){
##'   m <- 54.8104; s <- 51.9496
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   xx <- 15 * x[,1] - 5; y <- 15 * x[,2]
##'   f <- (y - 5.1 * xx^2/(4 * pi^2) + 5 * xx/pi - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(xx) + 10
##'   f <- (f - m)/s
##'   return(f)
##' }
##' 
##' ftest <- function(x, sd = 0.1){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1)
##'   return(apply(x, 1, branin) + rnorm(nrow(x), sd = sd))
##' }
##' 
##' ngrid <- 101; xgrid <- seq(0, 1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' Zgrid <- ftest(Xgrid)
##' 
##' n <- 20
##' N <- 500
##' X <- Xgrid[sample(1:nrow(Xgrid), n),]
##' X <- X[sample(1:n, N, replace = TRUE),]
##' Z <- ftest(X)
##' model <- mleHetGP(X, Z, lower = rep(0.001,2), upper = rep(1,2))
##' 
##' critgrid <- apply(Xgrid, 1, crit_MCU, model = model)
##' 
##' filled.contour(matrix(critgrid, ngrid), color.palette = terrain.colors, main = "MEE criterion")
##' 
crit_MCU <- function(x, model, thres = 0, gamma = 2, preds = NULL){
  
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  

  ## TP case
  if(class(model) %in% c("homTP", "hetTP")){
    return(-abs(preds$mean - thres) + gamma * sqrt(preds$sd2))
  }
  
  ## GP case
  return(-abs(preds$mean - thres) + gamma * sqrt(preds$sd2))
}
