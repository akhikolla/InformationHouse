##' Performs Gaussian process regression with heteroskedastic noise following 
##' Binois, M., Gramacy, R., Ludkovski, M. (2016) <arXiv:1611.05902>. 
##' The input dependent noise is modeled as another Gaussian process. 
##' Replicated observations are encouraged as they yield computational savings. 
##' Sequential design procedures based on the integrated mean square prediction error and lookahead heuristics are provided,
##' and notably fast update functions when adding new observations.
##' @title Package hetGP
##' @author Mickael Binois, Robert B. Gramacy
##' @docType package
##' @name hetGP-package
##' @references 
##' M. Binois, Robert B. Gramacy, M. Ludkovski (2018), Practical heteroskedastic Gaussian process modeling for large simulation experiments,
##' Journal of Computational and Graphical Statistics, 27(4), 808--821.\cr 
##' Preprint available on arXiv:1611.05902. \cr \cr
##' M. Binois, J. Huang, R. B. Gramacy, M. Ludkovski (2018+), Replication or exploration? Sequential design for stochastic simulation experiments,
##' Technometrics (to appear).\cr 
##' Preprint available on arXiv:1710.03206.
##' @details 
##' Important functions: \cr
##' \code{\link[hetGP]{mleHetGP}} as the main function to build a model. \cr
##' \code{\link[hetGP]{mleHomGP}} the equivalent for homoskedastic modeling. \cr
##' \code{\link[hetGP]{crit_IMSPE}} for adding a new design based on the Integrated Mean Square Prediction Error. \cr
##' \code{\link[hetGP]{IMSPE_optim}} for augmenting a design, possibly based on a lookahead heuristic to bias the search toward replication. \cr
##' \code{\link[hetGP]{crit_optim}} is similar to \code{IMSPE_optim} but for the optimization or contour finding criterion available.
##' @note 
##' The authors are grateful for support from National Science Foundation grant DMS-1521702 and DMS-1521743.
##' @examples 
##' ##------------------------------------------------------------
##' ## Example 1: Heteroskedastic GP modeling on the motorcycle data
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
##' ## Model fitting
##' settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training
##' model <- mleHetGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar),
##'                   covtype = "Matern5_2", settings = settings)
##' 
##' ## A quick view of the fit                  
##' summary(model)
##' 
##' ## Create a prediction grid and obtain predictions
##' xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 
##' predictions <- predict(x = xgrid, object =  model)
##' 
##' ## Display averaged observations
##' points(model$X0, model$Z0, pch = 20)
##' 
##' ## Display mean predictive surface
##' lines(xgrid, predictions$mean, col = 'red', lwd = 2)
##' ## Display 95% confidence intervals
##' lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
##' lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
##' ## Display 95% prediction intervals
##' lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
##'   col = 3, lty = 2)
##' lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
##'   col = 3, lty = 2)
##'   
##' ## Comparison with homoskedastic fit
##' predictions2 <- predict(x = xgrid, object = model$modHom)
##' lines(xgrid, predictions2$mean, col = 4, lty = 2, lwd = 2)
##' lines(xgrid, qnorm(0.05, predictions2$mean, sqrt(predictions2$sd2)), col = 4, lty = 3)
##' lines(xgrid, qnorm(0.95, predictions2$mean, sqrt(predictions2$sd2)), col = 4, lty = 3)
##' 
##' 
##' ##------------------------------------------------------------
##' ## Example 2: Sequential design
##' ##------------------------------------------------------------
##' \dontrun{
##' library(DiceDesign)
##'
##' ## Design configuration / Parameter settings
##' N_tot <- 500 # total number of points
##' n_init <- 10 # number of unique designs
##'
##' ## HetGP options
##' nvar <- 1 # number of variables
##' lower <- rep(0.001, nvar)
##' upper <- rep(1, nvar)
##'
##' ### Problem definition
##'
##' ## Mean function
##' forrester <- function(x){
##'  return(((x*6-2)^2)*sin((x*6-2)*2))
##' }
##'
##' ## Noise field via standard deviation
##' noiseFun <- function(x, coef = 1.1, scale = 1){
##'  if(is.null(nrow(x)))
##'    x <- matrix(x, nrow = 1)
##'  return(scale*(coef + sin(x * 2 * pi)))
##' }
##'
##' ### Test function defined in [0,1]
##' ftest <- function(x){
##'  if(is.null(nrow(x)))
##'    x <- matrix(x, ncol = 1)
##'  return(forrester(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
##' }
##'
##' ## Predictive grid
##' ngrid <- 51
##' xgrid <- seq(0,1, length.out = ngrid)
##' Xgrid <- matrix(xgrid, ncol = 1)
##'
##' par(mar = c(3,3,2,3)+0.1)
##' plot(xgrid, forrester(xgrid), type = 'l', lwd = 1, col = "blue", lty = 3,
##'     xlab = '', ylab = '', ylim = c(-8,16))
##'
##' set.seed(42)
##'
##' # Initial design
##' X <- maximinSA_LHS(lhsDesign(n_init, nvar, seed = 42)$design)$design
##' Z <- apply(X, 1, ftest)
##'
##' points(X, Z)
##'
##' model <- model_init <- mleHetGP(X = X, Z = Z, lower = lower, upper = upper)
##'
##' for(ii in 1:(N_tot - n_init)){
##'  ##Precalculations
##'  Wijs <- hetGP:::Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
##'  
##'  ## Adapt the horizon based on the training rmspe/score
##'    current_horizon <- horizon(model = model, Wijs = Wijs)
##'
##'  if(current_horizon == -1){
##'    opt <- IMSPE_optim(model = model, h = 0, Wijs = Wijs)
##'  }else{
##'    opt <- IMSPE_optim(model = model, h = current_horizon, Wijs = Wijs)
##'  }
##'  
##'  Xnew <- opt$par
##'  Znew <- apply(Xnew, 1, ftest)
##'  X <- rbind(X, Xnew)
##'  Z <- c(Z, Znew)
##'  points(Xnew, Znew)
##'  
##'  ## Update of the model
##'  model <- update(object = model, Xnew = Xnew, Znew = Znew, lower = lower, upper = upper)
##'  if(ii %% 25 == 0 || ii == (N_tot - n_init)){
##'    model_test <- mleHetGP(X = list(X0 = model$X0, Z0 = model$Z0, mult = model$mult), Z = model$Z,
##'                         lower = lower, upper = upper, maxit = 1000)
##'    model <- compareGP(model, model_test)
##'  }
##'}
##' ### Plot result
##' preds <- predict(x = Xgrid, model)
##' lines(Xgrid, preds$mean, col = 'red', lwd = 2)
##' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
##' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
##' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
##' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
##' par(new = TRUE)
##' plot(NA,NA, xlim = c(0, 1), ylim = c(0,max(model$mult)), axes = FALSE, ylab = "", xlab = "")
##' segments(x0 = model$X, x1 = model$X, y0 = rep(0, nrow(model$X)), y1 = model$mult, col = 'grey')
##' axis(side = 4)
##' mtext(side = 4, line = 2, expression(a[i]), cex = 0.8)
##' mtext(side = 2, line = 2, expression(f(x)), cex = 0.8)
##' mtext(side = 1, line = 2, 'x', cex = 0.8)
##' }
NULL