#' Rank Based Fixed Effect Regression
#' 
#' Computes rank based regression estimates for fixed effect models.
#' 
#' 
#' @param f A model formula
#' @param data Data to use for model fitting
#' @param se Boolean indicating whether or not to calculate standard errors for
#' intercept and slope estimates
#' @param method Optimization method to use. Will accept any method usable by
#' optim, e.g. one of c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#' "Brent"). "BFGS" or "L-BFGS-B" are reccomended. "L-BFGS-B" should be used
#' for large datasets to conserve memory.
#' @return
#'   \item{fixed.effects}{ Fixed effect estimates }
#'   \item{ehat}{Residuals from model }
#'   
#' @author Herb Susmann
#' 
#' @seealso rlme, optim
#' 
#' @examples
#' 
#' # load schools data
#' data(schools)
#' 
#' # Fit fixed effects model with lmr
#' lmr.fit = lmr(y ~ age + sex, data=schools)
#' 
#' summary(lmr.fit)
#' 
#' # Fit with lmr and calculate standard errors
#' lmr.fit = lmr(y ~ age + sex, data=schools, se=TRUE)
#' 
#' summary(lmr.fit)
#' 
#' @export
lmr <- function(f, data, se=FALSE, method='L-BFGS-B') {
  #for large data, use method='L-BFGS-B'
  
  response_name = as.character(as.formula(f)[[2]])
  covariate_names = attributes(terms(f))$term.labels
  
  x = as.matrix(data[covariate_names])
  x = apply(x, 2, function(x) {
    x - mean(x)
  })
  
  y = as.matrix(data[[response_name]])
  
  # Fit with minimize_dispersion
  
  estimate = minimize_dispersion(x, y, se = se, method=method)
  #for large data, use method='L-BFGS-B'
  
  fit = list()
  fit$num.obs = nrow(y)
  fit$y = y
  fit$formula = f
  fit$method = "lmr"
  fit$ehat = estimate$ehat
  
  if(se == TRUE) {
    intercept.se = taustar(estimate$ehat, ncol(x)) / sqrt(nrow(x))
    
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), 
                                   Estimate = estimate$theta,
                                   StdError = c(intercept.se, estimate$standard.error))
  } else {
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), 
                                   Estimate = estimate$theta)
  }
  
  class(fit) = "rlme"
  
  return(fit)
}