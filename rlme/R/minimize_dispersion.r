#' Minimize Dispersion Function
#' 
#' Uses optim to find regression estimates which minimize dispersion function
#' on X and Y input matrices
#' 
#' 
#' @param X Input matrix
#' @param Y Response vector
#' @param method Method optim should use - one of "Nelder-Mead", "BFGS", "CG",
#' "L-BFGS-B", "SANN", or "Brent".
#' @param init.guess How to calculate the first regression estimate. Defaults
#' to using quantile regression.
#' @param verbose Whether to print out verbose messages.
#' @param se Whether or not to calculate standard errors of regression
#' estimates.
#' 
#' @return
#'   \item{theta }{Regression parameter estimates}
#'   \item{ehat}{Regression residuals}
#'   
#' @author Herb Susmann
#' 
#' @importFrom quantreg rq
#' 
#' @export
#' 
minimize_dispersion = function(X, Y, method='BFGS', init.guess = 'quantreg', verbose=FALSE, se = TRUE) {
  #for large data, use method='L-BFGS-B'
  # install MASS and quantreg
  # se=FALSE if you want to get only R estimates
  
  rows = nrow(Y)
  
  model = Y ~ X
  
  if(init.guess == 'quantreg') {
    # Initial guess using quantreg
    beta.not = rq(model)$coefficients[-1] 
  }
  else {
    beta.not = rep(0, ncol(X))
  }
  
  if(verbose == TRUE) {
    cat("Minimize dispersion: Calling optim")
  }
  
  # Optimize the gradient function
  beta.optim = optim(beta.not, Y=Y, X=X, n=rows, D.beta, gr=S.beta, method=method)
  #for large data, use method='L-BFGS-B'
  
  if(verbose == TRUE) {
    cat("Minimize dispersion: optim complete")
  }
  
  beta = unname(beta.optim$par)
  
  # calculate residuals
  resd = calculate.residuals(X, Y, beta)
  
  inv = solve(t(X) %*% X)
  
  if(verbose == TRUE) {
    cat("Minimize dispersion: estimating tau")
  }
  
  # Calculate intercept
  intercept = median(resd)
  
  # Prepend intercept to beta list
  beta = c(intercept, beta)
  
  if(se == TRUE) {    
    # estimate tau
    tau = wilcoxontau(resd, ncol(X), verbose=verbose)
    
    # estimate tau star
    tau.star = taustar(resd, p = ncol(X), conf = 0.95)
    
    # calculate standard errors
    standard.error = sqrt((tau^2) * diag(inv))
    
    return(list(
      theta = beta,
      ehat = resd,
      tauhat = tau,
      taus = tau.star,
      standard.error = standard.error
    ))
  }

  return(list(
    theta = beta,
    ehat = resd
  ))
  
}

# Weights
a <- function(i, n) {
  sqrt(12) * (i / (n+1) - 0.5)
}

# Dispersion function given residuals
D.residuals <- function(e, n) {
  sum(a(rank(e), n) * e)
}

# Dispersion function given only beta
D.beta <- function(X, Y, beta, n) {
  e = calculate.residuals(X, Y, beta)
  D.residuals(e, n)
}

# Calculate residuals from estimated beta
calculate.residuals <- function(X, Y, beta) {
  Y - X %*% beta
}

# Dispersion gradient function given residuals
S <- function(X, e, n) {
  -1 * t(X) %*% a(rank(e), n)
}

# Dispersion gradient function given only beta
S.beta <- function(X, Y, beta, n) {
  e = calculate.residuals(X, Y, beta)
  S(X, e, n)
}
