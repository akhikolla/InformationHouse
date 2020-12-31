#' Updates \eqn{\beta} coefficients.
#' 
#' @param x input data matrix of size \eqn{n \times p}, n - number of observations, p - number of covariates
#' @param y response vector or size \eqn{n \times 1}
#' @param beta initial value for \eqn{\beta}. default - zero vector of size \eqn{n \times 1}
#' @param lambda1 lasso penalty coefficient
#' @param lambda2 network penalty coefficient
#' @param M1 penalty matrix
#' @param n.iter maximum number of iterations for \eqn{\beta} updating. default - 1e5
#' @param iscpp binary choice for using cpp function in coordinate updates. 1 - use C++ (default), 0 - use R.
#' @param tol convergence tolerance level. default - 1e-6
#' @return beta updated \eqn{\beta} vector
#' @return convergence binary variable for convergence
#' @return steps number of steps until convergence
beta.update.net <- function(x ,y ,beta=array(0,length(y)) ,lambda1 ,lambda2 ,M1 ,n.iter = 1e5 ,iscpp=TRUE ,tol=1e-6 ){
  
  beta <- as.matrix(beta)
  x    <- as.matrix(x)
  M1    <- as.matrix(M1)
  
  # initializing the variables
  iii    <- length(beta)
  beta.previous           <- numeric(length(beta))
  k <- 1
  conv <- FALSE
  xy <- crossprod(x,y)
  xx <- crossprod(x)
  
  if (iscpp==TRUE){ # using cpp
    # access cpp function
    est   <- betanew_lasso_cpp(xx, xy, beta, M1, y, lambda1, lambda2, n.iter, tol)
    # store beta
    beta  <- est$beta
    # store steps until beta convergence
    k    <- est$steps
    # store value for convergence
    if(k < n.iter){
      conv <- TRUE
    }
    k <- k+1 # cpp starts at 0, thus the output gives one step less than R
  } else { # using R 
    
    while ((conv == FALSE) && (k-1 < n.iter) ) {
      # store previous beta
      beta.previous   <-  beta
      # update each at a time each coordinates
      for (j in 1:iii){
        #  sum_i=1^n [ 2*{n*beta.tilde(j) + <xj,y> - sum_k=1^p [ <xj,xk> beta.tilde(k) ]} - 2*lambda2 sum_j!v [ M(jv) beta.tilde(v) ]
        numerator   <- 2*((length(y)*beta[j] + xy[j] - xx[,j]%*%beta)) - 2*lambda2*M1[j,-j]%*%beta[-j]
        #  2*n + 2*lambda2 M(jj)
        denominator <- 2*length(y)+2*lambda2*M1[j,j]
        # update beta(j) = numerator/denominator
        beta[j] <- soft.thresh(numerator,lambda1)/(denominator)
      }
      #update beta
      conv <- sum(abs(beta-beta.previous))/iii<tol
      k <- k+1
    }
  }
  return(list(beta = beta, convergence = conv, steps = k-1))
}
