#' Estimates coefficients and connection signs over the grid values (\eqn{\lambda}1, \eqn{\lambda}2) where connection signs are kept unchanged.
#'
#' @param x \eqn{n \times p} input data matrix
#' @param y response vector or size  \eqn{n \times 1}
#' @param beta.0 initial value for \eqn{\beta}. default - zero vector of size \eqn{n \times 1}
#' @param lambda1 lasso penalty coefficient
#' @param lambda2 network penalty coefficient
#' @param M1 penalty matrix
#' @param n.iter maximum number of iterations for \eqn{\beta} updating. default - 1e5
#' @param iscpp binary choice for using cpp function in coordinate updates. 1 - use C++ (default), 0 - use R.
#' @param tol convergence in \eqn{\beta} tolerance level. default - 1e-6
#' @return beta matrix of \eqn{\beta} coefficients
#' @return mse Mean squared error value
#' @return update.steps matrix with stored number of steps for \eqn{\beta} updates to converge 
#' @return convergence.in.grid matrix with stored values for convergence in \eqn{\beta} coefficients
lasso.net.fixed <- function(x ,y ,beta.0 = array(0,length(y)) ,lambda1=c(0,1) ,lambda2=c(0,1) ,M1 ,n.iter = 1e5 ,iscpp=TRUE ,tol=1e-6) {
  # initializing beta
  beta <- beta.0
  
  n1 <- length(lambda1)
  n2 <- length(lambda2)
  # convert to relevant data structure
  beta   <- as.matrix(beta)
  x      <- as.matrix(x)
  M1     <- as.matrix(M1)
  
  # storage
  beta.store  <- array(0,c(length(beta), n2, n1))   
  mse         <- array(0,c(n2,n1))
  n.store     <- array(0,c(n2,n1))
  conv.store2 <- array(0,c(n2,n1))
  # storage for warm starts
  beta.l1w <- numeric(length(beta)) # warm start for lambda1 
  beta.l2w <- numeric(length(beta)) # warm start for lambda2
  g        <- 1
  
  for (j in 1:n1) {
    # loop through lambda1
    
    #message("Looping through lambda1 - ", j, "\n") 
    for (i in 1:n2) {
      
      # loop through lambda2
      #message("Looping through lambda2 - ", i, "\n") 
      
      betaconv <- FALSE
      
      
      while(betaconv!=TRUE){
        est <- beta.update.net(x,y,beta,lambda1[j],lambda2[i],M1,n.iter,iscpp)
        # store coordinate descent convergence output
        if(est$convergence == TRUE) {
          betaconv <- TRUE
        } 
        beta <- est$beta
        # store beta
        beta.store[,i,j]  <- beta
        # compute MSE 
        mse[i,j]  <- mean((y - x%*%beta)^2)
        # store number of steps until convergence
        n.store[i,j]     <- est$steps
        if (j == 1){
          beta.l1w <- est$beta
        }
        if(betaconv == FALSE) { 
          #message("coordinate descent did not converge", beta, "/n")
        }
        if(betaconv == FALSE) { 
          conv.store2[i,j] <- FALSE 
        }
        else {
          conv.store2[i,j] <- betaconv 
        }  
      } 
    } # end of loop (n2)
    beta  <- beta.l1w # warm start for lasso
  }  # end of loop (n1)
  
  return(list(beta = beta.store, mse = mse, update.steps = n.store, convergence.in.grid = conv.store2))
  
}
