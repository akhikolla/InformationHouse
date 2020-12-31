#' Estimates coefficients and connection signs over the grid values (\eqn{\lambda}1, \eqn{\lambda}2).
#' 
#' @param x \eqn{n \times p} input data matrix
#' @param y response vector or size  \eqn{n \times 1}
#' @param beta.0 initial value for \eqn{\beta}. default - zero vector of size \eqn{n \times 1}
#' @param lambda1 lasso penalty coefficient
#' @param lambda2 network penalty coefficient
#' @param M1 penalty matrix
#' @param m.iter maximum number of iterations for sign matrix updating. default - 100
#' @param n.iter maximum number of iterations for \eqn{\beta} updating. default - 1e5
#' @param iscpp binary choice for using cpp function in coordinate updates. 1 - use C++ (default), 0 - use R.
#' @param tol convergence in \eqn{\beta} tolerance level. default - 1e-6
#' @param alt.num In case connection sign matrix alternates, alt.num last iterataions are stored default - 12 
#' @return beta matrix of \eqn{\beta} coefficients
#' @return mse Mean squared error value
#' @return M Array of connection signs
#' @return iterations matrix with stored number of steps for sign matrix to converge
#' @return update.steps matrix with stored number of steps for \eqn{\beta} updates to converge 
#' @return convergence.in.M matrix with stored values for convergence in sign matrix
#' @return convergence.in.grid matrix with stored values for convergence in \eqn{\beta} coefficients
#' @return xi.conv array with stored values about changes in connection signs
#' @return beta.alt vector for \eqn{\beta} value in the case when M1 matrix alternates, zero otherwise
lasso.net.grid <- function(x ,y ,beta.0 = array(0,length(y)) ,lambda1=c(0,1) ,lambda2=c(0,1) ,M1 ,m.iter = 100, n.iter = 1e5 , iscpp=TRUE, tol=1e-6,  alt.num=12){
  # initializing beta
  beta <- beta.0
  
  # define penalty grid space
  n1 <- length(lambda1)
  n2 <- length(lambda2)
  # convert to relevant data structure
  beta    <- as.matrix(beta)
  x       <- as.matrix(x)
  M1      <- as.matrix(M1)
  
  # initialize values
  signs   <- FALSE
  xi      <- sign(t(x)%*%x)
  M       <- matrix.M.update(M1, xi)
  k       <- 0
  alt.id  <- 0
  
  
  # storage
  beta.store  <- array(0,c(length(beta), n2, n1))   
  M.storej     <- NULL
  M.store      <- NULL
  k.store     <- array(0,c(n2,n1))
  n.store     <- array(0,c(n2,n1))
  conv.store  <- array(0,c(n2,n1)) 
  conv.store2 <- array(0,c(n2,n1))
  alt.store   <- array(0,c(dim(M1),alt.num))
  alt.beta    <- array(0,c(length(beta),alt.num))
  alt.xi.store<- array(0,c(dim(xi),alt.num))
  mse         <- array(0,c(n2,n1))
  
  alt.save    <- NULL
  xi.conv     <- NULL
  xi.conv1    <- NULL
  xi.conv2    <- NULL
  
  
  # storage for warm starts
  beta.l1w <- numeric(length(beta)) # warm start for lambda1 
  beta.l2w <- numeric(length(beta)) # warm start for lambda2
  
  # store matrices Bx and By for xi update before loops 
  BxBy <- get.BxBy(x,y,M1)
  Bx   <- BxBy$Bx
  By   <- BxBy$By
  
  
  
  for (j in 1:n1) {
    # loop through lambda1
    
    #message("Looping through lambda1 - ", j, "\n") 
    for (i in 1:n2) {
      
      # loop through lambda2
      #message("Looping through lambda2 - ", i, "\n") 
      # store check value for beta convergence given updated M 
      betaconv <- FALSE
      # while loop for xi convegence
      while ((signs != TRUE)  && (k < m.iter) && (betaconv != TRUE)) {
        
        # estimates lasso with network given lambda1, lambda2 and M matrix
        est <- beta.update.net(x,y,beta,lambda1[j],lambda2[i],M,n.iter,iscpp)
        
        # update xi with "covariance" type updates. sec 2.2.
        xi.u <- get.xi(Bx,By,est$beta,xi,M1)
        
        # store coordinate descent convergence output
        if(est$convergence == TRUE) {
          betaconv <- TRUE
        } 
        
        # check for alternating signs if iterator is close to maximum number of iteration
        if(k >= m.iter-alt.num){
          id <- m.iter-k-alt.num+1
          alt.xi.store[,,id] <- xi
          alt.beta[,id]     <- est$beta
        } 
        if(k == m.iter) 
        {
          alt.id                   <- alt.id + 1 
          alt.save$beta[[alt.id]]  <- alt.beta
          alt.save$xi[[alt.id]]    <- alt.xi.store
          alt.save$lams[[alt.id]]  <- c(lambda1[j],lambda2[i])
        }
        
        # update iterations number
        k <- k + 1
        # check xi matrix with all previous matrices
        xi.conv[[k]] <- cbind(which(xi.u!=xi, arr.ind = T),as.matrix(xi.u[which(xi.u!=xi, arr.ind = T)]))
        # check for convergence or stop once interations reaches n.iter
        signs <- isTRUE(all.equal(xi.u,xi))
        
        # update xi
        xi <- xi.u
        # storing warm starts for beta lambda1 (the first converged beta)
        if (k == 1 && j == 1){
          beta.l1w <- est$beta
        }
        # updating beta for lambda2 warm start
        beta <- est$beta
        if(betaconv == FALSE) { 
          #message("coordinate descent did not converge", beta, "/n")
        }
      } # Iterations through sign Matrix finishes here. [either converged or signs alternate]
      
      
      xi.conv1[[i]] <- xi.conv
      
      k.store[i,j]     <- k 
      n.store[i,j]     <- est$steps
      conv.store[i,j]  <- signs
      
      if(betaconv == FALSE) { 
        conv.store2[i,j] <- FALSE 
      }
      else {
        conv.store2[i,j] <- betaconv 
      }
      
      k     <- 0
      signs <- FALSE
      
      
      
      beta.store[,i,j]  <- beta
      # compute MSE 
      mse[i,j]  <- mean((y - x%*%beta)^2)
      # store connections
      M.storej[[i]]   <- get.signs.M(M)
      # upate the penalty matrix 
      M   <- matrix.M.update(M1,xi)
      
    } # end of for loop n2  (through lambda2)
    beta <- beta.l1w # warm start for lasso
    M.store[[j]] <- M.storej
    # store xi check
    xi.conv2[[j]] <- xi.conv1
    # update beta for lambda1 warm start
    beta <- beta.l2w
  }   # end of for loop n1  (through lambda1)
  
  return(list(beta = beta.store, mse = mse, signMat = M.store, iterations = k.store, update.steps = n.store, convergence.in.M = conv.store, convergence.in.grid = conv.store2, xi.conv = xi.conv2, alt.info = alt.save))
}
