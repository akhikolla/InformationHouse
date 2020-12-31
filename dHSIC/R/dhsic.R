dhsic <- function(X, Y, K, kernel = "gaussian", bandwidth=1, matrix.input=FALSE){
# Copyright (c) 2010 - 2017  Jonas Peters  [peters@stat.math.ethz.ch]
#                            Niklas Pfister [pfisteni@student.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

  # outputs the empirical dHSIC

  ###
  # Prerequisites
  ###

  # Check whether gram matrices have been specified in K and compute them if not
  if(missing(K)){
    # Check if Y has been passed as argument (implies two variable case)
    if(!missing(Y)){
      X <- list(X,Y)
    }
    
    # Deal with matrix input
    if(matrix.input){
      if(!is.matrix(X)){
        stop("X needs to be a matrix if matrix.input=TRUE")
      }
      else{
        X <- split(X, rep(1:ncol(X), each = nrow(X)))
      }
    }    
      
    # Set d, convert to matricies and set len
    d <- length(X)
    for(j in 1:d){
      if(is.matrix(X[[j]])==FALSE){
        X[[j]] <- as.matrix(X[[j]])
      }
    }     
    len <- nrow(X[[1]])
    
    # Case: len<2d
    if(len<2*d){
      warning("Sample size is smaller than twice the number of variables. dHSIC is trivial.")
      result=list(dHSIC = 0,
                  time = c(GramMat=0,HSIC=0))
      return(result)
    }
  
    
    # Check if enough kernels where set given the number of variables (else only used first kernel)
    if(length(kernel)<d){
      kernel <- rep(kernel[1],d)
    }
    if(length(bandwidth)<d){
      bandwidth <- rep(bandwidth[1],d)
    }
  
    # Define median heuristic bandwidth function
    median_bandwidth <- function(x){
      bandwidth <- median_bandwidth_rcpp(x[sample(1:len),,drop=FALSE],len,ncol(x))
      if(bandwidth==0){
        bandwidth <- 0.001
      }
      return(bandwidth)
    }
  
    # Define custom kernel Gram-matrix
    custom_grammat <- function(x,fun){
      KX <- matrix(nrow=len,ncol=len)
      for (i in 1:len){
        for (j in i:len){
          KX[i,j] <- match.fun(fun)(x[i,],x[j,])
          KX[j,i] <- KX[i,j]
        }
      }
      return(KX)
    }
  
    ###
    # Compute Gram-matrix (using specified kernels)
    ###
    K <- vector("list", d)
    ptm <- proc.time()
    for(j in 1:d){
      if(kernel[j]=="gaussian"){
        bandwidth[j] <- median_bandwidth(X[[j]])
        K[[j]] <- gaussian_grammat_rcpp(X[[j]],bandwidth[j],len,ncol(X[[j]]))
      }
      else if(kernel[j]=="gaussian.fixed"){
        K[[j]] <- gaussian_grammat_rcpp(X[[j]],bandwidth[j],len,ncol(X[[j]]))
      }
      else if(kernel[j]=="discrete"){
        bandwidth[j] <- NA
        K[[j]] <- discrete_grammat_rcpp(X[[j]],len,ncol(X[[j]]))
      }
      else{
        bandwidth[j] <- NA
        K[[j]] <- custom_grammat(X[[j]],kernel[j])
      }
    }
    timeGramMat <- as.numeric((proc.time() - ptm)[1])
  }
  else{
    # check if K is a list
    if(is.list(K)){
      d <- length(K)
      len <- nrow(K[[1]])
      timeGramMat <- NA
    }
    else{
      stop("K needs to be a list of matrices")
    }
  }
  
  ###
  # Compute dHSIC
  ###
  ptm <- proc.time()
  term1 <- 1
  term2 <- 1
  term3 <- 2/len
  for (j in 1:d){
    term1 <- term1*K[[j]]
    term2 <- 1/len^2*term2*sum(K[[j]])
    term3 <- 1/len*term3*colSums(K[[j]])
  }
  term1 <- sum(term1)
  term3 <- sum(term3)
  dHSIC=1/len^2*term1+term2-term3
  timeHSIC <- as.numeric((proc.time() - ptm)[1])

  ###
  # Collect result
  ###
  result=list(dHSIC = dHSIC,
              time = c(GramMat=timeGramMat,HSIC=timeHSIC),
              bandwidth=bandwidth)
  
  return(result)

}
