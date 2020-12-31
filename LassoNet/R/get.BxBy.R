#' Computes decomposition elements, p.g. 14 of Weber et. al (2014).
#' 
#' @param x Input data matrix of size \eqn{n \times p}, n - number of observations, p - number of covariates
#' @param y Response vector or size \eqn{n \times 1}
#' @param M Penalty matrix
#' @return Bx array of \eqn{B^{ij}_X} stored matrices. \eqn{Bx[,,k]} are the k-th combination of i and j non zero entry in the penalty matrix M
#' @return By array of \eqn{B^{ij}_y} stored matrices. \eqn{By[,k]} are the k-th combination of i and j non zero entry in the penalty matrix M
#' @examples
#' p<-200
#' n<-100
#' x<-matrix(rnorm(n*p),n,p)
#' y<-rnorm(n,mean=0,sd=1)
#' M<-diag(p)
#' get.BxBy(x, y, M)
get.BxBy   <- function(x , y ,M){
  # check only lower triangular values, since M is symmetric
  M[upper.tri(M)] <- 0
  non.zero <- which(M!=0,arr.ind = T)
  # in case diagonal has a non zero entry, we exlude that
  non.diag <- non.zero[,1]!=non.zero[,2]
  # store indicies of lower triagular signs being not zero
  e.con    <- non.zero[non.diag,]
  # number of iterations needed 
  k        <- nrow(e.con)
  # number of covariates
  p        <- ncol(x)
  # initializing
  By       <- array(0, c(2,k))
  Bx       <- array(0, c(2,p-2,k))
  if (k!=0){ # if there are off diagonal elements (i.e. if we are not at the ridge case)
    # loop through all indicies
    for (d in 1:k){
      elem    <- e.con[d,]
      
      x.tilde <- cbind(x[,elem[1]],x[,elem[2]])
      x.exc   <- cbind(x[,-c(elem[1],elem[2])])
      
      Byij    <- fastols(y    , x.tilde)
      Bxij    <- fastols(x.exc, x.tilde)
      
      By[,d]  <- Byij
      Bx[,,d] <- Bxij
    }
  } else{  
    Bx <- 0
    By <- 0
  }
  return(list(Bx = Bx, By = By))
}