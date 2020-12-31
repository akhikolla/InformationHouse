#' The reduced data matrix
#'
#' The function returns the reduced data matrix of the original data. The reduced data matrix
#' is obtained by the original data multiplied by the dimension reduction directions of given
#' dimension.
#'
#' @param dr the object returned by \code{\link{mave}} or \code{\link{mave.dim}}
#' @param x the original data matrix of p dimensions
#' @param dim the dimension of the reduced data matrix.
#' @seealso \code{\link{coef.mave}} for obtaining the dimension reduction directions
#' @export
#' @examples
#'
#' x <- matrix(rnorm(400),100,4)
#' y <- x[,1]+x[,2]+as.matrix(rnorm(100))
#' dr <- mave(y~x)
#' x.reduced <- mave.data(dr,x,3)

mave.data<-function(dr, x, dim = NULL){

  if(class(dr)=='mave'){
    if(!is.null(dim)){
      if(dim<=dr$max.dim){
         return(as.matrix(x%*%coef.mave(dr,dim)))
      }
      else{
        arg='central space'
        if(substr(dr$method,1,2)!='CS') arg='central mean space'
        stop(c('The ', arg, ' of dimension ',dim,' haven\'t been computed, please compute it first'))
      }
    }
    else{
      stop('dim should be given.')
    }
  }
  if(class(dr)=='mave.dim'){
    if(is.null(dim)){
      return(as.matrix(x%*%coef.mave.dim(dr)))
    }
    else{
      return(as.matrix(x%*%coef.mave(dr,dim)))
    }
  }

}
