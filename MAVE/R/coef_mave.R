#' Directions of CS or CMS of given dimension
#'
#' This function returns the basis matrix of CS or CMS of given dimension
#' @aliases coef.mave
#' @aliases coef.mave.dim
#' @rdname coef.mave
#' @param object the output of \code{\link{mave}} or the output of \code{\link{mave.dim}}
#' @param dim the dimension of CS or CMS. The value of dim should be given when the class of the
#' argument dr is mave. When the class of the argument dr is mave.dim and dim is not given, the
#' function will return the basis matrix of CS or CMS of dimension selected by \code{\link{mave.dim}}.
#' Note that the dimension should be > 0.
#' @param ... no use.
#'
#' @return dir the matrix of CS or CMS of given dimension
#' @seealso \code{\link{mave.data}} for obtaining the reduced data
#' @export
#' @method coef mave
#' @method coef mave.dim
#' @examples
#' x <- matrix(rnorm(400),100,4)
#' y <- x[,1]+x[,2]+as.matrix(rnorm(100))
#' dr <- mave(y~x)
#' dir3 <- coef(dr,3)
#'
#' dr.dim <- mave.dim(dr)
#' dir3 <- coef(dr.dim,3)
#' dir.best <- coef(dr.dim)
#'

coef.mave<-function(object, dim,...){
  return(as.matrix(object$dir[[dim]]))
}

#' @rdname coef.mave
#' @export

coef.mave.dim <- function(object, dim='dim.min',...){
  if(dim=='dim.min'){
    dim = which.min(object$cv)
  }
  return(coef.mave(object,dim))
}
