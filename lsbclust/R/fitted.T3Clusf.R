#' Extract Fitted Values for T3Clusf
#' 
#' An S3 method for \code{\link{fitted}} for class \code{"T3Clusf"}.
#' 
#' @param object An object of class \code{"T3Clusf"}
#' @param \dots Unimplemented
#' 
#' @return An array approximating the original data
#' @seealso \code{\link{T3Clusf}}
#' 
#' @method fitted T3Clusf
#' @export
fitted.T3Clusf <- function(object, ...) {
  return(object$fitted)  
}
