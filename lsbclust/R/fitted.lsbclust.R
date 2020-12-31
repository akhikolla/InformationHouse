#' Extract Fitted Values for LSBCLUST
#' 
#' An S3 method for \code{\link{fitted}} for class \code{"lsbclust"}.
#' 
#' @param object An object of class \code{"lsbclust"}
#' @param \dots Unimplemented
#' 
#' @return An array approximating the original data
#' @seealso \code{\link{lsbclust}}
#' 
#' @method fitted lsbclust
#' @export
fitted.lsbclust <- function(object, ...) {
  return(object$fitted)  
}