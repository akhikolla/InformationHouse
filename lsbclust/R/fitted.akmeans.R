#' Extract Fitted Values for akmeans
#' 
#' An S3 method for \code{\link{fitted}} for class \code{"akmeans"}.
#' 
#' @param object An object of class \code{"akmeans"}
#' @param \dots Unimplemented
#' 
#' @return An array approximating the original data
#' @seealso \code{\link{akmeans}}
#' 
#' @method fitted akmeans
#' @export
fitted.akmeans <- function(object, ...) {
  return(object$fitted)  
}
