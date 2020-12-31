#' GLcenfa-class
#'
#' An S4 object of class \code{GLcenfa} is created using the \code{\link{GLcenfa}}
#' function on a Raster* object. It is best used for making comparisons between
#' species in the same study area. It speeds up the computation of multiple CNFAs
#' or ENFAs by calculating the global covariance matrix as a first step, which
#' can then be fed into the \code{\link{cnfa}} or \code{\link{enfa}} functions
#' as their first argument. This saves the user from having to calculate the
#' global covariance matrix for each species, which can take quite a bit of time.
#'
#' @slot global_ras Raster* object \code{x} with p layers
#' @slot cov matrix. Global p x p covariance matrix
# @slot center numeric. Layer means of \code{x} before scaling
# @slot sd numeric. Layer standard deviations of \code{x} before scaling
#' @export

setClass("GLcenfa", slots = list(global_ras = "Raster", cov = "matrix"))#, center = "numeric", sd = "numeric"))

setMethod ("show", "GLcenfa", function(object){
  if (!inherits(object, "GLcenfa"))
    stop("Object of class 'GLcenfa' expected")
  cat("GLcenfa")
  cat("\n")
  cat("\n")
  print(round(object@cov, 2))
}
)
