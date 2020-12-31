#' GLdeparture-class
#'
#' An object of class \code{GLdeparture} is created by the \code{\link{GLdeparture}}
#' function. It is best used for making comparisons between
#' species in the same study area. It speeds up the computation of multiple
#' departures by calculating the global covariance matrix as a first step, which
#' can then be fed into the \code{\link{departure}} function as a first argument.
#' This saves the user from having to calculate the global covariance matrix for
#' each species, which can take quite a bit of time.
#'
#' @slot global_difras Raster* object of absolute differences between historical
#'   \code{x} and future \code{y} climate values
#' @slot cov matrix. Global covariance matrix
#' @export

setClass("GLdeparture", slots = list(global_difras = "Raster", cov = "matrix"))

setMethod ("show", "GLdeparture", function(object){
  if (!inherits(object, "GLdeparture"))
    stop("Object of class 'GLdeparture' expected")
  cat("GLdeparture")
  cat("\n")
  cat("\n")
  print(round(object@cov, 2))
}
)

