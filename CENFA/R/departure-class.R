#' departure-class
#'
#' An object of class \code{departure} is created by the \code{\link{departure}}
#' departure function, which quantifies the amount of change between historical
#' and future climate conditions inside a species' habitat.
#'
#' @slot call Original function call
#' @slot df Departure factor
#' @slot departure Magnitude of the departure factor
#' @slot g.cov historical global covariance matrix
#' @slot ras Raster* object of transformed climate values
#' @slot weights Raster layer of weights used for departure calculation
#' @export

setClass("departure", slots = list(call = "call", df = "numeric", departure = "numeric", g.cov = "matrix", ras = "Raster", weights = "Raster"))

setMethod("show",
          signature = "departure",
          function(object){
            if (!inherits(object, "departure"))
              stop("Object of class 'departure' expected")
            cat("CLIMATIC DEPARTURE")
            cat("\n\nDeparture factor: \n")
            print(round(object@df, 2))
            cat("\nOverall departure: ")
            cat(signif(object@departure, 4))
          }
)
