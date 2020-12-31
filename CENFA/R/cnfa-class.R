#' cnfa-class
#'
#' An object of class \code{cnfa} is created by performing climate-niche factor
#' analysis on species presence data using the \code{cnfa} function.
#'
#' @slot call Original function call
#' @slot mf numeric. Named vector representing the marginality factor, describing
#'   the location of the species niche relative to the global niche
#' @slot marginality numeric. Magnitude of the marginality factor \code{mf}, scaled
#'   by the global covariance matrix
#' @slot sf numeric. Named vector representing the sensitivity factor
#' @slot sensitivity numeric. The magnitude of the sensitivity factor \code{sf},
#'   scaled by the global covariance matrix
#' @slot eig numeric. Named vector representing the eigenvalues of specialization,
#'   reflecting the amount of variance on each factor
#' @slot co p x p matrix of standardized variable loadings
#' @slot cov p x p species covariance matrix
#' @slot g.cov p x p global covariance matrix
#' @slot ras RasterBrick of transformed climate values, with p layers
#' @slot weights Raster layer of weights used for CNFA calculation
#'
#' @export

setClass("cnfa", slots = list(call = "call", mf = "numeric", marginality = "numeric", sf = "numeric",
                           sensitivity = "numeric", eig = "numeric", co = "matrix", cov = "matrix",
                           g.cov = "matrix", ras = "Raster", weights = "Raster"))

setMethod ("show", "cnfa", function(object){
  if (!inherits(object, "cnfa"))
    stop("Object of class 'cnfa' expected")
  cat("CNFA\n")
  cat("\nOriginal function call: ")
  print(object@call)
  cat("\nMarginality factor: \n")
  print(round(object@mf, 2))
  cat("\nSensitivity factor: \n")
  print(round(object@sf, 2))
  cat("\nPercentage of specialization contained in CNFA factors: \n")
  p.spec <- 100 * object@eig / sum(object@eig)
  print(round(p.spec, 2))
  cat("\nOverall marginality: ", round(object@marginality, 3), "\n")
  cat("\nOverall sensitivity: ", round(object@sensitivity, 3), "\n")
  cat("\nSignificant CNFA factors: \n")
  n <- brStick(object@eig[-1])
  co <- as.data.frame(object@co)
  #co <- as.data.frame(object@co[order(abs(object@co[,1]), decreasing = T), ])
  print(round(co[, 1:(n+1)], 2))
}
)
