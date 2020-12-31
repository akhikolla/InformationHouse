#' enfa-class
#'
#' An object of class \code{enfa} is created from performing ecological-niche
#' factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot call Original function call
#' @slot mf numeric. Named vector representing the marginality factor, describing
#'   the location of the species niche relative to the global niche
#' @slot marginality numeric. Magnitude of the marginality factor \code{mf},
#' scaled by the global covariance matrix
#' @slot sf numeric. Named vector representing the specialization factor,
#'   equivalent to the eigenvalues of specialization
#' @slot specialization numeric. The square root of the sum of eigenvalues, divided
#'   by the length of \code{sf}
#' @slot sf.prop numeric. Named vector representing the proportion of
#'   specialization found on each factor
#' @slot co p x p matrix of standardized variable loadings
#' @slot cov p x p species covariance matrix
#' @slot ras RasterBrick of transformed climate values, with p layers
#' @slot weights Raster layer of weights used for ENFA calculation
#' @export

setClass("enfa", slots = list(call = "call", mf = "numeric", marginality = "numeric",
                              sf = "numeric", specialization = "numeric", sf.prop = "numeric",
                              co = "matrix", cov = "matrix", ras = "Raster", weights = "Raster"))

setMethod ("show", "enfa", function(object){
  if (!inherits(object, "enfa"))
    stop("Object of class 'enfa' expected")
  cat("ENFA\n")
  cat("\nOriginal function call: ")
  print(object@call)
  cat("\nMarginality factor: \n")
  print(round(object@mf, 2))
  cat("\nEigenvalues of specialization: \n")
  print(round(object@sf, 2))
  cat("\nPercentage of specialization contained in ENFA factors: \n")
  print(round(100*object@sf.prop, 2))
  cat("\nOverall marginality: ", round(object@marginality, 3), "\n")
  cat("\nOverall specialization: ", round(object@specialization, 3), "\n")
  cat("\nSignificant ENFA factors: \n")
  n <- brStick(object@sf.prop[-1])
  co <- as.data.frame(object@co)
  #co <- as.data.frame(object@co[order(abs(object@co[,1]), decreasing = T), ])
  print(round(co[, 1:(n+1)], 2))
}
)
