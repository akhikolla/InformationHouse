#' Predict methods
#'
#' Make a RasterLayer with predictions from a fitted model object.
#'
#' @param object model object
#' @param newdata optional new data
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Number of CPU cores to utilize for parallel processing
#' @param filename character. Optional filename to save the RasterBrick output
#'   to file. If this is not provided, a temporary file will be created for large
#'   \code{x}
#' @param ... Additional arguments for \code{\link[raster]{writeRaster}}
#'
#' @include CENFA.R vulnerability-class.R
#' @name predict
NULL

#' @rdname predict
setMethod("predict",
          signature(object = "cnfa"),
          function(object, newdata, filename = "", parallel = FALSE, n = 1, ...){

            x <- get(as.character(object@call$x))
            if (is(x, "GLcenfa")) y <- raster(x) else {
              if (is(x, "Raster")) y <- x}
            nm <- names(y)
            if (!missing(newdata)){
              if (!all.equal(nm, names(newdata))) stop("layer names of newdata do not match layer names of model")
            }
            if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
              center <- cellStats(y, 'mean', na.rm = TRUE)
              sd <- cellStats(y, 'sd', na.rm = TRUE)
              if (missing(newdata) && is(x, "Raster")) y <- parScale(y, center = center, scale = sd, parallel = parallel, n = n)
            }

            if (!missing(newdata)) {
              if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
                y <- parScale(newdata, center = center, scale = sd, parallel = parallel, n = n)
              } else {
                y <- newdata
              }
            }

            m <- object@mf
            s <- object@sf
            filename <- trim(filename)
            if (!canProcessInMemory(y) && filename == '') {
              filename <- rasterTmpFile()
            }

            f1 <- function(x) (abs(x - m) %*%  s) / length(s)
            if(parallel) {
              beginCluster(n)
              ras <- clusterR(y, fun = .calc, args = list(fun = f1, forceapply = T, names = "Sensitivity"), filename = filename, ...)
              endCluster()
            } else {
              ras <- .calc(y, fun = f1, forceapply = T, filename = filename, names = "Sensitivity", ...)
            }

            return(ras)
          }
)

#' @rdname predict
setMethod("predict",
          signature(object = "enfa"),
          function(object, newdata, filename = "", parallel = FALSE, n = 1, ...){

            x <- get(as.character(object@call$x))
            if (!is(x, "Raster")) x <- raster(x)
            nm <- names(x)
            if (!missing(newdata)){
              if (!all.equal(nm, names(newdata))) stop("layer names of newdata do not match layer names of model")
            }
            U <- object@co
            if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
              center <- cellStats(x, 'mean', na.rm = TRUE)
              sd <- cellStats(x, 'sd', na.rm = TRUE)
              if (missing(newdata)) x <- parScale(x, center = center, scale = sd, parallel = parallel, n = n)
            }

            if (!missing(newdata)) {
              if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
                x <- parScale(newdata, center = center, scale = sd, parallel = parallel, n = n)
              } else {
                x <- newdata
              }
            }

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            f1 <- function(x) x %*% U
            if(parallel) {
              beginCluster(n)
              ras <- clusterR(x, fun = .calc, args = list(fun = f1, forceapply = T, names = nm), filename = filename, ...)
              endCluster()
            } else {
              ras <- .calc(x, fun = f1, forceapply = T, filename = filename, names = nm, ...)
            }

            return(ras)
          }
)

#' @rdname predict
setMethod("predict",
          signature(object = "departure"),
          function(object, filename = "", parallel = FALSE, n = 1, ...){

            x <- get(as.character(object@call$x))
            if (is(x, "GLdeparture")) {
              x <- raster(x)
            } else if (is(x, "Raster")) {
              y <- get(as.character(object@call$y))
              if (is.null(object@call$center)) center <- TRUE else center <- as.logical(as.character(object@call$center))
              if (is.null(object@call$scale)) scale <- TRUE else scale <- as.logical(as.character(object@call$scale))
              gld <- GLdeparture(x = x, y = y, center = center, scale = scale, parallel = parallel, n = n)
              x <- raster(gld)
            }

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            d <- object@df
            f1 <- function(x) (x %*% d) #/ length(d)
            if(parallel) {
              beginCluster(n)
              ras <- clusterR(x, fun = .calc, args = list(fun = f1, forceapply = T, names = "Exposure"), filename = filename, ...)
              endCluster()
            } else {
              ras <- .calc(x, fun = f1, forceapply = T, filename = filename, names = "Exposure", ...)
            }

            return(ras)
          }
)

#' @rdname predict
setMethod("predict",
          signature(object = "vulnerability"),
          function(object, newdata, filename = "", parallel = FALSE, n = 1, ...){

            x <- get(as.character(object@call$cnfa))
            y <- get(as.character(object@call$dep))
            if (is.null(object@call$w)) {
              w <- c(1, 1)
            } else {
              w <- as.numeric(object@call$w)
            }

            if (is.null(object@call$method)) {
              method <- "geometric"
            } else {
              method <- as.character(object@call$method)
            }

            if (missing(newdata)) {
              s.map <- predict(x, parallel = parallel, n = n)
            } else {
              s.map <- predict(x, newdata = newdata, parallel = parallel, n = n)
            }
            e.map <- predict(y, parallel = parallel, n = n)

            if (method == "arithmetic") {
              f1 <- function(x,y) (x*w[1] + y*w[2]) / sum(w)
              ras <- overlay(s.map, e.map, fun = f1, filename = filename, ...)
            } else if (method == "geometric") {
              if(w[1] == w[2]) {
                w <- c(1, 1)
              } else {
                w <- w / sum(w)
              }
              f1 <- function(x,y) (x^w[1] * y^w[2])^(1 / sum(w))
              ras <- overlay(s.map, e.map, fun = f1, filename = filename, ...)

            }
            return(ras)
          }
)
