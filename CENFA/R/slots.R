#' Accessing CENFA slots
#'
#' Functions for extracting data from slots of objects of classes \code{cnfa} and \code{enfa}.
#'
#' @aliases m.factor s.factor marginality specialization sensitivity cov raster
#'
#' @param x cnfa or enfa object
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' m.factor(mod1)
#'
#' @name slot-access
NULL

#' @rdname slot-access
#' @export
m.factor <- function(x){
  return(x@mf)
}

#' @rdname slot-access
#' @export
s.factor <- function(x){
  return(x@sf)
}

#' @rdname slot-access
#' @export
marginality <- function(x){
  if (!inherits(x, c("cnfa", "enfa"))) stop("Object of class 'cnfa' or 'enfa' expected")
  return(x@marginality)
}

#' @rdname slot-access
#' @export
specialization <- function(x){
  if (!inherits(x, "enfa"))
    stop("Object of class 'enfa' expected")
  return(x@specialization)
}

#' @rdname slot-access
#' @export
sensitivity <- function(x){
  if (!inherits(x, "cnfa"))
    stop("Object of class 'cnfa' expected")
  return(x@sensitivity)
}

# "cov" <- function(x, ...) UseMethod("cov")

#' @rdname slot-access
#' @export
cov.cnfa <- function(x){
  if (!inherits(x, "cnfa"))
    stop("Object of class 'cnfa' expected")
      x@cov
}

#' @rdname slot-access
#' @export
cov.enfa <- function(x){
  if (!inherits(x, "enfa"))
    stop("Object of class 'enfa' expected")
  x@cov
}

#' @rdname slot-access
#' @export
cov.GLcenfa <- function(x){
  if (!inherits(x, "GLcenfa"))
    stop("Object of class 'GLcenfa' expected")
  x@cov
}

#' @rdname slot-access
setMethod("raster",
          signature(x = "enfa"),
          function(x){
            temp <- x@ras
            return(temp)
          }
)

#' @rdname slot-access
setMethod("raster",
          signature(x = "cnfa"),
          function(x){
            temp <- x@ras
            return(temp)
          }
)

#' @rdname slot-access
setMethod("raster",
          signature(x = "GLcenfa"),
          function(x){
            temp <- x@global_ras
            return(temp)
          }
)

#' @rdname slot-access
setMethod("raster",
          signature(x = "GLdeparture"),
          function(x){
            temp <- x@global_difras
            return(temp)
          }
)

#' @rdname slot-access
setMethod("names",
          signature(x = "GLcenfa"),
          function(x){
            temp <- names(x@global_ras)
            return(temp)
          }
)

#' @rdname slot-access
setMethod("names",
          signature(x = "GLdeparture"),
          function(x){
            temp <- names(x@global_difras)
            return(temp)
          }
)

#' @rdname slot-access
setMethod("names",
          signature(x = "cnfa"),
          function(x){
            temp <- names(x@ras)
            return(temp)
          }
)

#' @rdname slot-access
setMethod("names",
          signature(x = "enfa"),
          function(x){
            temp <- names(x@ras)
            return(temp)
          }
)

#' @rdname slot-access
setMethod("names",
          signature(x = "departure"),
          function(x){
            temp <- names(x@ras)
            return(temp)
          }
)


# comment out
# #' @rdname covmat
# #' @export
# setMethod("covmat",
#           signature(x = "enfa"),
#           function(x) x@cov
# )
#
# #' @rdname covmat
# #' @export
# setMethod("covmat",
#           signature(x = "GLcenfa"),
#           function(x) x@cov
# )
