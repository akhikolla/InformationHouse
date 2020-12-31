#' names
#'
#' Names of raster layers.
#'
#' @description
#' Get or set the names of the layers of a Raster* object
#'
#' @param x Raster* object
#' @param value character (vector)
#'
#' @return Character
#'
#' @examples
#' names(climdat.hist)
#'
#' @name names
#'
#' @importFrom methods .hasSlot
# @importFrom raster nlayers
# @importClassesFrom raster RasterLayer RasterStack
# @importMethodsFrom raster raster
NULL

#' @rdname names
setMethod('names', signature(x='Raster'),
          function(x) {
            if (.hasSlot(x@data, 'names')) {
              ln <- x@data@names
            } else {
              ln <- x@layernames
            }
            ln <- ln[1:nlayers(x)]
            validNames(as.vector(ln))
          }
)

#' @rdname names
setMethod('names', signature(x='RasterStack'),
          function(x) {
            ln <- sapply(x@layers, function(i) i@data@names)
            ln <- ln[1:nlayers(x)]
            validNames(as.vector(ln))
          }
)

#' @rdname names
setMethod('names<-', signature(x='Raster'),
          function(x, value)  {
            nl <- nlayers(x)
            if (is.null(value)) {
              value <- rep('', nl)
            } else if (length(value) != nl) {
              stop('incorrect number of layer names')
            }
            value <- validNames(value)

            if (inherits(x, 'RasterStack')){

              x@layers <- sapply(1:nl, function(i){
                r <- x@layers[[i]]
                r@data@names <- value[i]
                r
              })

            } else {
              if (.hasSlot(x@data, 'names')) {
                x@data@names <- value
              } else {
                x@layernames <- value
              }
            }

            return(x)
          }
)

#' @keywords internal
.uniqueNames <- function(x, sep='.') {
  y <- as.matrix(table(x))
  y <- y[y[,1] > 1, ,drop=F]
  if (nrow(y) > 0) {
    y <- rownames(y)
    for (i in 1:length(y)) {
      j <- which(x==y[i])
      x[j] <- paste(x[j], sep, 1:length(j), sep='')
    }
  }
  x
}

#' @keywords internal
.goodNames <- function(ln, prefix='layer') {
  validNames(ln, prefix)
}

#' @keywords internal
validNames <- function(x, prefix='layer') {
  x <- trim(as.character(x))
  x[is.na(x)] <- ""
  if (.standardnames()) {
    x[x==''] <- prefix
    x <- make.names(x, unique=FALSE)
  }
  .uniqueNames(x)
}

#' @keywords internal
.standardnames <- function(..., standardnames) {
  if (missing(standardnames)) {
    standardnames <- getOption('rasterStandardNames')
    if (is.null(standardnames)) {
      return(TRUE)  # the default
    } else {
      try (todisk <- as.logical(standardnames))
      if (is.logical(standardnames)) {
        return(standardnames)
      } else {
        return(TRUE)
      }
    }
  } else {
    if (is.logical(todisk)) {
      return(todisk)
    } else {
      return(TRUE)
    }
  }
}



