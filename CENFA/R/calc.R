#' @keywords internal
#'
# @importFrom raster brick trim canProcessInMemory setValues writeRaster pbStep
#   nlayers getValues writeStart blockSize pbCreate writeValues writeStop pbClose
#' @importFrom foreach %dopar% foreach
#' @importFrom stats cov na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar

# This function is the exact same as raster::calc, with the addition of a
# 'names' argument, which allows the user to specify a character vector of
# layer names for the Raster* object that gets created.
#
# Internal functions from the raster package are included here, as
# necessary for raster::calc to run.

.calc <- function(x, fun, filename='', na.rm, forcefun=FALSE, forceapply=FALSE, names, ...) {

  nl <- nlayers(x)

  test <- .calcTest(x[1:5], fun, na.rm, forcefun, forceapply)
  doapply <- test$doapply
  makemat <- test$makemat
  trans <- test$trans

  if (test$nlout == 1) {
    out <- raster(x)
  } else {
    out <- brick(x, values=FALSE)
    out@data@nlayers <- test$nlout
  }

  names(out) <- names

  fun <- .makeTextFun(fun)
  if (class(fun) == 'character') {
    doapply <- FALSE
    fun <- .getRowFun(fun)
  }

  filename <- trim(filename)

  if (canProcessInMemory(x, max(nlayers(x), nlayers(out)) * 2)) {
    x <- getValues(x)
    if (makemat) {
      x <- matrix(x, ncol=1)
    }
    if (missing(na.rm)) {
      if (! doapply ) {
        x <- fun(x )
      } else {
        x <- apply(x, 1, fun )
      }
    } else {
      if ( ! doapply ) {
        x <- fun(x, na.rm=na.rm )
      } else {
        x <- apply(x, 1, fun, na.rm=na.rm)
      }
    }
    if (trans) {
      x <- t(x)
    }
    x <- setValues(out, x)
    if (filename != '') {
      x <- writeRaster(x, filename, ...)
    }
    return(x)
  }

  out <- writeStart(out, filename=filename, ...)
  tr <- blockSize(out)
  pb <- pbCreate(tr$n, label='calc', ...)

  if (missing(na.rm)) {
    for (i in 1:tr$n) {
      v <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
      if ( ! doapply ) {
        v <- fun(v)
      } else {
        if (makemat) {
          v <- matrix(v, ncol=1)
        }
        v <- apply(v, 1, fun)
        if (trans) {
          v <- t(v)
        }
      }
      out <- writeValues(out, v, tr$row[i])
      pbStep(pb)
    }
  } else {
    for (i in 1:tr$n) {
      v <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
      if ( ! doapply ) {
        v <- fun(v, na.rm=na.rm)
      } else {
        if (makemat) {
          v <- matrix(v, ncol=1)
        }
        v <- apply(v, 1, fun, na.rm=na.rm)
        if (trans) {
          v <- t(v)
        }
      }
      out <- writeValues(out, v, tr$row[i])
      pbStep(pb)
    }
  }
  out <- writeStop(out)
  pbClose(pb)
  return(out)
}

.calcTest <- function (tstdat, fun, na.rm, forcefun = FALSE, forceapply = FALSE) {
  if (forcefun & forceapply) {
    forcefun <- FALSE
    forceapply <- FALSE
  }
  trans <- FALSE
  doapply <- FALSE
  makemat <- FALSE
  nl <- NCOL(tstdat)
  if (nl == 1) {
    if (forceapply) {
      doapply <- TRUE
      makemat <- TRUE
      tstdat <- matrix(tstdat, ncol = 1)
      if (missing(na.rm)) {
        test <- try(apply(tstdat, 1, fun), silent = TRUE)
      }
      else {
        test <- try(apply(tstdat, 1, fun, na.rm = na.rm),
                    silent = TRUE)
      }
      if (length(test) < length(tstdat) | class(test) ==
          "try-error") {
        stop("cannot forceapply this function")
      }
      if (is.matrix(test)) {
        if (ncol(test) > 1) {
          trans <- TRUE
        }
      }
    }
    else {
      if (!missing(na.rm)) {
        test <- try(fun(tstdat, na.rm = na.rm), silent = TRUE)
        if (class(test) == "try-error") {
          test <- try(apply(tstdat, 1, fun, na.rm = na.rm),
                      silent = TRUE)
          doapply <- TRUE
          if (class(test) == "try-error") {
            stop("cannot use this function. Perhaps add '...' or 'na.rm' to the function arguments?")
          }
          if (is.matrix(test)) {
            if (ncol(test) > 1) {
              trans <- TRUE
            }
          }
        }
      }
      else {
        test <- try(fun(tstdat), silent = TRUE)
        if (length(test) < length(tstdat) | class(test) ==
            "try-error") {
          doapply <- TRUE
          makemat <- TRUE
          tstdat <- matrix(tstdat, ncol = 1)
          test <- try(apply(tstdat, 1, fun), silent = TRUE)
          if (class(test) == "try-error") {
            stop("cannot use this function")
          }
          if (is.matrix(test)) {
            if (ncol(test) > 1) {
              trans <- TRUE
            }
          }
        }
      }
    }
  }
  else {
    if (forcefun) {
      doapply <- FALSE
      test <- fun(tstdat)
    }
    else {
      doapply <- TRUE
      if (!missing(na.rm)) {
        test <- try(apply(tstdat, 1, fun, na.rm = na.rm),
                    silent = TRUE)
        if (class(test) == "try-error") {
          doapply <- FALSE
          test <- try(fun(tstdat, na.rm = na.rm), silent = TRUE)
          if (class(test) == "try-error") {
            stop("cannot use this function. Perhaps add '...' or 'na.rm' to the function arguments?")
          }
        }
        else if (is.matrix(test)) {
          trans <- TRUE
        }
      }
      else {
        test <- try(apply(tstdat, 1, fun), silent = TRUE)
        if (class(test) == "try-error") {
          doapply <- FALSE
          test <- try(fun(tstdat), silent = TRUE)
          if (class(test) == "try-error") {
            stop("cannot use this function")
          }
        }
        else if (is.matrix(test)) {
          trans <- TRUE
        }
      }
    }
  }
  if (trans) {
    test <- t(test)
    test <- ncol(test)
  }
  else {
    test <- length(test)/5
  }
  nlout <- as.integer(test)
  list(doapply = doapply, makemat = makemat, trans = trans,
       nlout = nlout)
}

.getRowFun <- function (fun) {
  if (fun == "mean") {
    return(rowMeans)
  }
  else if (fun == "sum") {
    return(rowSums)
  }
  else if (fun == "min") {
    return(.rowMin)
  }
  else if (fun == "max") {
    return(.rowMax)
  }
  else {
    stop("unknown fun")
  }
}

.rowMin <- function (x, na.rm = TRUE) {
  .doRowMin(x, narm = na.rm)
}

.rowMax <- function (x, na.rm = TRUE) {
  .doRowMax(x, narm = na.rm)
}

.makeTextFun <- function (fun) {
  if (class(fun) != "character") {
    if (is.primitive(fun)) {
      test <- try(deparse(fun)[[1]], silent = TRUE)
      if (test == ".Primitive(\"sum\")") {
        fun <- "sum"
      }
      else if (test == ".Primitive(\"min\")") {
        fun <- "min"
      }
      else if (test == ".Primitive(\"max\")") {
        fun <- "max"
      }
    }
    else {
      test1 <- isTRUE(try(deparse(fun)[2] == "UseMethod(\"mean\")",
                          silent = TRUE))
      test2 <- isTRUE(try(fun@generic == "mean", silent = TRUE))
      if (test1 | test2) {
        fun <- "mean"
      }
    }
  }
  return(fun)
}

.fullFilename <- function (x, expand = FALSE) {
  x <- trim(x)
  if (identical(basename(x), x)) {
    x <- file.path(getwd(), x)
  }
  if (expand) {
    x <- path.expand(x)
  }
  return(x)
}
