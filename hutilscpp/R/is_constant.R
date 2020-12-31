#' @title Is a vector constant?
#' @name is_constant
#' @description Efficiently decide whether an atomic vector is constant; that is,
#' contains only one value.
#'
#' Equivalent to
#'
#' \code{data.table::uniqueN(x) == 1L}
#'
#' or
#'
#' \code{forecast::is.constant(x)}
#'
#' @param x An atomic vector. Only logical, integer, double, and character
#' vectors are supported. Others may work but have not been tested.
#'
#' @param nThread \describe{
#' \item{\code{integer(1)}}{Number of threads to use in \code{is_constant}.}
#' }
#'
#' @return
#' Whether or not the vector \code{x} is constant:
#' \describe{
#' \item{\code{is_constant}}{\code{TRUE} or \code{FALSE}. Missing values are considered to
#' be the same as each other, so a vector entirely composed of missing values is
#' considered constant. Note that \code{is_constant(c(NA_real_, NaN))} is \code{TRUE}.}
#' \item{\code{isntConstant}}{If constant, \code{0L}; otherwise, the first integer position at
#' which \code{x} has a different value to the first.
#'
#' This has the virtue of \code{!isntConstant(x) == is_constant(x)}.}
#' }
#'
#' Multithreaded \code{is_constant(x, nThread)} should only be used if
#' \code{x} is expected to be true. It will be faster when
#' \code{x} is constant but much slower otherwise.
#'
#' Empty vectors are constant, as are length-one vectors.
#'
#'
#' @examples
#' library(hutilscpp)
#' library(data.table)
#' N <- 1e9L
#' N <- 1e6  # to avoid long-running examples on CRAN
#'
#' ## Good-cases
#' nonconst <- c(integer(1e5), 13L, integer(N))
#' bench_system_time(uniqueN(nonconst) == 1L)
#' #> process    real
#' #> 15.734s  2.893s
#' bench_system_time(is_constant(nonconst))
#' #> process    real
#' #>   0.000   0.000
#' bench_system_time(isntConstant(nonconst))
#' #> process    real
#' #>   0.000   0.000
#'
#' ## Worst-cases
#' consti <- rep(13L, N)
#' bench_system_time(uniqueN(consti) == 1L)
#' #> process    real
#' #>  5.734s  1.202s
#' bench_system_time(is_constant(consti))
#' #>   process      real
#' #> 437.500ms 437.398ms
#' bench_system_time(isntConstant(consti))
#' #>   process      real
#' #> 437.500ms 434.109ms
#'
#' nonconsti <- c(consti, -1L)
#' bench_system_time(uniqueN(nonconsti) == 1L)
#' #> process    real
#' #> 17.812s  3.348s
#' bench_system_time(is_constant(nonconsti))
#' #>   process      real
#' #> 437.500ms 431.104ms
#' bench_system_time(isntConstant(consti))
#' #>   process      real
#' #> 484.375ms 487.588ms
#'
#' constc <- rep("a", N)
#' bench_system_time(uniqueN(constc) == 1L)
#' #> process    real
#' #> 11.141s  3.580s
#' bench_system_time(is_constant(constc))
#' #> process    real
#' #>  4.109s  4.098s
#'
#' nonconstc <- c(constc, "x")
#' bench_system_time(uniqueN(nonconstc) == 1L)
#' #> process    real
#' #> 22.656s  5.629s
#' bench_system_time(is_constant(nonconstc))
#' #> process    real
#' #>  5.906s  5.907s
#'
#'
#' @export is_constant isntConstant

is_constant <- function(x, nThread = getOption("hutilscpp.nThread", 1L)) {
  if (!is.atomic(x)) {
    stop("`x` was not atomic. ",
         "Such objects are not supported.")
  }
  if (is.double(x) && is.na(x[1])) {
    return(all_na_real(x, nThread = nThread))
  }
  do_is_constant(x, nThread = nThread)
}

#' @rdname is_constant
isntConstant <- function(x) {
  if (!is.atomic(x)) {
    stop("`x` was not atomic. ",
         "Such objects are not supported.")
  }
  if (length(x) <= 1L) {
    return(0L)
  }

  x1 <- x[1L]
  if (is.logical(x)) {
    wmin <- which.min(x)
    if (anyNA(x1)) {
      if (length(wmin)) {
        wmax <- which.max(x[seq_len(wmin)])
        return(min(wmax, wmin))
      } else {
        return(0L)
      }
    } else {
      if (wmin == 1L) {
        # still unknown whether constant or not
        wmax <- which.max(x)
        if (wmax == 1L) {
          return(0L)
        } else {
          return(wmax)
        }
      } else {
        return(wmin)
      }
    }
  }

  do_isntConstant(x)
}
