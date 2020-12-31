#' Coerce from double to integer if safe
#' @description The same as \code{as.integer(x)} but only if \code{x} consists only of
#' whole numbers and is within the range of integers.
#'
#' @param x A double vector. If not a double vector, it is simply returned without any coercion.
#'
#' @examples
#'
#' N <- 1e6  # run with 1e9
#' x <- rep_len(as.double(sample.int(100)), N)
#' alt_as_integer <- function(x) {
#'   xi <- as.integer(x)
#'   if (isTRUE(all.equal(x, xi))) {
#'     xi
#'   } else {
#'     x
#'   }
#' }
#' bench_system_time(as_integer_if_safe(x))
#' #> process    real
#' #>  6.453s  6.452s
#' bench_system_time(alt_as_integer(x))
#' #> process    real
#' #> 15.516s 15.545s
#' bench_system_time(as.integer(x))
#' #> process    real
#' #>  2.469s  2.455s
#'
#' @export

as_integer_if_safe <- function(x) {
  if (is.double(x)) {
    if (is_altrep(x)) {
      # Implies out of range
      return(x)
    }
    status <- is_safe2int(x)
    if (status) {
      return(force_as_integer(x, status))
    }
  }
  x
}


