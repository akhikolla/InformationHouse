## usethis namespace: start
#' @useDynLib cort, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("cort", libpath)
}
