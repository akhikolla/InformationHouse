#' Evaluate time of computation
#' @description (Used for examples and tests)
#' @param expr Passed to \code{\link[bench]{system_time}}.
#' @export
#'

bench_system_time <- function(expr) {
  # bench::system_time only evaluates expr in the
  # parent environment (reasonably). So we copy
  # the names of everything in expr so that
  # bench::system_time can run.
  for (`***` in all.vars(substitute(expr))) {
    assign(`***`, value = get(`***`, envir = parent.frame()))
  }
  if (requireNamespace("bench", quietly = TRUE)) {
    bench::system_time(expr)
  }
}
