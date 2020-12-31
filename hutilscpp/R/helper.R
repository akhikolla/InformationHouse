#' Helper
#' @param expr An expression
#' @return The expression evaluated.
#' @export
#'
#' @examples
#' x6 <- 1:6
#' helper(x6 + 1)
#'

# nocov start
helper <- function(expr) {
  sexpr <- substitute(expr)
  if (is.call(sexpr)) {
    cat("ok\n")
  }
  lhs <- sexpr[[2L]]
  cat(as.character(lhs), "\n")
  eval.parent(lhs)
}
# nocov end





