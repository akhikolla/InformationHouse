#' Separated which
#' @description Same as \code{which(exprA)} where \code{exprA} is a binary
#' expression.
#'
#' @param exprA An expression. Useful when of the form \code{a <op> b} for
#' \code{a} an atomic vector. Long expressions are not supported.
#' @param .env The environment in which \code{exprA} is to be evaluated.
#' @param nThread Number of threads to use.
#'
#' @return Integer vector, the indices of \code{exprA} that return \code{TRUE}.
#'
#' @export

whichs <- function(exprA, .env = parent.frame(), nThread = getOption("hutilscpp.nThread", 1L)) {
  sexprA <- substitute(exprA)
  is_binary_sexpr <-
    is.call(sexprA) &&
    length(sexprA) == 3L

  if (!is_binary_sexpr) {
    return(which(exprA))
  }
  opc <- as.character(sexprA[[1L]])
  opm <- op2M(opc)
  lhs <- eval(sexprA[[2L]], envir = .env, enclos = baseenv())
  rhs <- eval(sexprA[[3L]], envir = .env, enclos = baseenv())
  rhs_good <-
    is.integer(rhs) ||
    AND(is.double(rhs) && length(rhs) == 1, rhs == as.integer(rhs))
  if (opm < 1 || opm > 6 || is.double(lhs)) {
    return(which(exprA))
  }
  do_whichs_16(opm,
               lhs,
               rhs,
               nThread = nThread)
}


