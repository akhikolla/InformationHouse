#' Sum of logical expressions
#' @name sum_and3s
#' @aliases sum_or3s
#' @param exprA,exprB,exprC,... Expressions of the form \code{x <op> y}.
#' with \code{<op>} one of the standard binary operators.
#' @param nThread \describe{
#' \item{\code{integer(1)}}{Number of threads to use.}
#' }
#' @param .env The environment in which the expressions are to be evaluated.
#'
#' @return
#' Equivalent to \code{sum(exprA & exprB & exprC)} or
#' \code{sum(exprA | exprB | exprC)} as desired.
#'
#' @export sum_and3s


sum_and3s <- function(exprA, exprB, exprC, ...,
                      nThread = getOption("hutilscpp.nThread", 1L),
                      .env = parent.frame()) {
  sexprA <- substitute(exprA)
  sexprB <- substitute(exprB)
  sexprC <- substitute(exprC)
  missingA <- missing(exprA)
  missingB <- missing(exprB)
  missingC <- missing(exprC)

  d <- decompose_expr(sexprA, sexprB, sexprC,
                      missingA, missingB, missingC,
                      .env = .env,
                      nThread = nThread)

  do_sum3s_par(
    d[[1]],
    d[[2]],
    d[[3]],
    d[[4]],
    d[[5]],
    d[[6]],
    d[[7]],
    d[[8]],
    d[[9]],
    d[[10]],
    d[[11]],
    d[[12]],
    d[[13]],
    d[[14]],
    d[[15]],
    d[[16]],
    d[[17]],
    d[[18]],
    d[[19]],
    d[[20]],
    d[[21]],
    d[[22]],
    d[[23]],
    d[[24]],
    ampersand = TRUE,
    nThread = nThread
  )
}

#' @rdname sum_and3s
#' @export

sum_or3s <- function(exprA, exprB, exprC, ...,
                     .env = parent.frame(),
                     nThread = getOption("hutilscpp.nThread", 1L)) {
  sexprA <- substitute(exprA)
  sexprB <- substitute(exprB)
  sexprC <- substitute(exprC)
  missingA <- missing(exprA)
  missingB <- missing(exprB)
  missingC <- missing(exprC)

  d <- decompose_expr(sexprA, sexprB, sexprC,
                      missingA, missingB, missingC,
                      .env = .env,
                      nThread = nThread)

  do_sum3s_par(
    d[[1]],
    d[[2]],
    d[[3]],
    d[[4]],
    d[[5]],
    d[[6]],
    d[[7]],
    d[[8]],
    d[[9]],
    d[[10]],
    d[[11]],
    d[[12]],
    d[[13]],
    d[[14]],
    d[[15]],
    d[[16]],
    d[[17]],
    d[[18]],
    d[[19]],
    d[[20]],
    d[[21]],
    d[[22]],
    d[[23]],
    d[[24]],
    ampersand = FALSE,
    nThread = nThread
  )
}

