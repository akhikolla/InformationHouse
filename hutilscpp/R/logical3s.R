#' Complex logical expressions
#' @name logical3s
#' @description Performant implementations of \code{&} et \code{or}.
#' Performance is high when the expressions are long (i.e. over 10M elements)
#' and in particular when they are of the form \code{lhs <op> rhs} for binary
#' \code{<op>}.
#' @param exprA,exprB,exprC,... Expressions of the form \code{x <op> y}.
#' with \code{<op>} one of the standard binary operators.
#'
#' Only \code{exprA} is required, all following expressions are optional.
#'
#' @param .parent_nframes \describe{
#' \item{\code{integer(1)}}{For internal use. Passed to \code{eval.parent}.}
#' }
#' @param nThread \describe{
#' \item{\code{integer(1)}}{Number of threads to use.}
#' }
#'
#'
#'
#'
#'
#' @return
#'
#' \code{and3s} and \code{or3s} return \code{exprA & exprB & exprC} and
#' \code{exprA | exprB | exprC} respectively. If any expression is missing
#' it is considered \code{TRUE} for \code{and3s} and \code{FALSE} for \code{or3s};
#' in other words only the results of the other expressions count towards the result.
#'
#'
#'
NULL


op2M <- function(operator) {
  switch(operator,
         "!=" = 1L,
         "==" = 2L,
         ">=" = 3L,
         "<=" = 4L,
         ">"  = 5L,
         "<"  = 6L,
         "%in%" = 7L,
         "%between%" = 8L,
         "%(between)%" = 9L,
         "%]between[%" = 10L)
}

is_binary_sexp <- function(sexprA, .parent_nframes = 2L) {
  # e.g.
  #      x > -1L
  isBinary <-
    is.call(sexprA) &&
    length(sexprA) == 3L &&
    (M <- op2M(as.character(sexprA[[1L]]))) &&
    is.name(lhs <- sexprA[[2L]])


  if (isBinary) {
    # c(1, 2)  gives rhs a language object, but rhs_eval is length-2 numeric
    # lhs <- sexprA[[2]]
    rhs <- sexprA[[3]]

    # lhs_eval <- eval.parent(lhs)


    attr(isBinary, "M") <- M
    attr(isBinary, "lhs") <- lhs
    attr(isBinary, "rhs") <- rhs
    attr(isBinary, "rhs_eval") <- rhs_eval <- eval.parent(rhs, n = .parent_nframes)

    if (OR(is.integer(rhs),
           AND(is.integer(rhs_eval),
               OR(length(rhs_eval) == 1L,
                  OR(M == op2M("%between%") && length(rhs_eval) == 2L,
                     M == op2M("%in%")))))) {
      attr(isBinary, "rhs_eval") <- rhs_eval
      return(isBinary)
    }
  }
  FALSE
}

is_lgl_sym <- function(sexpr, expr) {
  # Symbol first since is.logical will evaluate
  # Keep expr from evaluating so that the call can
  # be separated
  is.symbol(sexpr) && is.logical(expr)
}

is_lgl_negation <- function(sexpr, expr) {
  is.call(sexpr) && as.character(sexpr)[[1]] == "!" &&
    length(sexpr) > 1L &&
    is.symbol(sexpr[[2]])
}

#' @rdname logical3s
#' @export
and3s <- function(exprA, exprB, exprC, ..., .parent_nframes = 1L,
                  nThread = getOption("hutilscpp.nThread", 1L)) {
  nThread <- check_omp(nThread)

  sexprA <- substitute(exprA)
  sexprB <- substitute(exprB)
  sexprC <- substitute(exprC)

  X3 <- Y3 <- Z3 <- integer(0)
  A <- B <- C <- logical(0)
  x <- y <- z <- integer(0)
  ox <- oy <- oz <- -1L
  x1 <- y1 <- z1 <- 0L
  x2 <- y2 <- z2 <- 0L

  is_seq <- function(A) {
    identical(A, seq.int(A[1L], along.with = A))
  }

  isBinaryA <- is_binary_sexp(sexprA, .parent_nframes = .parent_nframes + 1L)
  if (isBinaryA) {
    x <- eval.parent(sexprA[[2]], n = .parent_nframes)
    ox <- attr(isBinaryA, "M")
    rhs_eval <- attr(isBinaryA, "rhs_eval")
    if (ox == op2M("%between%")) { # between so two elements
      x1 <- rhs_eval[[1]]
      x2 <- rhs_eval[[2]]
    } else if (ox == op2M("%in%")) {
      X3 <- rhs_eval

      x1 <- X3[1]
      if (is_seq(X3)) {
        x2 <- X3[length(X3)]
        ox <- op2M("%between%")
      } else {
        if (length(X3) > 100L) {
          A <- x %in% X3
        } else {
          A <- do_par_in(x, X3, nThread = nThread)
        }
        x <- integer(0)
        ox <- -1L
        x2 <- X3[2]
      }
    } else {
      x1 <- x2 <- attr(isBinaryA, "rhs_eval")
    }

  } else if (is_lgl_negation(sexprA, exprA)) {
    A <- eval.parent(sexprA[[2]], n = .parent_nframes)
    ox <- 1L
  } else {
    A <- exprA
  }
  if (!missing(exprB)) {
    isBinaryB <- is_binary_sexp(sexprB, .parent_nframes = .parent_nframes + 1L)
    if (isBinaryB) {
      y <- eval.parent(sexprB[[2]], n = .parent_nframes)
      oy <- attr(isBinaryB, "M")
      rhs_eval <- attr(isBinaryB, "rhs_eval")
      if (oy == op2M("%between%")) { # between so two elements
        y1 <- rhs_eval[[1]]
        y2 <- rhs_eval[[2]]
      } else if (oy == op2M("%in%")) {
        Y3 <- rhs_eval
        stopifnot(length(Y3) > 1)
        y1 <- Y3[1L]
        if (is_seq(Y3)) {
          # Treat x %in% a:b as x %between% c(a, b)
          y2 <- Y3[length(Y3)]
          oy <- op2M("%between%")
        } else {
          if (length(Y3) > 100L) {
            B <- y %in% Y3
          } else {
            B <- do_par_in(y, Y3, nThread = nThread)
          }
          y <- integer(0)
          oy <- -1L
          y2 <- Y3[2]
        }
      } else {
        y1 <- y2 <- rhs_eval
      }

    } else if (is_lgl_negation(sexprB, exprB)) {
      B <- eval.parent(sexprB[[2]], n = .parent_nframes)
      oy <- 1L
    } else {
      B <- exprB
    }
  }


  if (!missing(exprC)) {
    isBinaryC <- is_binary_sexp(sexprC, .parent_nframes = .parent_nframes + 1L)

    if (isBinaryC) {
      z <- eval.parent(sexprC[[2]], n = .parent_nframes)
      oz <- attr(isBinaryC, "M")
      rhs_eval <- attr(isBinaryC, "rhs_eval")
      if (oz == op2M("%between%")) { # between so two elements
        z1 <- rhs_eval[[1]]
        z2 <- rhs_eval[[2]]
      } else if (oz == op2M("%in%")) {
        Z3 <- rhs_eval
        stopifnot(length(Z3) > 1)
        z1 <- Z3[1]
        if (is_seq(Z3)) {
          z2 <- Z3[length(Z3)]
          oz <- op2M("%between%")
        } else {
          if (length(Z3) > 100L) {
            C <- z %in% Z3
          } else {
            C <- do_par_in(z, Z3, nThread = nThread)
          }
          z <- integer(0L)
          oz <- -1L
        }
      } else {
        z1 <- z2 <- rhs_eval
      }

    } else if (is_lgl_negation(sexprC, exprC)) {
      C <- eval.parent(sexprC[[2]], n = .parent_nframes)
      oz <- 1L
    } else {
      C <- exprC
    }
  }

  ans <-
    do_and3_par(x, ox, x1, x2,
                y, oy, y1, y2,
                z, oz, z1, z2,
                A, B, C,
                deparse(sexprA),
                nThread = nThread)

  if (missing(..1)) {
    return(ans)
  }

  and3s(ans, ..., .parent_nframes = .parent_nframes + 1L, nThread = nThread)
}

#' @rdname logical3s
#' @export

or3s <- function(exprA, exprB, exprC, ..., .parent_nframes = 1L,
                 nThread = getOption("hutilscpp.nThread", 1L)) {
  nThread <- check_omp(nThread)

  sexprA <- substitute(exprA)
  sexprB <- substitute(exprB)
  sexprC <- substitute(exprC)

  X3 <- Y3 <- Z3 <- integer(0)
  A <- B <- C <- logical(0)
  x <- y <- z <- integer(0)
  ox <- oy <- oz <- -1L
  x1 <- y1 <- z1 <- 0L
  x2 <- y2 <- z2 <- 0L

  is_seq <- function(A) {
    identical(A, seq.int(A[1L], along.with = A))
  }

  isBinaryA <- is_binary_sexp(sexprA)
  if (isBinaryA) {
    x <- eval.parent(sexprA[[2]], n = .parent_nframes)
    ox <- attr(isBinaryA, "M")
    rhs_eval <- attr(isBinaryA, "rhs_eval")
    if (ox == op2M("%between%")) { # between so two elements
      x1 <- rhs_eval[[1]]
      x2 <- rhs_eval[[2]]
    } else if (ox == op2M("%in%")) {
      X3 <- rhs_eval

      x1 <- X3[1]
      if (is_seq(X3)) {
        x2 <- X3[length(X3)]
        ox <- op2M("%between%")
      } else {
        if (length(X3) > 100L) {
          A <- x %in% X3
        } else {
          A <- do_par_in(x, X3, nThread = nThread)
        }
        x <- integer(0)
        ox <- -1L
        x2 <- X3[2]
      }
    } else {
      x1 <- x2 <- attr(isBinaryA, "rhs_eval")
    }

  } else if (is_lgl_negation(sexprA, exprA)) {
    A <- eval.parent(sexprA[[2]], n = .parent_nframes)
    ox <- 1L
  } else {
    A <- exprA
  }
  if (!missing(exprB)) {
    isBinaryB <- is_binary_sexp(sexprB)
    if (isBinaryB) {
      y <- eval.parent(sexprB[[2]], n = .parent_nframes)
      oy <- attr(isBinaryB, "M")
      rhs_eval <- attr(isBinaryB, "rhs_eval")
      if (oy == op2M("%between%")) { # between so two elements
        y1 <- rhs_eval[[1]]
        y2 <- rhs_eval[[2]]
      } else if (oy == op2M("%in%")) {
        Y3 <- rhs_eval
        stopifnot(length(Y3) > 1)
        y1 <- Y3[1L]
        if (is_seq(Y3)) {
          # Treat x %in% a:b as x %between% c(a, b)
          y2 <- Y3[length(Y3)]
          oy <- op2M("%between%")
        } else {
          if (length(Y3) > 100L) {
            B <- y %in% Y3
          } else {
            B <- do_par_in(y, Y3, nThread = nThread)
          }
          y <- integer(0)
          oy <- -1L
          y2 <- Y3[2]
        }
      } else {
        y1 <- y2 <- rhs_eval
      }

    } else if (is_lgl_negation(sexprB, exprB)) {
      # Make coverage explicit
      if (.parent_nframes > 1L) {
        B <- eval.parent(sexprB[[2]], n = .parent_nframes)
      } else {
        B <- eval.parent(sexprB[[2]], n = .parent_nframes)
      }
      oy <- 1L
    } else {
      B <- exprB
    }
  }


  if (!missing(exprC)) {
    isBinaryC <- is_binary_sexp(sexprC)

    if (isBinaryC) {
      z <- eval.parent(sexprC[[2]], n = .parent_nframes)
      oz <- attr(isBinaryC, "M")
      rhs_eval <- attr(isBinaryC, "rhs_eval")
      if (oz == op2M("%between%")) { # between so two elements
        z1 <- rhs_eval[[1]]
        z2 <- rhs_eval[[2]]
      } else if (oz == op2M("%in%")) {
        Z3 <- rhs_eval
        stopifnot(length(Z3) > 1)
        z1 <- Z3[1]
        if (is_seq(Z3)) {
          z2 <- Z3[length(Z3)]
          oz <- op2M("%between%")
        } else {
          if (length(Z3) > 100L) {
            C <- z %in% Z3
          } else {
            C <- do_par_in(z, Z3, nThread = nThread)
          }
          z <- integer(0L)
          oz <- -1L
        }
      } else {
        z1 <- z2 <- rhs_eval
      }

    } else if (is_lgl_negation(sexprC, exprC)) {
      C <- eval.parent(sexprC[[2]], n = .parent_nframes)
      oz <- 1L
    } else {
      C <- exprC
    }
  }
  ans <-
    do_or3_par(x, ox, x1, x2,
               y, oy, y1, y2,
               z, oz, z1, z2,
               A, B, C,
               deparse(sexprA),
               nThread = nThread)
  if (missing(..1)) {
    return(ans)
  }
  if (missing(..2)) {
    return(ans | ..1)
  }
  or3s(ans, ..., .parent_nframes = .parent_nframes + 1L, nThread = nThread)
}


