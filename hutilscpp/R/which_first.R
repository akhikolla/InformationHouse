#' Where does a logical expression first return \code{TRUE}?
#' @description A faster and safer version of \code{which.max} applied
#' to simple-to-parse logical expressions.
#'
#' @param expr An expression, such as \code{x == 2}.
#' @return The same as \code{which.max(expr)} or \code{which(expr)[1]} but returns \code{0L}
#'  when \code{expr} has no \code{TRUE} values.
#'
#'
#' @details
#' If the \code{expr} is of the form \code{LHS <operator> RHS}
#' and \code{LHS} is a single symbol, \code{operator} is one of
#' \code{==},  \code{!=}, \code{>}, \code{>=}, \code{<}, \code{<=},
#' \code{\%in\%},
#' or
#'  \code{\%between\%},
#' and \code{RHS} is numeric, then \code{expr} is not
#' evaluated directly; instead, each element of \code{LHS} is compared
#' individually.
#'
#' If \code{expr} is not of the above form, then \code{expr} is evaluated
#' and passed to \code{which.max}.
#'
#' Using this function can be significantly faster than the alternatives
#' when the computation
#' of \code{expr} would be expensive, though the difference is only likely to
#' be clear when \code{length(x)} is much larger than 10 million.
#' But even for smaller vectors, it has the benefit of returning
#' \code{0L} if none of the values in \code{expr} are \code{TRUE}, unlike
#' \code{which.max}.
#'
#' Compared to \code{\link[base:funprog]{Position}} for an appropriate
#'  choice of \code{f} the speed of \code{which_first} is not much faster
#'  when the expression is \code{TRUE} for some position. However, \code{which_first}
#'  is faster when all elements of \code{expr} are \code{FALSE}.
#'  Thus \code{which_first} has a smaller worst-case time than the
#'  alternatives for most \code{x}.
#'
#' Missing values on the RHS are handled specially.
#' \code{which_first(x \%between\% c(NA, 1))} for example is equivalent to
#' \code{which_first(x <= 1)}, as in \code{\link[data.table:between]{data.table::between}}.
#'
#'
#'
#'
#'
#' @param verbose \describe{
#' \item{\code{logical(1)}, default: \code{FALSE}}{If \code{TRUE} a message is emitted
#' if \code{expr} could not be handled in the advertised way.}
#' }
#'
#' @param reverse \describe{
#'   \item{\code{logical(1)}, default: \code{FALSE}}{Scan \code{expr} in reverse.}
#' }
#'
#' @param sexpr Equivalent to \code{substitute(expr)}. For internal use.
#' @param eval_parent_n Passed to \code{eval.parent}, the environment in which
#' \code{expr} is evaluated.
#'
#' @param suppressWarning Either a \code{FALSE} or \code{TRUE}, whether or not
#' warnings should be suppressed. Also supports a string input which suppresses a
#' warning if it matches as a regular expression.
#'
#' @param use.which.max If \code{TRUE}, \code{which.max} is dispatched immediately,
#' even if \code{expr} would be amenable to separation. Useful when evaluating
#' many small \code{expr}'s when these are known in advance.
#'
#'
#'
#' @examples
#'
#' N <- 1e5
#' # N <- 1e8  ## too slow for CRAN
#'
#' # Two examples, from slowest to fastest,
#' # run with N = 1e8 elements
#'
#'                                        # seconds
#' x <- rep_len(runif(1e4, 0, 6), N)
#' bench_system_time(x > 5)
#' bench_system_time(which(x > 5))        # 0.8
#' bench_system_time(which.max(x > 5))    # 0.3
#' bench_system_time(which_first(x > 5))  # 0.000
#'
#' ## Worst case: have to check all N elements
#' x <- double(N)
#' bench_system_time(x > 0)
#' bench_system_time(which(x > 0))        # 1.0
#' bench_system_time(which.max(x > 0))    # 0.4  but returns 1, not 0
#' bench_system_time(which_first(x > 0))  # 0.1
#'
#' x <- as.character(x)
#' # bench_system_time(which(x == 5))     # 2.2
#' bench_system_time(which.max(x == 5))   # 1.6
#' bench_system_time(which_first(x == 5)) # 1.3
#'
#' @export which_first


which_first <- function(expr,
                        verbose = FALSE,
                        reverse = FALSE,
                        sexpr,
                        eval_parent_n = 1L,
                        suppressWarning = getOption("hutilscpp_suppressWarning", FALSE),
                        use.which.max = FALSE) {
  if (use.which.max) {
    return(.which_first(expr, verbose = verbose, reverse = reverse))
  }
  rhs <- NULL
  if (missing(sexpr)) {
    sexpr <- substitute(expr)
  }
  if (!is.call(sexpr) ||
      length(sexpr) != 3L ||
      !(op <- do_op2M(operator <- as.character(sexpr[[1L]])))) {
    o <- .which_first(expr, verbose = verbose, reverse = reverse)
    return(o)
  }
  lhs <- sexpr[[2L]]

  if (!is.name(lhs)) {
    o <- .which_first(expr, verbose = verbose, reverse = reverse)
    return(o)
  }

  rhs <- sexpr[[3L]]

  isValidExpr <-
    is.numeric(rhs) ||  # bare doubles and integers
    AND(is.call(rhs) && is.numeric(rhs[[2L]]),
        # negatives
        as.character(rhs[[1L]]) == "-") ||
    is.logical(rhs_eval <- eval.parent(rhs, n = eval_parent_n)) ||

    # c(0, 1, 2) is not numeric but it is when evaluated
    AND(op >= 7 && op <= 10,
        is.numeric(rhs_eval)) ||

    # x != y  or  x == y

    AND(is.symbol(rhs),
        # Must be both integers or both doubles otherwise will fall through too late.
        # Too expensive to check or coerce to doubles in case.
        AND(is.numeric(lhs_eval <- eval.parent(lhs, n = eval_parent_n)),
            is.numeric(rhs_eval))) ||

    AND(is.character(rhs_eval),
        # Must be both integers or both doubles otherwise will fall through too late.
        # Too expensive to check or coerce to doubles in case.
        all(is.character(lhs_eval <- eval.parent(lhs, n = eval_parent_n)),
            OR(length(lhs_eval) == length(rhs_eval),
               length(rhs_eval) == 1L),
            OR(operator == "!=",
               operator == "=="),
            na.rm = TRUE))



  if (!isValidExpr) {
    o <- .which_first(expr, verbose = verbose, reverse = reverse)
    return(o)
  }

  lhs_eval <- eval.parent(lhs, n = eval_parent_n)

  if (length(lhs_eval) <= 1L) {
    if (isTRUE(expr)) {
      return(1L)
    } else {
      return(0L)
    }
  }
  rhs_eval <- eval.parent(rhs, n = eval_parent_n)


  if (is.logical(lhs_eval)) {
    if (operator == "%in%") {
      if (reverse) {
        o <- do_which_last_in_lgl(lhs_eval,
                                  anyNA_ = anyNA(rhs_eval),
                                  any_ = any(rhs_eval, na.rm = TRUE),
                                  nall_ = !all(rhs_eval, na.rm = TRUE))
      } else {
        o <- do_which_first_in_lgl(lhs_eval,
                                   anyNA_ = anyNA(rhs_eval),
                                   any_ = any(rhs_eval, na.rm = TRUE),
                                   nall_ = !all(rhs_eval, na.rm = TRUE))
      }
      return(R_xlen_t(o))
    }

    if (length(rhs_eval) == length(lhs_eval) ||
        operator == "%in%" ||
        operator == "%between%" ||
        operator == "%(between)%") {
      if (is.logical(rhs_eval)) {
        # is.logical(rhs_eval)  necessary to ensure integers don't falsely resemble TRUE
        # e.g. which_first(c(TRUE, TRUE) != c(1L, 2L)) should be 2 not 0.
        # Quicker to do this, apparently.
        if (identical(as.logical(lhs_eval),
                      as.logical(rhs_eval))) {
          if (operator == "==" || operator == "<=" || operator == ">=" || operator == "%in%") {
            return(1L)
          } else {
            return(0L)
          }
        } else {
          return(do_which_first_lgl_lgl_op(lhs_eval, rhs_eval, op, reverse = reverse))
        }
      }
    }

    if (anyNA(rhs_eval)) {
      warn_msg <- "`rhs` appears to be logical NA. Treating as\n\twhich_first(is.na(.))"
      if (operator == "!=") {
        warn_msg <- "`rhs` appears to be logical NA. Treating as\n\twhich_first(!is.na(.))"
      }
      suppressThisWarn <-
        isTRUE(suppressWarning) ||
        AND(is.character(suppressWarning),
            any(vapply(suppressWarning, grepl, x = warn_msg, perl = TRUE, FUN.VALUE = FALSE),
                na.rm = TRUE))

      switch(operator,
             "==" = {
               if (!suppressThisWarn) {
                 warning(warn_msg)
               }
               o <- .which_first(is.na(lhs_eval), verbose = verbose, reverse = reverse)
             },
             "!=" = {
               if (!suppressThisWarn) {
                 warning(warn_msg)
               }
               o <- .which_first(!is.na(lhs_eval), verbose = verbose, reverse = reverse)
             },
             stop("`rhs` appears to be logical NA. This is not supported for operator '", operator, "'."))
    } else {
      o <- .which_first_logical(lhs_eval, as.logical(rhs_eval), operator = operator, rev = reverse)
    }
    return(R_xlen_t(o))
  }

  if (is.numeric(lhs_eval) && is.numeric(rhs_eval)) {
    if (op == 7L) {
      o <-
        if (reverse) {
          switch(typeof(lhs_eval),
                 "integer" = switch(typeof(rhs_eval),
                                    "integer" = do_which_last_xi_ini(lhs_eval, rhs_eval),
                                    "double" =  do_which_last_xi_ind(lhs_eval, rhs_eval)),
                 "double" = do_which_last_xd_ind(lhs_eval, rhs_eval, anyNA(rhs_eval)))
        } else {
          switch(typeof(lhs_eval),
                 "integer" = switch(typeof(rhs_eval),
                                    "integer" = do_which_first_xi_ini(lhs_eval, rhs_eval),
                                    "double" =  do_which_first_xi_ind(lhs_eval, rhs_eval)),
                 "double" = do_which_first_xd_ind(lhs_eval, rhs_eval, anyNA(rhs_eval)))
        }
      return(R_xlen_t(o))
    }

    nx <- length(lhs_eval)
    ny <- length(rhs_eval)

    # op2M("%in%") = 7
    # op2M("%between%") > 7
    if (op > 7L && op <= 10L && ny != 2) {
      stop("Expression in `which_first` was of the form:\n\t which_first(x %between% rhs)\n",
           "yet `length(rhs) = ", ny, "`. Ensure the object passed to RHS is an atomic vector ",
           "of length two.")
    }
    if (op < 7L) {
      if (nx == 0L || ny == 0L) {
        return(0L)
      }
      # simple binary op
      if (nx != ny && ny != 1L) {
        stop("In `expr ~ <lhs> <op> <rhs>`, length of lhs = ", nx, " but length of rhs = ", ny, ". ",
             "With operator '", operator, "' the only permitted lengths of the RHS are 1 and length(lhs).")
      }
      if (ny == 2L) {
        # Otherwise will be redirected to between switch
        ny <- 3L
      }
    }

    if (is.integer(rhs_eval)) {
      y1d <- y1i <- rhs_eval[1L]
      y2d <- y2i <- rhs_eval[2L]

    } else {
      y1d <- y1i <- rhs_eval[1L]
      y2d <- y2i <- rhs_eval[2L]
    }

    if (reverse) {
      o <- do_which_last__(lhs_eval, op, rhs_eval,
                           ny = qd2i(ny),
                           y1i = qd2i(y1i),
                           y2i = qd2i(y2i),
                           y1d = y1d,
                           y2d = y2d)
    } else {
      o <- do_which_first__(lhs_eval, op, rhs_eval,
                            ny = qd2i(ny),
                            y1i = qd2i(y1i),
                            y2i = qd2i(y2i),
                            y1d = y1d,
                            y2d = y2d)
      # o <- .which_first(expr)
    }
    return(R_xlen_t(o))
  }

  if (is.character(lhs_eval)) {
    oc <-
      switch(operator,
             "==" = AnyCharMatch(lhs_eval, as.character(rhs_eval)),
             "!=" = AnyCharMatch(lhs_eval, as.character(rhs_eval), opposite = TRUE),
             {
               o <- .which_first(expr, verbose = verbose, reverse = reverse)
             })
    if (oc <= .Machine$integer.max) {
      oc <- as.integer(oc)
    }
    return(oc)
  }

  o <- .which_first(expr, reverse = reverse)

  R_xlen_t(o)
}

.which_first <- function(expr, verbose = FALSE, reverse = FALSE) {
  if (verbose) {
    message("Falling back to `which.max(expr)`.")
  }
  if (reverse) {
    return(do_which_last(expr))
  }
  o <- which.max(expr)
  if (length(o) == 0L) {
    # i.e. all expr NA
    return(0L)
  }
  # o == 1L is wrong if all expr are FALSE
  if (o == 1L && !expr[1L]) {
    o <- 0L
  }
  o
}

.which_first_logical <- function(lhs, rhs, operator = "==", verbose = FALSE, rev = FALSE) {
  stopifnot(length(lhs) >= 1L,
            is.logical(rhs), length(rhs) == 1L, !anyNA(rhs),
            is.logical(verbose))

  #    operator   rhs   first/last
  # 1:       != FALSE   notFALSE
  # 2:       !=  TRUE   notTRUE
  switch({
    operator
  },
  "!=" = {
    if (rhs) {
      if (rev) {
        return(do_which_last_false(lhs))
      } else {
        return(do_which_first_false(lhs))
      }
    } else {
      if (rev) {
        return(do_which_last(lhs))
      } else {
        return(do_which_first(lhs))
      }
    }
  },
  "==" = {
    # 0:       ==  TRUE   TRUE
    # 9:       == FALSE   FALSE
    if (rhs) {
      if (rev) {
        return(do_which_last(lhs))
      } else {
        return(do_which_first(lhs))
      }
    } else {
      if (rev) {
        return(do_which_last_false(lhs))
      } else {
        return(do_which_first_false(lhs))
      }
    }
  },
  "<" = {
    # 5:        < FALSE   0L
    # 6:        <  TRUE   FALSE
    if (rhs) {
      if (rev) {
        return(do_which_last_false(lhs))
      } else {
        return(do_which_first_false(lhs))
      }
    } else {
      return(0L)
    }
  },
  "<=" = {
    # 8:       <=  TRUE   1L/length
    # 7:       <= FALSE   FALSE
    if (rhs) {
      if (rev) {
        return(length(lhs))
      } else {
        return(1L)
      }
    } else {
      if (rev) {
        return(do_which_last_false(lhs))
      } else {
        return(do_which_first_false(lhs))
      }
    }
  },
  ">" = {
    # 2:        >  TRUE   0L
    # 1:        > FALSE   TRUE
    if (rhs) {
      return(0L)
    } else {
      if (rev) {
        return(do_which_last(lhs))
      } else {
        return(do_which_first(lhs))
      }
    }
  },
  ">=" = {
    # 4:       >=  TRUE   TRUE
    # 3:       >= FALSE   TRUE/FALSE
    if (rhs) {
      if (rev) {
        return(do_which_last(lhs))
      } else {
        return(do_which_first(lhs))
      }
    } else {
      if (rev) {
        return(length(lhs))
      } else {
        return(1L)
      }
    }
  })
  # nocov start
  print(ls.str())
  stop("Internal error: 559:2019:11:18. ")
  # nocov end
}


#' @rdname which_first
#' @export which_last
which_last <- function(expr,
                       verbose = FALSE,
                       reverse = FALSE,
                       suppressWarning = getOption("hutilscpp_suppressWarning", FALSE)) {
  force(reverse)
  which_first(expr,
              verbose = verbose,
              reverse = isFALSE(reverse),
              sexpr = substitute(expr),
              eval_parent_n = 2L,
              suppressWarning = suppressWarning)
}








