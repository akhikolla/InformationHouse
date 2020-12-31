

# returns the safe integer for a single double
# iop2safe_int <- function(op, dbl) {
#   is_wholer <- function(dbl) dbl >= -2147483647 && dbl <= 2147483647 && dbl == as.integer(dbl)
#
#   switch(op,
#          # "!=" true iff dbl is whole number, so != a if a whol numebr and falsey
#          {
#            if (dbl >= -2147483647 && dbl <= 2147483647 && dbl == as.integer(dbl)) {
#              return(list(op, as.integer(dbl)))
#            } else {
#              return(list(-1L, 0L))
#            }
#          },
#          # "==" true
#          list())
# }


decompose_expr <- function(sexprA,
                           sexprB,
                           sexprC,
                           missingA,
                           missingB,
                           missingC,
                           .parent_nframes = 1L,
                           .env = parent.frame(.parent_nframes),
                           nThread = 1L) {
  X3 <- Y3 <- Z3 <- integer(0)
  A <- B <- C <- logical(0)
  x <- y <- z <- integer(0)
  ox <- oy <- oz <- -1L
  x1 <- y1 <- z1 <- 0L
  x2 <- y2 <- z2 <- 0L

  is_binary_sexp <- function(sexprA, .env) {
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
      attr(isBinary, "rhs_eval") <- rhs_eval <- eval(rhs, envir = .env)

      if (OR(is.numeric(rhs),
             AND(is.numeric(rhs_eval),
                 OR(length(rhs_eval) == 1L,
                    OR(M == op2M("%between%") && length(rhs_eval) == 2L,
                       M == op2M("%in%")))))) {
        attr(isBinary, "rhs_eval") <- rhs_eval
        return(isBinary)
      }
    }
    FALSE
  }

  is_seq <- function(A) {
    identical(A, seq.int(A[1L], along.with = A))
  }

  isBinaryA <- is_binary_sexp(sexprA, .env)
  if (isBinaryA) {
    x <- eval(sexprA[[2]], envir = .env)
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

  } else if (is_lgl_negation(sexprA)) {
    A <- eval(sexprA[[2]], envir = .env)
    ox <- 1L
  } else {
    A <- eval(sexprA, envir = .env)
  }
  if (!missingB) {
    isBinaryB <- is_binary_sexp(sexprB, .env)
    if (isBinaryB) {
      y <- eval(sexprB[[2]], envir = .env)
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

    } else if (is_lgl_negation(sexprB)) {
      B <- eval(sexprB[[2]], envir = .env)
      oy <- 1L
    } else {
      B <- eval(sexprB, envir = .env)
    }
  }


  if (!missingC) {
    isBinaryC <- is_binary_sexp(sexprC, .env)

    if (isBinaryC) {
      z <- eval(sexprC[[2]], envir = .env)
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

    } else if (is_lgl_negation(sexprC)) {
      C <- eval(sexprC[[2]], envir = .env)
      oz <- 1L
    } else {
      C <- eval(sexprC, envir = .env)
    }
  }

  # Handle double vectors in expression
  # Issues:
  # make sure the correct type is not coerced
  # make sure if lhs double, the integer is 0-length
  #
  # if it's lhs int then the rhs needs care
  # e.g. 1:5 >= 2.6 !=  1:5 >= 2L


  if (is.double(x)) {
    xd <- x
    x <- integer(0)
    # these are likely doubles, and could exceed range, causing spurious warnings
    # when passed (coercively) to Rcpp. Set to an integer and assume downstream
    # ignore
    x1d <- x1
    x2d <- x2
    x1 <- x2 <- NA_integer_
  } else {
    xd <- double(0)
    x1d <- x2d <- NA_real_
    if (is_wholer(x1) && is_wholer(x2)) {
      x1 <- as.integer(x1)
      x2 <- as.integer(x2)
    } else {
      A <- eval(sexprA, envir = .env)
      x <- integer(0)
    }
  }
  if (is.double(y)) {
    yd <- y
    y <- integer(0)
    y1d <- y1
    y2d <- y2
    y1 <- y2 <- NA_integer_
  } else {
    yd <- double(0)
    y1d <- y2d <- NA_real_
    if (is_wholer(y1) && is_wholer(y2)) {
      y1 <- as.integer(y1)
      y2 <- as.integer(y2)
    } else {
      B <- eval(sexprB, envir = .env)
      y <- integer(0)
    }
  }

  if (is.double(z)) {
    zd <- z
    z <- integer(0)
    z1d <- z1
    z2d <- z2
    z1 <- z2 <- NA_integer_
  } else {
    zd <- double(0)
    z1d <- z2d <- NA_real_
    if (is_wholer(z1) && is_wholer(z2)) {
      z1 <- as.integer(z1)
      z2 <- as.integer(z2)
    } else {
      C <- eval(sexprC, envir = .env)
      z <- integer(0)
    }
  }



  list(x, ox, x1, x2,
       xd, x1d, x2d,
       y, oy, y1, y2,
       yd, y1d, y2d,
       z, oz, z1, z2,
       zd, z1d, z2d,
       A, B, C)
}










