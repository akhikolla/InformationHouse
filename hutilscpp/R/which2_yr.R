which2_Year <- function(expr1, expr2, Year = NULL, yr = NULL) {

  opx <- as.character(substitute(expr1)[[1]])
  eqx <- opx == "==" || opx == "<=" || opx == ">="
  gtx <- opx == ">=" || opx == ">"
  ltx <- opx == "<=" || opx == "<"

  opy <- as.character(substitute(expr2)[[1]])
  eqy <- opy == "==" || opy == "<=" || opy == ">="
  gty <- opy == ">=" || opy == ">"
  lty <- opy == "<=" || opy == "<"

  lhs1 <- eval.parent(substitute(expr1)[[2]])
  rhs1 <- eval.parent(substitute(expr1)[[3]])

  lhs2 <- eval.parent(substitute(expr2)[[2]])
  rhs2 <- eval.parent(substitute(expr2)[[3]])

  stopifnot(length(rhs1) == 1L,
            is.integer(rhs1),
            !anyNA(rhs1),
            length(rhs2) == 1L,
            is.integer(rhs2),
            !anyNA(rhs2))

  consider_yr <- !is.null(Year)
  if (is.null(Year)) {
    if (!is.null(yr)) {
      warning("`Year = NULL`, but `yr` has been provided, ",
              "and will be ignored.")
    }
    Year <- integer(0L)
    yr <- 0L
  } else {
    if (is.null(yr)) {
      stop("`Year` was provided but `yr` remained NULL.")
    }
    if (!is.integer(Year)) {
      stop("`Year` was type ", typeof(Year), ", but must be (strictly) type integer.")
    }
    if (length(Year) != length(lhs2)) {
      stop("`Year` had length ", length(Year), ", yet `expr1` evaluates to ",
           " a length-", length(lhs1), " vector.")
    }
  }

  do_which2_yr(Year, yr, consider_yr = consider_yr,
               lhs1, rhs1, eqx, gtx, ltx,
               lhs2, rhs2, eqy, gty, lty)

}
