
# function selectCases:
# =====================
# Subsets a data frame x based on a condition cond
# Arguments:
#   x     data frame
#   cond  character string specifying a condition on the variables in x
#   type
# Value: Reduced data frame
selectCases <- function(cond, x = full.ct(cond), type, cutoff = 0.5, 
                        rm.dup.factors = FALSE, rm.const.factors = FALSE){
  if (inherits(x, "configTable")){
    type <- attr(x, "type")
  } else {
    if (missing(type)) type <- "cs"    # "hidden" Default value!
    x <- configTable(x, type = type, rm.dup.factors = FALSE, 
                     rm.const.factors = FALSE, verbose = FALSE)
  }
  stopifnot(length(cond) == 1)
  co <- condition.default(cond, x, force.bool = TRUE)[[1]]
  if (inherits(co, "invalidCond")){
    stop("The condition is invalid (", reason(co), "): ", format.condString(cond))
  }
  x[co[[1]] >= cutoff, , rm.dup.factors = rm.dup.factors, rm.const.factors = rm.const.factors]
}


# selectCases1: Mit Moeglichkeit der Vorgabe von con- und cov-Werten
# macht keinen Unterschied zwischen "->" und "<->"!!
#   selectCases("A->B", ...) entspricht selectCases1 mit con=1 und cov=0
#   selectCases("A<->B", ...) entspricht selectCases1 mit con=1 und cov=1
selectCases1 <- function(cond, x = full.ct(cond), type, con = 1, cov = 1, 
                         rm.dup.factors = FALSE, rm.const.factors = FALSE){
  if (inherits(x, "configTable")){
    type <- attr(x, "type")
  } else {
    if (missing(type)) type <- "cs"    # "hidden" Default value!
      x <- configTable(x, type = type, rm.dup.factors = FALSE, 
                       rm.const.factors = FALSE, verbose = FALSE)
  }
  # Check inputs
  if (length(con) != 1 || con < 0 || con > 1 || 
      length(cov) != 1 || cov < 0 || cov > 1){
    stop("Invalid input for 'con' or 'cov'")
  }
  if (length(cond) != 1){
    stop("'cond' must have length 1.")
  }

  # Slightly reduce con and cov values to avoid failing to find conditions due to rounding issues
  d.eps <- nrow(x) * .Machine$double.eps
  con <- con - d.eps
  cov <- cov - d.eps
  
  a <- condition.default(cond, x, rm.parentheses = TRUE)[[1]]
  if (inherits(a, "invalidCond"))
    stop("The condition is invalid (", reason(a), "): ", format.condString(cond))
  if (!inherits(a, "atomicCond")) stop("selectCases1 only works with a condition of type 'atomic'.")
  a <- ct2df(a)
  d <- a[[1L]] - a[[2L]]        # differences
  r <- rank(d, ties.method = "random") # ranks (random within ties)
  # lower and upper limiting values for potential subsets
  loVal <- sort(r[d<0], decreasing = TRUE)
  upVal <- sort(r[d>0])
  if (any(d==0)){
    loVal <- c(loVal[1]+1L, loVal)
    upVal <- c(upVal[1]-1L, upVal)
  }
  if ((length(loVal) == 1L) && is.na(loVal)) loVal <- 1L
  if ((length(upVal) == 1L) && is.na(upVal)) upVal <- length(r)
  # Build up matrices of con and cov
  # First define the auxiliary function .cc
  .cc <- function(subset = TRUE){
    x <- a[[1]][subset]
    y <- a[[2]][subset]
    Sminxy <- sum(pmin(x, y))
    Sminxy / c(con = sum(x), cov = sum(y))
  }
  ccmat <- array(NA_real_, dim = c(length(loVal), length(upVal), 2L))
  for (i in seq_along(loVal)){
    for (j in seq_along(upVal)){
      subs <- d == 0 | (r >= loVal[i] & r <= upVal[j])
      ccmat[i, j, ] <- .cc(subset = subs)
    }
  }
  rm(i, j, subs)
  conMat <- array(ccmat[, , 1], dim = dim(ccmat)[1:2])
  covMat <- array(ccmat[, , 2], dim = dim(ccmat)[1:2])
  # Identify subset that meet the con- and cov-requirement
  okMat <- !is.na(conMat) & conMat >= con & !is.na(covMat) & covMat >= cov
  if (!any(okMat)) stop("Not possible!")
  rcMat <- row(okMat)+col(okMat)
  rcmax <- max(rcMat[okMat])
  wmax <- which(rcMat==rcmax & okMat)
  # return the corresponding subset
  imax <- row(okMat)[wmax]
  jmax <- col(okMat)[wmax]
  subs <- d == 0 | (r >= loVal[imax] & r <= upVal[jmax])
  configTable(ct2df(x)[subs, ], type = attr(x, "type"), 
              rm.dup.factors = rm.dup.factors, rm.const.factors = rm.const.factors)
}

