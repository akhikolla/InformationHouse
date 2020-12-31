
# Function is.submodel
# checks "submodel"-condition
#   x      is a character vector containing asf and csf components
#   y      is a character expresion (vector of length 1) containing a asf or csf
#   strict Logical: apply "strict" submodel condition?
# Returns a logical vectror of the same length as x
is.submodel <- function(x, y, strict = FALSE){
  if (length(y) != 1){
    stop("y must be a character vector of length 1 expressing an asf or csf")
  }
  x <- noblanks(x)
  y <- noblanks(y)

  # Check input conditions [identical code in cyclic()]
  xy <- unique(c(x, y))
  px <- lapply(xy, tryparse)
  ok <- !vapply(px, is.null, logical(1)) 
  ct_type <- if (any(grepl("=", xy, fixed = T))) "mv" else "cs"
  if (ct_type == "cs"){
    vals <- unique.default(unlist(lapply(px, all.vars)))
  } else {
    vals <- rapply(px, .call2list, how = "unlist",
                   stopOps = c("==", "<", ">", "<=", ">="), 
                   validOps = c("<-", "<<-", "=", "&", "|", "(", "-"))
    vals <- unique.default(vals[!vapply(vals, is.symbol, FUN.VALUE = TRUE)])
    vals <- sub(" == ", "=", vapply(vals, deparse, character(1)))
  }
  cond_type <- .qcondType(xy, values = vals, ct_type = ct_type, stdComplex.multiple.only = FALSE)
  ok <- ok & cond_type %in% c("stdAtomic", "stdComplex")
  if (any(!ok)) 
    stop("Invalid input to is.submodel:\n", paste0("  ", xy[!ok], collapse = "\n"),
         call. = FALSE)
  
  yy <- extract_asf(y)[[1]]
  xx <- extract_asf(x)

  # check and prepare y
  lhsy <- lhs(yy)
  rhsy <- rhs(yy)
  if (!all(nzchar(c(lhsy, rhsy))) || length(lhsy) != length(rhsy) || length(lhsy) <= 0)
    stop("Check the structure of y: should be a asf or csf.")
  lhsy <- setNames(lhs(yy), rhs(yy))

  # check and prepare x
  if (length(x) == 0) return(logical(0))
  ux <- unlist(xx)
  llx <- lengths(xx)
  lux <- lhs(ux)
  rux <- rhs(ux)
  if (!all(nzchar(c(lux, rux))) || length(lux) != length(rux) || length(lux) <= 0)
    stop("Check the structure of x: should contain asf or csf.")

  clx <- hstrsplit(lux, c("+", "*"), relist = FALSE)
  cly <- hstrsplit(lhsy, c("+", "*"), relist = FALSE)
  nms <- unique.default(c(clx, cly))
  
  ilx <- match(clx, table = nms)
  lx2 <- attr(clx, "lengths")[[2]]
  ilx <- ilx[order(rep(seq_along(lx2), lx2), ilx) ]
  rilx <- hrelist(ilx, attr(clx, "lengths"), relist = FALSE)

  ily <- match(cly, table = nms)
  ly2 <- attr(cly, "lengths")[[2]]
  ily <- ily[order(rep(seq_along(ly2), ly2), ily) ]
  rily <- hrelist(ily, attr(cly, "lengths"), relist = FALSE)
  names(rily) <- names(lhsy)
  
  ok <- rep(NA, sum(lengths(xx)))
  # test for non-strict submodels
  for (r in unique(rux)){
    ok[rux == r] <- C_is_submodel(rilx[rux == r], rily[[r]], strict = FALSE)
  }
  out <- m_all(C_relist_Log(ok, llx))
  # test for strict submodels if necessary
  # THIS COULD PROBABLY BE DONE MORE EFFICIENTLY!
  if (strict && any(out)){
    ok <- rep(NA, sum(lengths(xx)))
    for (r in unique(rux)){
      #ok[rux == r] <- C_is_submodel(rilx[rux == r], rily[[r]], strict = TRUE)
      ok[rux == r] <- vapply(rilx[rux == r], intList_equal, rily[[r]], FUN.VALUE = logical(1))
    }
    out <- out & !m_all(C_relist_Log(ok, llx))
  }
  names(out) <- x
  attr(out, "target") <- y
  out
}

# Function identical.model
# checks identity of two csf
identical.model <- function(x, y){
  if (length(x) != 1){
    stop("x must be a character vector of length 1 expressing an asf or csf")
  }
  if (length(y) != 1){
    stop("y must be a character vector of length 1 expressing an asf or csf")
  }
  x <- noblanks(x)
  y <- noblanks(y)
  is.submodel(x, y, strict = FALSE) && is.submodel(y, x, strict = FALSE)
}
