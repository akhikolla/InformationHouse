
# check if conditions are minimal wrt a reference set of conds
#   cond, ref     character strings
#   nms           optional: names of all factor values (either A, a or A=1, A=2, ...)
#   ignore.equals ignore conds in ref that are identical to the specific cond.
#                 If ref=cond and ignore.equals=FALSE, the result will always be a vector of FALSEs
#   sorted        Logical: Are x and y sorted in ascending order?
# Value:  logical vector with same length as cond
is.minimal <- function(cond, ref = cond, nms = NULL, ignore.equals = missing(ref), 
                       use.names = TRUE, sorted = FALSE){
  cond <- noblanks(cond)
  missref <- missing(ref)
  if (!missref) ref <- noblanks(ref)
  force(ignore.equals)
  
  rcl <- hstrsplit(cond, c("+", "*"), relist = FALSE)
  rcl_ref <- if (!missref){
    hstrsplit(ref, c("+", "*"), relist = FALSE)
  } else {
    rcl
  }
  if (is.null(nms)){
    nms <- if (missref){
      unique.default(rcl)
    } else {
      unique.default(c(rcl, rcl_ref))
    }
  }
  ril <- match(rcl, table = nms)
  if (!sorted){
    ll <- attr(rcl, "lengths")[[2]]
    ril <- ril[order(rep.int(seq_along(ll), ll), ril)]
  }
  if (missref){
    ril_ref <- ril
  } else {
    ril_ref <- match(rcl_ref, table = nms)
    if (!sorted){
      ll_ref <- attr(rcl_ref, "lengths")[[2]]
      ril_ref <- ril_ref[order(rep.int(seq_along(ll_ref), ll_ref), ril_ref)]
    }
  }
  out <- C_minimal(hrelist(ril, attr(rcl, "lengths")),
                   hrelist(ril_ref, attr(rcl_ref, "lengths")),
                   strict = ignore.equals)
  if (use.names){
    setNames(out, cond)
  } else {
    out
  }
}
