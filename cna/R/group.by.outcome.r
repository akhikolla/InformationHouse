
# function group.by.outcome
#   group a list of atomicCond's by their outcome
group.by.outcome <- function(condlst, cases = TRUE){
  stopifnot(is.list(condlst))
  if (!all(isAtomic <- vapply(condlst, inherits, logical(1), "atomicCond"))){
    condlst <- condlst[isAtomic]
    message("group.by.outcomes only considers conditions of type 'atomic' - other are ignored.")
  }
  outc <- rhs(names(condlst))
  outc[!nzchar(outc)] <- "(No outcome)"
  out <- lapply(split(condlst, outc), grbyout1, cases = cases)
  structure(out, class = c("groupedConds", "listof"), cases = cases)
}
grbyout1 <- function(x, cases){
  secondCols <- do.call(cbind, lapply(x, "[", 2))
  outName <- names(secondCols)[1]
  if (!all(secondCols == secondCols[, 1]))
    stop("Response ", outName, "is not identical in all condition tables.")
  ncases <- do.call(cbind, lapply(x, "attr", "n"))
  if (!all(ncases == ncases[, 1]))
    stop("n is not identical in all condition tables.")
  ncases <- ncases[, 1, drop = TRUE]
  if (cases){
    Cases <- unname(lapply(x, "attr", "cases"))
    if (!all(duplicated(Cases)[-1]))
      stop("Cases are not identical in all condition tables.")
    Cases <- C_mconcat(Cases[[1]], sep = ",")
  }
  out <- do.call(data.frame, c(lapply(x, function(a) as.data.frame(a)[1]), list(check.names = FALSE)))
  out[[outName]] <- x[[1]][[outName]]
  out$n.obs <- ncases
  if (cases) rownames(out) <- Cases
  class(out) <- c("groupedC", "data.frame")
  out
}

print.groupedConds <- function(x, cases, ...)
  print.listof(x, row.names = attr(x, "cases"), ...)

print.groupedC <- function(x, ...){
  printDfwithBars(x, ..., bar.after = ncol(x) - 1)
}
