
# Auxiliary function .redCsf
# takes a logical matrix x, as returned by redundant()
# Returns a character vector of conditions resulting from reducing the conds in
# rownames(x) according to x
.redCsf <- function(x){
  nms <- names(x)
  names(x) <- NULL
  which_asf_red <- lapply(x, which)
  orig <- mapply(rep, seq_along(which_asf_red), lengths(which_asf_red), SIMPLIFY = FALSE)
  orig <- unlist(orig)
  
  asfList <- extract_asf(nms)
  asfList <- mapply("[", asfList[orig], -unlist(which_asf_red), SIMPLIFY = FALSE)
  
  out <- C_mconcat(happly(asfList, function(x) paste0("(", x, ")")), sep = "*")
  structure(out, originalCond = orig)
}  

# recursive function .minCsf
# takes cond and (either cti or full.cti)
# returns cond without redundancies
.minCsf <- function(cond, cti = NULL, full.cti = NULL, verbose = FALSE){
  ids <- attr(cond, "id")
  if (is.null(attr(cond, "id")))
    ids <- as.list(seq_along(cond))

  if (is.null(full.cti)){
    full.cti <- full.ct(cti)
  }
  rr <- redundant(cond, full.cti, simplify = FALSE)
  anyRed <- m_any(rr)
  n.asf <- lengths(rr, use.names = FALSE)
  if (verbose)
    cat("csf-lengths ", n.asf, ": ",
        sum(anyRed), " of ", length(cond), " csf are reducible ",
        sep = "")
  noRed <- names(rr)[!anyRed]
  attr(noRed, "id") <- ids[!anyRed]
  attr(noRed, "n.asf") <- n.asf[!anyRed]
  if (!any(anyRed)){
    if (verbose) cat("\n")
    return(noRed)
  }
  withRed <- rr[anyRed]
  condNew <- .redCsf(withRed)
  idsNew <- ids[anyRed][attr(condNew, "originalCond")]
  splitId <- split(idsNew, condNew)
  dups <- duplicated(condNew)
  if (verbose) cat("(", sum(!dups), " uniques)\n", sep = "")
  condNew <- condNew[!dups]
  condNew <- structure(condNew, 
                       id = unname(splitId[condNew]),
                       n.asf = n.asf[anyRed][!dups] - 1L)
  recurs <- .minCsf(condNew, NULL, full.cti = full.cti, verbose = verbose)
  structure(c(recurs, noRed), 
            id = lapply(c(attr(recurs, "id"), attr(noRed, "id")),
                        function(...) unique.default(unlist(...)),
                        recursive = FALSE, use.names = FALSE),
            n.asf = c(attr(recurs, "n.asf"), attr(noRed, "n.asf")))
}

# Generic function minimalizeCsf
minimalizeCsf <- function(x, ...){
  UseMethod("minimalizeCsf")
}
# Default method (for character vector)
minimalizeCsf.default <- function(x, ct = full.ct(x), verbose = FALSE, ..., data){
  
    # Ensure backward compatibility of argument data
    if (!missing(data)){
      warning("Argument 'data' is deprecated in minimalizeCsf(); use 'ct' instead.", 
              call. = FALSE)
      if (missing(ct)) ct <- data
    }

  if (length(x) == 0){
    return(emptyMinimalizeCsf())
  }
  x <- noblanks(x)
  ct <- configTable(ct, verbose = verbose, 
                    rm.dup.factors = FALSE, rm.const.factors = FALSE)
  minim <- .minCsf(x, ctInfo(ct), verbose = verbose)
  redundantParts <- vector("list", length(minim))
  for (i in seq_along(minim)){
    redasf <- lapply(extract_asf(x[attr(minim, "id")[[i]]]), 
                     setdiff, extract_asf(minim[[i]])[[1]])
    redasf <- happly(redasf, function(x) if (length(x)) paste0("(", x, ")") else character(0))
    redundantParts[[i]] <- C_mconcat(redasf, "*")
  }
  out <- qcondTbl_csf(minim, 
                      ctInfo(ct)$scores, 
                      attr(ct, "n"))
  out$n.asf <- attr(minim, "n.asf")
  out$redundantParts <- redundantParts 
  class(out) <- c("minimalizeCsf", "data.frame")
  out
}
# Method for class cna
minimalizeCsf.cna <- function(x, n = 20, verbose = FALSE, ...){
  csfs <- csf(x, n.init = max(n), inus.only = FALSE)$condition
  if (length(csfs) == 0){
    return(emptyMinimalizeCsf())
  }
  if (length(n)>1){
    n <- n[n<length(csfs)]
    csfs <- csfs[n]
  }
  minimalizeCsf.default(csfs, ct = x$configTable, verbose = verbose, ...)
}

# print method
print.minimalizeCsf <- function(x, subset = 1:5, ...){
  n <- nrow(x)
  subset <- intersect(subset, seq_len(n))
  cat("Object of class 'minimalizeCsf' containing", n, "solution(s)\n")
  if (length(subset) < n) 
    cat("Printing a subset of", length(subset), "solution(s)\n")
  cat("\n")
  for (i in subset){
    cat("=== Solution ", i, " ===\nCondition:\n", sep = "")
    writeLines(c(x$condition[i], ""))
    print.data.frame(x[i, -match(c("condition", "redundantParts"), names(x))])
    writeLines(c("\nredundant parts:",
                 unlist(x[i, 6])))
    cat("\n")
  }
  invisible()
}

# Output object for input with zero csf
emptyMinimalizeCsf <- function() structure(list(
  outcome = c("JSR", "JSR", "JSR", "JSR"), condition = character(0), 
  con = numeric(0), 
  cov = numeric(0), 
  n.asf = integer(0),
  redundantParts = vector("list", 0)), 
  class = c("minimalizeCsf", "data.frame")
)