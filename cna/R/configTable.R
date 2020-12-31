
# New version of configTable
# ==========================
# Modifications:
#   - types "cs", "mv", und "fs"
#   - Behaviour if x already is a configTable
configTable <- function(x, type = c("cs", "mv", "fs"), frequency = NULL,
                        case.cutoff = 0, rm.dup.factors = TRUE, rm.const.factors = TRUE,
                        .cases = NULL, verbose = TRUE){
  nm.x <- deparse(substitute(x))
  stopifnot(is.data.frame(x) || is.matrix(x), length(case.cutoff) == 1)
  if (is.matrix(x)) {
    if (is.null(colnames(x)))
      colnames(x) <- make.unique(rep(LETTERS, length.out = ncol(x)))
    x <- as.data.frame(x)
  } else if (inherits(x, "configTable")){
    if (!is.null(frequency))
      stop(nm.x, " already is a configTable - argument frequency is not admissable in this case!")
    frequency <- attr(x, "n", exact = TRUE)
    .cases <- attr(x, "cases")
    type <- attr(x, "type")
    x <- as.data.frame(x)
  }
  type <- tolower(type)
  type <- match.arg(type)
  if (is.null(.cases)){
    .cases <- try(as.character(rownames(x)))
    #if (inherits(.cases, "try-error")){
    if (anyNA(iconv(.cases, "", "ASCII"))){
      warning("The row names contain special (non-ASCII) characters and are therefore not used.")
      .cases <- seq_len(nrow(x))
    }
    .cases <- as.list(.cases)
  } else {
    if (length(.cases) != nrow(x)) stop("length(.cases) must be the same as nrow(x)")
    if (is.atomic(.cases)) .cases <- as.character(.cases)
  }
    
  # check if data entries are correct
  datatype <- switch(type, cs = "integer",
                           mv = "integer",
                           fs = "numeric",
                           mv2cs = "integer")
  xvalues <- unlist(x, use.names = FALSE, recursive = FALSE)
  if (type == "cs" & !all(xvalues %in% 0:1))
    stop("Only values 0 and 1 are allowed in x with type=\"cs\".")
  if (type == "mv" & !all(xvalues %% 1 == 0 & xvalues >= 0))
    stop("Only integer values >=0 are allowed in x with type=\"mv\".")
  if (type == "fs" & !all(xvalues >= 0 & xvalues <= 1))
    stop("Only values >=0 and <=1 are allowed in x with type=\"fs\".")
  for (i in seq_along(x)) mode(x[[i]]) <- datatype
  cx <- do.call(paste, c(x, list(sep = "\r")))
  cx <- factor(cx, levels = unique(cx))
  ct <- x[!duplicated(cx), , drop = FALSE]
  rownames(ct) <- NULL
  if (is.null(frequency)){
    f <- tabulate(cx)
  } else {
    # frequencies
    if (length(frequency) != nrow(x) || any(frequency < 0) || any(is.na(frequency)))
      stop("Inadmissible frequency argument")
    f <- as.integer(as.vector(rowsum(frequency, cx)))
  }
  .cases <- lapply(split(.cases, as.integer(cx)), unlist, recursive = FALSE, use.names = FALSE)
  ll <- lengths(.cases)
  if (!isTRUE(all.equal(f, ll, check.attributes = FALSE)) && length(ll) > 0){
    .cases <- mapply(rep_len, .cases, f, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    .cases <-
      split.default(make.unique(unlist(.cases, use.names = FALSE)),
                    rep(seq_along(.cases), f))
  }
  # case.cutoff
  del.cases <- f < case.cutoff
  if (any(del.cases) && verbose){
    message(sum(f[del.cases]), " of ", sum(f),
        " cases are removed due to case.cutoff = ", case.cutoff,
        ".")
    ct <- ct[!del.cases, , drop = FALSE]
    f <- f[!del.cases]
    .cases <- .cases[!del.cases]
  }
  # remove rows with frequency 0
  if (any(rm.rows <- (f == 0L))) {
    ct <- ct[!rm.rows, , drop = FALSE]
    f <- f[!rm.rows]
  }
  # constant columns
  constcols <- vapply(ct, isConstant, logical(1))
  if (any(constcols)){
    nms.constcols <- names(ct)[constcols]
    if (rm.const.factors){
      ct <- ct[, !constcols, drop = FALSE]
      if (verbose) message("The following factors are constant and therefore eliminated: ",
                           C_concat(nms.constcols, sep = ", "))
    } else if (verbose){
    message("The following factors are constant: ",
            C_concat(nms.constcols, sep = ", "))
    }        
  }          
  # eliminate duplicated columns
  dupcols <- duplicated(t(ct))
  if (any(dupcols)){
    nms.dupcols <- names(ct)[dupcols]
    if (rm.dup.factors){
      ct <- ct[, !dupcols, drop = FALSE]
      if (verbose)
        message("The following factors are duplicates of factors in preceding columns of the data table: ",
                C_concat(nms.dupcols, sep = ", "),
                "\nThey are eliminated (because rm.dup.factors = TRUE).")
    } else if (verbose){
    message("The following factors are duplicates of factors in preceding columns of the data table: ",
            C_concat(nms.dupcols, sep = ", "))
    } 
  }
  # Warn if names are not syntactically valid
  nms <- names(ct) <- toupper(names(ct))
  if (!identical(nms, make.names(nms, unique = TRUE))){
    warning("configTable has syntactically invalid names. condition(), cna() and other functions may not work.",
            call. = FALSE)
  }
  # output
  class(ct) <- c("configTable", "data.frame", "truthTab")
  attr(ct, "n") <- as.vector(f)
  attr(ct, "cases") <- lapply(.cases, function(x) sort(unlist(x, use.names = FALSE, recursive = FALSE)))
  attr(ct, "type") <- type
  ct
}

# print-method for configTable
# ============================
# New: Type is displayed
# cases now shown as rownames
#   Default show.cases: TRUE, if all rownames shorter than 20 characters
#   ...  digits, quote, right passed to print.data.frame
print.configTable <- function(x, show.cases = NULL, ...){ 
  # check attribute 'n'
  attr.n <- attr(x, "n", exact = TRUE)
  n.ok <- is.numeric(attr.n) && length(attr.n) == nrow(x)
  if (is.null(attr.n)){
    warning("Attribute \"n\" is missing")
  } else if (!n.ok){
    warning("Number of case frequencies does not match nrow(x) - all case frequencies are set to 1.")
  } 
  if (!n.ok) attr(x, "n") <- rep(1L, nrow(x))
  # check attribute 'cases'
  attr.cases <- attr(x, "cases", exact = TRUE)
  cases.ok <- length(attr.cases) == nrow(x)
  if (is.null(attr.cases)){
    warning("Attribute \"cases\" is missing")
  } else if (!cases.ok){
    warning("Number of case labels does not match nrow(x) - case labels are ignored.")
    attr(x, "cases") <- NULL
  } else if (max(lengths(attr(x, "cases")))>100){
    # message("Case labels are suppressed.")
    cases.ok <- FALSE
  }
  # create data frame for printing
  df.args <- list(as.data.frame(x), n.obs = attr(x, "n", exact = TRUE), check.names = FALSE)
  ncols <- length(df.args[[1]])
  prntx <- do.call(data.frame, df.args)
  if (cases.ok){
    rownames(prntx) <- C_mconcat(lapply(attr(x, "cases"), sort), sep = ",")
  }
  cat(paste0("configTable of type \"", attr(x, "type"), "\"\n"))
  if (is.null(show.cases) && nrow(x) > 0){
    show.cases <- max(nchar(rownames(prntx))) <= 20
    labels.suppressed <- !show.cases
  }  
  printDfwithBars(prntx, row.names = show.cases, ..., bar.after = ncols)
  cat("Total no.of.cases:", sum(attr(x, "n", exact = TRUE)), "\n")
  if (exists("labels.suppressed", inherits = FALSE) && labels.suppressed)
    message("Printing of case labels is suppressed because some labels are longer than 20 characters.\n",
            "Use print with show.cases=TRUE to display them.")
  invisible(x)
}

# "["-method for configTable:
# ===========================
#  Note that default of argument 'drop' differs from that in "[.data.frame"
`[.configTable` <- function(x, i, j, drop = FALSE, rm.dup.factors = FALSE, rm.const.factors = FALSE){
  cl <- sys.call()
  cl$rm.dup.factors <- cl$rm.const.factors <- NULL
  len <- length(cl)
  nms <- names(cl)
  if ("drop" %in% nms && !isFALSE(drop)){
    warning("'drop=TRUE' is ignored")
    cl$drop <- NULL
    len <- len - 1L
  }
  noRowPos <- len == 3
  cl1 <- quote(`[.data.frame`(x, , , drop = FALSE))
  cl1[[2]] <- call("as.data.frame", cl[[2]])
  cl1[[3]] <- if (noRowPos || identical(cl[[3]], quote(expr = ))) TRUE else cl[[3]]
  cl1[[3]] <- i <- eval.parent(cl1[[3]])
  cl1[[4]] <- if (noRowPos) cl[[3]] else cl[[4]]
  xx <- eval.parent(cl1)
  configTable(xx,
              type = attr(x, "type"),
              frequency = attr(x, "n", exact = TRUE)[i],
              .cases = attr(x, "cases")[i],
              rm.dup.factors = rm.dup.factors, rm.const.factors = rm.const.factors, verbose = FALSE)
}


# More methods for class "configTable"
# ====================================
as.data.frame.configTable <- function(x, ...){
  class(x) <- "data.frame"
  attributes(x)[c("type", "cases", "n")] <- NULL
  x
}
`[<-.configTable` <- function(x, i, j, value){
  out <- as.data.frame(NextMethod("[<-", x))
  row.sel <- nargs() == 5 && !missing(i)
  rows <- if (row.sel) i else if (nrow(x)) TRUE else logical(0)
  out <- configTable(out,
           type = attr(x, "type"),
           frequency = attr(x, "n", exact = TRUE),
           .cases = attr(x, "cases")[rows],
           rm.dup.factors = FALSE, rm.const.factors = FALSE, verbose = FALSE)
}

`$<-.configTable` <- function(x, name, value){
  out <- as.data.frame(NextMethod("$<-", x))
  out <- configTable(out,
           type = attr(x, "type"),
           frequency = attr(x, "n", exact = TRUE),
           .cases = attr(x, "cases"),
           rm.dup.factors = FALSE, rm.const.factors = FALSE, verbose = FALSE)
}

`[[<-.configTable` <- function(x, i, j, value){
  out <- as.data.frame(NextMethod("[[<-", x))
  out <- configTable(out,
           type = attr(x, "type"),
           frequency = attr(x, "n", exact = TRUE),
           .cases = attr(x, "cases"),
           rm.dup.factors = FALSE, rm.const.factors = FALSE, verbose = FALSE)
}

subset.configTable <- function(x, subset, ...){
  r <- if (missing(subset))
    rep_len(TRUE, nrow(x))
  else {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r))
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  x[r, , ...]
}

head.configTable <- function (x, n = 6L, ...) 
{
  stopifnot(length(n) == 1L)
  n <- if (n < 0L) 
    max(nrow(x) + n, 0L)
  else min(n, nrow(x))
  x[seq_len(n), , ...]
}

tail.configTable <- function (x, n = 6L, ...) 
{
  stopifnot(length(n) == 1L)
  nrx <- nrow(x)
  n <- if (n < 0L) 
    max(nrx + n, 0L)
  else min(n, nrx)
  x[seq.int(to = nrx, length.out = n), , ...]
}

# combine:
# ========
# Neu ist der Parameter sep
# combine <- function (x, rows = rownames(x), sep = ""){
#   rows <- do.call(paste, c(x, sep = sep))
#   factor(rows, levels = unique(rows))
# }

################################################################################

# versions of configTable with fixed type
csct <- function(...){
  cl <- match.call(configTable, sys.call())
  stopifnot(is.null(cl$type))
  cl[[1]] <- quote(configTable)
  cl$type <- "cs"
  eval.parent(cl)
}
mvct <- function(...){
  cl <- match.call(configTable, sys.call())
  stopifnot(is.null(cl$type))
  cl[[1]] <- quote(configTable)
  cl$type <- "mv"
  eval.parent(cl)
}
fsct <- function(...){
  cl <- match.call(configTable, sys.call())
  stopifnot(is.null(cl$type))
  cl[[1]] <- quote(configTable)
  cl$type <- "fs"
  eval.parent(cl)
}
