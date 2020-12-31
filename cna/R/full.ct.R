
# ==== full.ct() ====
# Generic function
full.ct <- function(x, ...) UseMethod("full.ct")

.expandUniqvals <- function(valueList){
  do.call(expand.grid, 
          c(valueList, list(KEEP.OUT.ATTRS = FALSE)))
}

# ==== Method for class 'cti' ====
# builds 
#   x       output of ctInfo configTable
# value:    cti of configTable, mv if original is mv, cs else
full.ct.cti <- function(x, ...){
  expanded <- .expandUniqvals(x$uniqueValues)
  x$scores <- do.call(cbind, lapply(x$resp_nms, factMat(x$type), 
                                      ct = expanded))
  if (x$type == "fs") x$type <- "cs"
  vId <- switch(x$type,
    cs = valueId <- 2 - as.matrix(expanded),
    mv = mapply(match, expanded, x$uniqueValues, SIMPLIFY = TRUE))
  nVal <- lengths(x$uniqueValues)
  vId <- sweep(vId, 
               2, 
               c(if (length(nVal)) 0, cumsum(nVal[-length(nVal)])), 
               "+")
  mode(vId) <- "integer"
  x$valueId <- vId
  x$freq <- rep(1L, nrow(expanded))
  x
}

# ==== Method for class 'configTable' ====
# builds full configTable
#   x       configTable
#   value:  configTable, mv if original is mv, cs else
full.ct.configTable <- function(x, ...){
  cti <- ctInfo(x)
  expanded <- .expandUniqvals(cti$uniqueValues)
  configTable(expanded, type = if (cti$type == "mv") "mv" else "cs",
              rm.dup.factors = FALSE, rm.const.factors = FALSE, verbose = FALSE)
}

# ==== Default Method (for matrix or data.frame) ====
# builds full.ct
# cases: character, list, matrix/data.frame
full.ct.default  <- function(x, type = c("cs", "mv", "fs"), ...){
  if (is.numeric(x) && length(x) == 1 && x %in% 0:length(LETTERS)){
    x <- LETTERS[seq_len(x)]
  }
  if (is.character(x)){
    type <- if (any(grepl("=", x))) "mv" else "cs"
    px <- lapply(x, tryparse)
    if (!all(ok <- !vapply(px, is.null, logical(1))))
      stop("Invalid input to full.ct:\n", paste0("  ", x[!ok], collapse = "\n"),
           call. = FALSE)
    if (type == "cs"){
      nms <- unique(toupper(unlist(lapply(px, all.vars))))
      x <- setNames(rep(list(0:1), length(nms)), nms)
    } else {
      vals <- rapply(px, .call2list, how = "unlist",
                     stopOps = c("==", "<", ">", "<=", ">="), 
                     validOps = c("<-", "<<-", "=", "&", "|", "(", "-"))
      vals <- vals[!vapply(vals, is.symbol, FUN.VALUE = TRUE)]
      if (any(lengths(vals) != 3)) stop("Check the conditions")
      vals <- unique.default(vals)
      var <- vapply(vals, function(x) as.character(x[[2]]), character(1))
      val <- vapply(vals, "[[", 3, FUN.VALUE = numeric(1))
      x <- lapply(split(val, var), sort)
    }
  }
  if (!is.data.frame(x) && is.list(x) && !is.null(names(x))){
    x <- .expandUniqvals(x)
    if (!all(unlist(x, use.names = FALSE) %in% 0:1))
      type <- "mv"
  }
  if (is.data.frame(x) || is.matrix(x)){
    return(full.ct.configTable(
      configTable(x, type = type, 
                  rm.dup.factors = FALSE, 
                  rm.const.factors = FALSE, 
                  verbose = FALSE, ...)))
  }
  stop("Don't know how to apply 'full.ct' to this input.")
}
  

# Aux function rowID
# generates a unique integer identifier of the rows in the ct
rowID <- function(cti){
  # if (cti$type == "fs") warning("Be careful when using rowID() with 'fs' data.")
  v <- matrix(unlist(cti$uniqueValues, use.names = FALSE)[cti$valueId], 
              nrow = nrow(cti$valueId), 
              dimnames = list(NULL, names(cti$nVal)))
  maxVal <- vapply(cti$uniqueValues, max, 1L) + 1L
  m <- c(1L, cumprod(maxVal[-length(maxVal)]))
  drop(v %*% m)
}
