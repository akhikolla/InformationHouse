
# print a data frame with vertical bars after certain columns
printDfwithBars <- function (x, digits = NULL, quote = FALSE, right = TRUE, 
                  row.names = TRUE, bar.after = integer(0), ...){
  n <- length(row.names(x))
  if (length(x) == 0L) {
    cat(sprintf(ngettext(n, "data frame with 0 columns and %d row", 
      "data frame with 0 columns and %d rows"), n), "\n", 
      sep = "")
  }
  else if (n == 0L) {
    print.default(names(x), quote = FALSE)
    cat(gettext("<0 rows> (or 0-length row.names)\n"))
  }
  else {
    m <- as.matrix(format.data.frame(x, digits = digits, 
      na.encode = FALSE))
    if (length(bar.after)){
      shift <- rowSums(outer(seq_along(x), bar.after, ">"))
      colPos <- seq_along(x) + shift
      mm <- matrix("|", nrow = n, ncol = ncol(x)+length(bar.after))
      colnames(mm) <- rep("|", ncol(mm))
      mm[, colPos] <- m
      colnames(mm)[colPos] <- colnames(m)
      rownames(mm) <- rownames(m)
      m <- mm
    }
    if (!isTRUE(row.names)) 
      dimnames(m)[[1L]] <- if (identical(row.names, FALSE)) 
        rep.int("", n)
      else row.names
    print(m, quote = quote, right = right, ...)
  }
  invisible(x)
}

# -----------------------------------------------------------------------------

# new print method for cond - with parameter add.data
print.cond <- function (x, digits = 3, print.table = TRUE, 
                        show.cases = NULL, add.data = NULL, ...){
  condType <- if (inherits(x, "booleanCond")) {
    "boolean"
  } else if (inherits(x, "atomicCond")) {
    "atomic"
  }
  cat("type of condition:", condType, "\n")
  if (print.table){
    dfr <- x
    if (!is.null(n.obs <- attr(x, "n", exact = TRUE)))
      dfr$n.obs <- n.obs
    rownames(dfr) <- if (is.null(attr(x, "cases"))){
      NULL
    } else {
      C_mconcat(attr(x, "cases"), sep = ",")
    }
    if (any(longName <- nchar(colnames(dfr))>23)){
      cnms <- colnames(dfr)
      cnms[longName] <- paste0(substring(cnms[longName], 1, 20), "...")
      colnames(dfr) <- cnms
    }    
    if (is.null(add.data) & !is.null(attr(x, "ct")))
      add.data <- attr(x, "ct")
    if (!is.null(add.data)){
      if (!inherits(add.data, "configTable")){
        add.data <- configTable(add.data, type = attr(x, "type"))
      }
      stopifnot(attr(x, "type") == attr(add.data, "type"), 
                isTRUE(all.equal(attr(add.data, "cases"), attr(x, "cases"))))
      dfr <- cbind(dfr, as.data.frame(add.data))
      bar.after <- 1:2 + (condType == "atomic")
    } else {
      bar.after <- 1 + (condType == "atomic")
    }
    if (is.null(show.cases))
      show.cases <- max(0L, nchar(rownames(dfr))) <= 20
    printDfwithBars(dfr, row.names = show.cases, ..., bar.after = bar.after)
  }
  info <- attr(x, "info")
  if (condType == "atomic"){
    cat("Consistency: ", 
        formatC(info[["consistency"]], format = "f", digits = digits), 
        " (", info[["sxy"]], "/", info[["sx"]], ")\n", 
        "Coverage:    ", 
        formatC(info[["coverage"]], format = "f", digits = digits), 
        " (", info[["sxy"]], "/", info[["sy"]], ")\n", 
        "Total no. of cases: ", info[["sumf"]], "\n", sep = "")
  }
  else 
    cat("Frequency:   ", 
        formatC(info[["freq"]], format = "f", digits = digits), 
        " (", info[["sy"]], "/", info[["sumf"]], ")\n", 
        sep = "")
  invisible(x)
}

# ------------------------------------------------------------------------------

# print method for complexCond - with parameter add.data
print.complexCond <- function (x, digits = 3, print.table = TRUE, 
                               show.cases = NULL, add.data = NULL, ...){
  cat("type of condition:", "complex\n")
  if (print.table) {
    dfr <- do.call(data.frame, c(x, list(check.names = FALSE)))
    names(dfr) <- unlist(lapply(x, names), use.names = FALSE)
    if (!is.null(n.obs <- attr(x, "n", exact = TRUE)))
      dfr$n.obs <- n.obs
    rownames(dfr) <- if (is.null(attr(x, "cases"))){
      NULL
    } else {
      C_mconcat(attr(x, "cases"), sep = ",")
    }
    if (is.null(add.data) & !is.null(attr(x, "ct")))
      add.data <- attr(x, "ct")
    if (!is.null(add.data)){
      if (!inherits(add.data, "configTable")){
        add.data <- configTable(add.data, type = attr(x, "type"))
      }
      stopifnot(attr(x, "type") == attr(add.data, "type"), 
                isTRUE(all.equal(attr(add.data, "cases"), attr(x, "cases"))))
      dfr <- cbind(dfr, as.data.frame(add.data))
      bar.after <- c(seq_along(x) * 2, length(x)*2 + 1)
    } else {
      bar.after <- seq_along(x) * 2
    }
    if (is.null(show.cases)) 
      show.cases <- max(0L, nchar(rownames(dfr))) <= 20
    printDfwithBars(dfr, row.names = show.cases, ..., bar.after = bar.after)
  }
  info <- attr(x, "info")
  cat("Consistency: ", formatC(info[["consistency"]], format = "f", digits = digits), "\n", 
      "Coverage:    ", formatC(info[["coverage"]], format = "f", digits = digits), "\n", 
      "Total no. of cases: ", sum(attr(x, "n", exact = TRUE)), "\n", 
      sep = "")
  cat("Consistency and coverage of components:\n")
  conCov <- cbind(info[["asfCons"]][[1]], info[["asfCovs"]][[1]])
  dimnames(conCov) <- list(names(x), c("consistency", "coverage"))
  print.default(conCov, digits = digits)
  invisible(x)
}

# ------------------------------------------------------------------------------


print.invalidCond <- function (x, ...){
  cat(x, " (", reason(x), ")\n", sep = "")
  invisible(x)
}

reason <- function(x){
  stopifnot(inherits(x, "invalidCond"))
  switch(attr(x, "info")$condTypes, 
         invalidValues = "non-existing values", 
         invalidSyntax = "invalid syntax")
}
