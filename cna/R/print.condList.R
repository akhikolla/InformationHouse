
print.condList <- function (x, ...){
  nn <- names(x)
  ll <- length(x)
  if (length(nn) != ll)
    nn <- paste("Component", seq.int(ll))
  for (i in seq_len(ll)){
    x1 <- x[[i]]
    info <- attr(x1, "info")
    cat(addblanks(nn[i]), 
        if (info$reshaped) paste0(" [input: ", info$input, "]"),
        ":\n")
    print(x1, ...)
    cat("\n")
  }
  invisible(x)
}
`[.condList` <- function(x, ...){
  out <- NextMethod()
  attributes(out) <- c(
    attributes(out),
    attributes(x)[c("class", "type", "n", "cases", "ct")]
  )
  i <- eval.parent(sys.call()[[3]])
  if (is.character(i)) i <- match(i, names(x), 0L)
  attr(out, "info") <- attr(x, "info")[i, ]
  out
}
`[[.condList` <- function(x, ...){
  out <- NextMethod()
  if (identical(out, "Invalid condition")) return(out)
  attributes(out) <- c(
    attributes(out),
    attributes(x)[c("type", "n", "cases", "ct")])
  i <- eval.parent(sys.call()[[3]])
  if (is.character(i)) i <- match(i, names(x), 0L)
  attr(out, "info") <- attr(x, "info")[i, ]
  out
}
`$.condList` <- function(x, ...){
  out <- NextMethod()
  if (is.null(out) || identical(out, "Invalid condition")) return(out)
  attributes(out) <- c(
    attributes(out),
    attributes(x)[c("type", "n", "cases", "ct")])
  i <- sys.call()[[3]]
  i <- match(as.character(i), names(x), 0L)
  attr(out, "info") <- attr(x, "info")[i, ]
  out
}

summary.condList <- function(object, ...){
  print.condList(object, print.table = FALSE)
}
