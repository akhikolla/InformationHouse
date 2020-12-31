#' Print CGGP object
#' 
#' Default print as a list is bad since there's a lot of elements.
#'
#' @param x CGGP object
#' @param ... Passed to print
#'
#' @return String to be printed
#' @export
#'
#' @examples
#' SG = CGGPcreate(3,21)
#' print(SG)
#' f <- function(x) {x[1]+exp(x[2]) + log(x[3]+4)}
#' y <- apply(SG$design, 1, f)
#' SG <- CGGPfit(SG, y)
#' print(SG)
print.CGGP <- function(x, ...) {
  s <- paste0(
    c(
      "CGGP object\n",
      "   d = ", x$d, '\n',
      "   output dimensions = ", if (is.matrix(x$Y)) ncol(x$Y) else {1}, '\n',
      "   CorrFunc = ", x$CorrName, '\n',
      "   number of design points             = ", if (is.null(x$design) || length(x$design)==0) {"0"} else {nrow(x$design)}, '\n',
      "   number of unevaluated design points = ", if (is.null(x$design_unevaluated)) 0 else nrow(x$design_unevaluated), '\n',
      if (is.null(x$Xs)) {""} else {paste0("   number of supplemental points       = ", nrow(x$Xs), '\n')},
      "   Available functions:\n",
      "     - CGGPfit(CGGP, Y) to update parameters with new data\n",
      "     - CGGPpred(CGGP, xp) to predict at new points\n",
      "     - CGGPappend(CGGP, batchsize) to add new design points\n",
      "     - CGGPplot<name>(CGGP) to visualize CGGP model\n"
    )
  )
  cat(s, sep="")
}

#' S3 predict method for CGGP
#' 
#' Passes to CGGPpred
#' 
#' @param object CGGP object
#' @param ... Other arguments passed to `CGGPpred`
#'
#' @rdname CGGPpred
#' @export
predict.CGGP <- function(object, xp, ...) {
  CGGPpred(CGGP=object, xp=xp, ...)
}

#' S3 plot method for CGGP
#' 
#' There are a few different plot functions for CGGP objects:
#' `CGGPplotblocks`, `CGGPplotblockselection`, 
#' `CGGPplotcorr`, `CGGPplotheat`, `CGGPplothist`,
#' `CGGPvalplot`, 
#' `CGGPplotslice`, `CGGPplotslice`, and `CGGPplotvariogram`.
#' Currently `CGGPplotblocks` is the default plot object.
#'
#' @param x CGGP object
#' @param y Don't use
#' @param ... Passed to CGGPplotblocks
#'
#' @return Either makes plot or returns plot object
#' @export
#'
#' @examples
#' SG = CGGPcreate(3,100)
#' plot(SG)
plot.CGGP <- function(x, y, ...) {
  CGGPplotblocks(x, ...)
}