#' Inverse Run Length Encoding (Alternate Implementation)
#' 
#' Re-constructs vector compressed by \code{\link{rle2}}.
#' 
#' 
#' @param x Object returned by \code{\link{rle2}}.
#' 
#' 
#' @return
#' Integer or numeric vector.
#' 
#' 
#' @examples
#' # Create dummie vector x
#' x <- c(0, 0, 0, -1, -1, 10, 10, 4, 6, 6)
#' 
#' # Summarize x using rle2
#' x.summarized <- rle2(x)
#' 
#' # Reconstruct x
#' x.reconstructed <- inverse_rle2(x.summarized)
#' 
#' 
#' @export
inverse_rle2 <- function(x) {
  return(rep(x[, "value"], x[, "length"]))
}
