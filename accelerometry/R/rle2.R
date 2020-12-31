#' Run Length Encoding (Alternate Implementation)
#' 
#' Summarizes vector containing runs of repeated values. Very similar to 
#' \code{\link[base]{rle}}, but sometimes much faster, and with an option to 
#' return the start/end indices for each run.
#' 
#' 
#' @param x Vector (see \code{class}).
#' @param class Character string specifying class of \code{x}. If unspecified, 
#' function figures it out (at cost of slightly slower run time).
#' @param indices Logical value for whether to record start/stop positions in 
#' addition to values and lengths for each run. 
#' 
#' @return
#' Integer or numeric matrix.
#' 
#' 
#' @examples
#' # Create dummie vector x
#' x <- c(0, 0, 0, -1, -1, 10, 10, 4, 6, 6)
#' 
#' # Summarize x using rle2
#' x.summarized <- rle2(x)
#' 
#' # Repeat, but also record start/stop indices for each run
#' x.summarized <- rle2(x = x, indices = TRUE)
#' 
#' 
#' @export
rle2 <- function(x, class = NULL, indices = FALSE) {
  
  # If class unspecified, figure it out
  if (is.null(class)) {
    class <- class(x)
  }
  
  # Call C++ function depending on 'class'
  if (class == "integer") {
    out <- .Call(`_accelerometry_rle2_i`, x, indices)
    colnames(out) <- c("values", "lengths")
  } else {
    out <- .Call(`_accelerometry_rle2_n`, x, indices)
  }
  
  # Add column names
  if (indices) {
    colnames(out) <- c("value", "start", "stop", "length")
  } else {
    colnames(out) <- c("value", "length")
  }
  return(out)
  
}
