#' Run length encoding with indices
#'
#' @param x vector of values on which to perform run length encoding
#' @param zero_index logical. Should indices be indexed from zero (useful for
#'   Rcpp)?
#'
#' @return A data frame with information about the runs and start/stop indices
#' @export
#'
#' @examples
#'
#' x <- c(
#'   FALSE, TRUE, FALSE, FALSE, FALSE, TRUE,
#'   FALSE, TRUE, TRUE, FALSE, TRUE, FALSE,
#'   FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'   FALSE, TRUE
#' )
#' head(index_runs(x))
index_runs <- function(x, zero_index = FALSE) {

  result <- do.call(
    data.frame, rle(x)
  ) %>% {cbind(.,
    end_index = cumsum(.$lengths)
  )} %>% {cbind(.,
    start_index = .$end_index - .$lengths + 1
  )} %>% {
    .[ ,rev(names(.))]
  }

  if (!zero_index) return(result)

  sapply(
    result[ ,1:2], function(y) y - 1
  ) %>% {cbind(.,
    result[ ,3:4]
  )}

}
