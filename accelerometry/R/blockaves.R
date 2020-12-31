#' Block Averages
#' 
#' Calculates block averages (i.e. moving averages but for non-overlapping 
#' intervals) or maximum block average. For optimal speed, use 
#' \code{integer = TRUE} if \code{x} is an integer vector and 
#' \code{integer = FALSE} otherwise. If \code{length(x)} is not an exact 
#' multiple of \code{window}, the last partial segment is dropped.
#' 
#' @inheritParams movingaves
#' @inherit movingaves return
#' 
#' 
#' @examples
#' # Load accelerometer data for first 5 participants in NHANES 2003-2004
#' data(unidata)
#' 
#' # Get data from ID number 21005, Saturday only
#' counts.sat <- unidata[unidata[, "seqn"] == 21005 & unidata[, "paxday"] == 7, 
#'                       "paxinten"]
#'                       
#' # Calculate and plot hourly count averages
#' hourly.averages <- blockaves(x = counts.sat, window = 60, integer = TRUE)
#' plot(hourly.averages)
#' 
#' 
#' @export
blockaves <- function(x, window, integer = FALSE, max = FALSE) {
  
  # Call C++ function depending on 'integer' and 'max'
  if (integer) {
    if (! max) {
      return(.Call(`_accelerometry_blockaves_i`, x, window))
    }
    return(.Call(`_accelerometry_blockaves_i_max`, x, window))
  }
  if (! max) {
    return(.Call(`_accelerometry_blockaves_n`, x, window))
  }
  return(.Call(`_accelerometry_blockaves_n_max`, x, window))
  
}