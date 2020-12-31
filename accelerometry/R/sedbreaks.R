#' Sedentary Breaks
#' 
#' Identifies sedentary breaks in accelerometer count data. 
#' 
#' 
#' @inheritParams bouts
#' 
#' @param thresh Integer value specifying minimum count value to consider a 
#' break from sedentary time.
#' 
#' @param flags Logical value for whether to return a vector of 1's and 0's 
#' flagging the sedentary breaks (as opposed to the total number of sedentary 
#' breaks). 
#' 
#' 
#' @return
#' Integer value or vector depending on \code{flags}.
#' 
#' 
#' @examples
#' # Load accelerometer data for first 5 participants in NHANES 2003-2004
#' data(unidata)
#' 
#' # Get data from ID number 21005
#' id.part1 <- unidata[unidata[, "seqn"] == 21005, "seqn"]
#' counts.part1 <- unidata[unidata[, "seqn"] == 21005, "paxinten"]
#' 
#' # Identify periods of valid wear time
#' wear.part1 <- weartime(counts = counts.part1)
#' 
#' # Count number of sedentary breaks (over full week)
#' n.sedbreaks <- sedbreaks(counts = counts.part1, weartime = wear.part1)
#' 
#' # Flag sedentary breaks
#' sedbreaks.flagged <- sedbreaks(counts = counts.part1, weartime = wear.part1, 
#'                                flags = TRUE)
#' 
#' 
#' @export
sedbreaks <- function(counts, weartime = NULL, thresh = 100, flags = FALSE) {
  
  # If 'weartime' unspecified, make it a vector of 1's
  if (is.null(weartime)) {
    weartime <- as.integer(rep(1, length(counts)))
  }
  
  # Call C++ function depending on 'flags'
  if (flags) {
    out <- .Call(`_accelerometry_sedbreaks_flags`, counts, weartime, thresh)
  } else {
    out <- .Call(`_accelerometry_sedbreaks`, counts, weartime, thresh)
  }
  return(out)
  
}