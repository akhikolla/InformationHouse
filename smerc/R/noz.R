#' Determine non-overlapping zones
#' 
#' Determine non-overlapping zones from a list of
#' candidate zones.
#' 
#' The function takes a list of candidate zones. Each
#' element of the list contains a candidate zones. The
#' candidate zones are defined by the location indices of
#' the regions comprising the zones.  Starting with the
#' first candidate zone, the function excludes every
#' candidate zone that intersects the first (any other
#' candidate zone that shares indices with the first zone).
#' Moving onto the next non-overlapping candidate zone,
#' the process is repeated. 
#' The function returns the indices (in the list of
#' zones) of the zones that do not overlap.
#' 
#' @param x A list containing the candidate zones.
#' @return A vector with the list indices of the 
#'   non-overlapping zones.
#' @author Joshua French
#' @export
#' @examples 
#' x = list(1:2, 1:3, 4:5, 4:6, 7:8)
#' noz(x)
noz = function(x)  {
  if (!is.list(x)) {
    stop("x should be a list with location ids for the regions in the clusters")
  }
  remain_idx = seq_along(x)
  i = 1
  u = 1
  while (i < length(x)) {
    no_inter = sapply(remain_idx, function(j) !any(duplicated(unlist(x[c(i, j)]))))
    remain_idx = remain_idx[no_inter]
    if (length(remain_idx) > 0) {
      i = min(remain_idx)
      u = c(u, i)
    }
    else {
      i = length(x) + 1
    }
  }
  return(u)
}
