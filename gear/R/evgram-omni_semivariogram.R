#' Compute an omnidirectional empirical semivariogram
#'
#' @param d Vector of pairwise distances
#' @param diff Vector of pairwise response differences
#' @param cuts Cut points for bins
#' @param npmin Minimum pairs of differences
#' @param type Type of empirical semivariogram
#' @keywords internal
#' @return A data.frame with the computed omnidirectional semivariogram
omni_semivariogram = function(d, diff, cuts, npmin, type) {
  np = unlist(lapply(split(d, cuts), length), use.names = FALSE)

  if (type == "standard") {
    semivariance = unlist(lapply(split(diff^2, cuts), mean), use.names = FALSE)/2
  } else if (type == "cressie") {
    semivariance = unlist(lapply(split(sqrt(abs(diff)), cuts), mean), use.names = FALSE)^4/2/(0.457 + 0.494/np)
  } else if (type == "cloud") {
    non_na = !is.na(cuts)
    return(data.frame(distance = d[non_na], semivariance = diff[non_na]^2/2, np = 1))
  }

  distance  = unlist(lapply(split(d, cuts), mean), use.names = FALSE)
  which_small = which(np < npmin)
  if (length(which_small) > 0) {
    return(data.frame(distance = distance[-which_small],
                      semivariance = semivariance[-which_small],
                      np = np[-which_small]))
  } else {
    return(data.frame(distance = distance,
                      semivariance = semivariance,
                      np = np))
  }
}
