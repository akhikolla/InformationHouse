#' Binary Distance of a Network Path
#' 
#' Calculates the binary distance of a user-specified network path through a network,
#' if all edges exist. Otherwise, returns \code{Inf} to signify infinite distance.
#'
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param path an integer vector of node indices from \code{sociomatrix}.
#' 
#' @examples
#' ## Calculate binary distance along a path in a sociomatrix
#' binary_distance(YangKnoke01, path = c(1,2,5))
#' 
#' ## This path doesn't exist
#' binary_distance(YangKnoke01, path = c(1,2,4,5))
#' 
#' @export
binary_distance <- function(sociomatrix, path){
  check_input(sociomatrix, path = path)
  return(1/gpv(sociomatrix, path, p = 0))
}

#' Peay's Path Value Measure
#'
#' Calculates path value as defined in Peay (1980). That is, returns the
#' value of the weakest connection in the path, if all edges exist. 
#' Otherwise, returns 0.
#'
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param path an integer vector of node indices from \code{sociomatrix}.
#' 
#' @examples
#' ## Calculate Peay's Path Value along a path in a sociomatrix
#' peay_path_value(YangKnoke01, path = c(1,2,5))
#' 
#' ## This path doesn't exist
#' peay_path_value(YangKnoke01, path = c(1,2,4,5))
#' 
#' @seealso \code{\link{peay_average_path_value}}
#' @export
peay_path_value <- function(sociomatrix, path){
  check_input(sociomatrix, path = path)
  return(gpv(sociomatrix, path, p = Inf))
}

#' Flament's Path Length Measure
#'
#' Calculates path length as defined in Flament (1963). That is, sums the
#' values of each edge in the path, if all edges exist. Otherwise, returns \code{NA}.
#'
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param path an integer vector of node indices from \code{sociomatrix}.
#'
#' @examples
#' ## Calculate Flament's Path Length along a path in a sociomatrix
#' flament_path_length(YangKnoke01, path = c(1,2,5))
#' 
#' ## This path doesn't exist
#' flament_path_length(YangKnoke01, path = c(1,2,4,5))
#'
#' @seealso \code{\link{flament_average_path_length}}
#' @export
flament_path_length <- function(sociomatrix, path){
  check_input(sociomatrix, path = path)
  path_sum <- 0
  nsteps <- length(path) - 1
  for(i in 1:nsteps){
    on <- path[i]
    to <- path[i + 1]
    if(sociomatrix[on,to] == 0){return(NA)}
    path_sum <- path_sum + sociomatrix[on, to]
  }
  return(path_sum)
}

#' Yang and Knoke's Average Path Value
#' 
#' Calculates 'APV' (Average Path Value) as defined in Yang, Knoke (2001)
#' Called \code{peay_average_path_value} in homage to E.R. Peay, who defined 
#' path length in 1980.
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param path an integer vector of node indices from \code{sociomatrix}.
#' 
#' @examples
#' ## Calculate 'APV' of a path in a sociomatrix
#' peay_average_path_value(YangKnoke01, path = c(1,2,5))
#' 
#' ## This path doesn't exist
#' peay_average_path_value(YangKnoke01, path = c(1,2,4,5))
#' 
#' @seealso \code{\link{peay_path_value}}
#' @export
peay_average_path_value <- function(sociomatrix, path){
  #peay_path_value() will handle input checking
  return(peay_path_value(sociomatrix, path)/(length(path) - 1))
}

#' Yang and Knoke's Average Path Length
#' 
#' Calculates 'APL' (Average Path Length) as defined in Yang, Knoke (2001).
#' Called \code{flament_average_path_length} in homage to A.C. Flament, who defined
#' path length in 1963.
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param path an integer vector of node indices from \code{sociomatrix}.
#' 
#' @examples
#' ## Calculate 'APL' of a path in a sociomatrix
#' flament_average_path_length(YangKnoke01, path = c(1,2,5))
#' 
#' ## This path doesn't exist
#' flament_average_path_length(YangKnoke01, path = c(1,2,4,5))
#' 
#' @seealso \code{\link{flament_path_length}}
#' @export
flament_average_path_length <- function(sociomatrix, path){
  #flament_path_length() will handle input checking
  return(flament_path_length(sociomatrix, path)/(length(path) - 1))
}