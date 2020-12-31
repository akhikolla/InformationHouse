#' Color clusters
#'
#' \code{color.clusters} is a helper function to color
#' clusters of regions produced by an appropriate method,
#' e.g., \code{scan.test} or \code{uls.test}.  Regions that
#' are not part of any cluster have no color.
#'
#' @param x An object of class scan produced by a function
#'   such as \code{scan.test}.
#' @param col A vector of colors to color the clusters in
#'   \code{x}.  Should have same length as the number of
#'   clusters in \code{x}.
#' @return Returns a vector with colors for each
#'   region/centroid for the data set used to construct
#'   \code{x}.
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, alpha = 0.12, longlat = TRUE,
#'                 nsim = 9)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
color.clusters = function(x, col = grDevices::hcl.colors(length(x$clusters))) {
  if (class(x) != "scan" & class(x) != "smerc_cluster") {
    stop("x should be an object of class scan or smerc_cluster.")
  }
  if (length(x$clusters) != length(col)) {
    stop("The number of colors must match the number of clusters.")
  }

  mycol = numeric(nrow(x$coords))
  for (i in seq_along(x$clusters)) {
    mycol[x$clusters[[i]]$loc] = col[i]
  }
  return(mycol)
}
