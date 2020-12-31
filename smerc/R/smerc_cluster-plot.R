#' Plot object of class \code{smerc_cluster}.
#'
#' Plot clusters (the centroids of the regions in each
#' cluster) in different colors.  The most likely cluster is
#' plotted with solid red circles by default.  Points not in
#' a cluster are black open circles.  The other cluster
#' points are plotted with different symbols and colors.
#'
#' @param x An object of class scan to be plotted.
#' @param ... Additional graphical parameters passed to the
#'   \code{plot} function.
#' @param nclusters Number of clusters to plot.
#' @param ccol Fill color of the plotted points.  Default is
#' \code{grDevices::hcl.colors(nclusters, palette = "viridis")}.
#' @param cpch Plotting character to use for points in each
#'   cluster.  Default is NULL, indicating pch = 20 for the
#'   most likely cluster and then pch = 2, 3, .., up to the
#'   remaining number of clusters.
#' @param add A logical indicating whether results should be
#'   drawn on existing map.
#' @param usemap Logical indicating whether the maps::map
#'   function should be used to create a plot background for
#'   the coordinates.  Default is \code{FALSE}.  Use
#'   \code{TRUE} if you have longitude/latitude coordinates.
#' @param mapargs A list of arguments for the map function.
#' @export
#' @method plot smerc_cluster
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 0,
#'                 longlat = TRUE, alpha = 1)
#' plot(out, nclusters = 3)
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "county", region = "new york",
#'                xlim = range(out$coords[,1]),
#'                ylim = range(out$coords[,2]))
#' # needed for "county" database (unless you execute library(maps))
#' data(countyMapEnv, package = "maps")
#' plot(out, nclusters = 3, usemap = TRUE, mapargs = mapargs)
plot.smerc_cluster = function(x, ...,
                              nclusters = length(x$clusters),
                              ccol = NULL, cpch = NULL,
                              add = FALSE, usemap = FALSE,
                              mapargs = list()) {
  arg_check_nclusters(nclusters, length(x$clusters))
  # number of centroids
  nc = nclusters

  # set default values
  if (is.null(ccol)) {
    ccol = grDevices::hcl.colors(nclusters, palette = "viridis")
  }# seq_len(nclusters) + 1
  if (is.null(cpch)) cpch = rep(20, nclusters)
  if (nc > 1) cpch[2:nclusters] = 3:(nclusters + 1)

  # more sanity checking
  if (length(ccol) != nclusters) {
    stop("ccol must have length equal to nclusters")
  }
  if (length(cpch) != nclusters) {
    stop("cpch must have length equal to nclusters")
  }

  # extract coordinates and cluster coordinates
  coords = x$coords

  if (!add) {
    if (usemap) {
      if (requireNamespace("maps", quietly = TRUE)) {
        do.call(maps::map, mapargs)
      } else {
        message("Please install the maps package to enable this functionality")
        graphics::plot(coords, ...)
      }
    } else {
      graphics::plot(coords, ...)
    }
  }

  graphics::points(coords, ...)
  # plot clusters
  for (i in seq_len(nc)) {
    graphics::points(coords[x$clusters[[i]]$locids, 1],
                     coords[x$clusters[[i]]$locids, 2],
                     col = ccol[i], pch = cpch[i])

    if (!is.null(x$clusters[[i]]$w)) {
      if (length(x$clusters[[i]]$w) > 1) {
      seg = w2segments(x$clusters[[i]]$w,
                       coords[x$clusters[[i]]$locids, ])
      graphics::segments(seg[, 1], seg[, 2], seg[, 3], seg[, 4],
                         col = ccol[i], pch = cpch[i])
      }
    }
  }
}

#' Returns segments connecting neighbors
#'
#' w2segments takes a binary spatial proximity matrix and
#' associated coordinates.  It returns segments connecting
#' the neighbors for plotting purposes
#' @return NULL
#' @export
#' @keywords internal
w2segments = function(w, coords) {
  nnb = sum(w != 0)
  count = 0
  seg = matrix(0, nrow = nnb, ncol = 4)
  for (i in seq_len(nrow(w))) {
    nbi = which(w[i, ] != 0)
    nnbi = length(nbi)
    seg[count + seq_len(nnbi), 1] = coords[i, 1]
    seg[count + seq_len(nnbi), 2] = coords[i, 2]
    seg[count + seq_along(nbi), 3] = coords[nbi, 1]
    seg[count + seq_along(nbi), 4] = coords[nbi, 2]
    count = count + nnbi
  }
  return(seg)
}
