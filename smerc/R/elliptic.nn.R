#' Nearest neighbors for elliptic scan
#'
#' \code{elliptic.nn} computes the nearest neighbors
#' relationships for \code{elliptic.test}.  It will provide
#' a list of nearest neighbors, and a list of the associated
#' shape and angle.
#'
#' @inheritParams elliptic.test
#'
#' @return A list of nested nearest neighbors, the
#'   associated shapes and angles for each set of nn, and
#'   all of the shapes and angles you get for each zone
#'   constructed from the set of nearest neighbors.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' enn = elliptic.nn(coords, nydf$pop, 0.1,
#'                   shape = c(1, 1.5), nangle = c(1, 4))
elliptic.nn = function(coords, pop, ubpop = 0.5,
                       shape = c(1, 1.5, 2, 3, 4, 5),
                       nangle = c(1, 4, 6, 9, 12, 15)) {

  if (is.null(dim(coords))) {
    stop("coords must be matrix-like")
  }
  if (length(pop) != nrow(coords)) {
    stop("length(pop) != nrow(coords)")
  }
  if (length(ubpop) != 1 | ubpop <= 0 | ubpop > 1) {
    stop("ubpop must be in (0, 1]")
  }
  if (length(shape) != length(nangle)) {
    stop("length(shape) != length(nangle)")
  }

  coords = as.matrix(coords)

  # determine distances for each angle/shape
  d = lapply(seq_along(shape), function(i) {
    all_shape_dists(shape[i], nangle[i], coords)
  })
  # bind distances vertically
  if (length(shape) > 1) {
    d = do.call(rbind, d)
  } else {
    d = d[[1]]
  }

  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  nn = nnpop(d, pop, ubpop)

  # determine number of observations, length of zones
  N = nrow(coords)
  nnn = unlist(lapply(nn, length))

  # shape associated with each set of nn, all zones
  shape_nn =  rep(shape, times = N * nangle)
  shape_all = rep(shape_nn, times = nnn)

  # angle associated with each set of nn, all zones
  angle_nn = unlist(sapply(seq_along(nangle), function(i) {
    seq(90, 270, len = nangle[i] + 1)[seq_len(nangle[i])]
  }))
  angle_nn = rep(angle_nn, each = N)

  angle_all = rep(angle_nn, times = nnn)
  return(list(nn = nn, shape_nn = shape_nn, angle_nn = angle_nn,
              shape_all = shape_all, angle_all = angle_all))
}
