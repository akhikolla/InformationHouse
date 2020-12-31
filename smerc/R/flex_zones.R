#' Determine zones for flexibly shaped spatial scan test
#'
#' \code{flex_zones} determines the unique zones to consider
#' for the flexibly shaped spatial scan test of Tango and
#' Takahashi (2005).  The algorithm uses a breadth-first
#' search to find all subgraphs connected to each vertex
#' (region) in the data set of size \eqn{k} or less.
#'
#' @param cl Ignored, but retained for backwards compatibility
#' @inheritParams flex.test
#' @inheritParams rflex.zones
#' @return Returns a list of zones to consider for
#'   clustering.  Each element of the list contains a vector
#'   with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Tango, T., & Takahashi, K. (2005). A flexibly
#'   shaped spatial scan statistic for detecting clusters.
#'   International journal of health geographics, 4(1), 11.
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$x, nydf$y)
#' zones = flex_zones(coords, w = nyw, k = 3)
#' \dontrun{
#' # see what happens when verbose = TRUE
#' zones = flex_zones(coords, w = nyw, k = 3, verbose = TRUE)
#' }
flex_zones = function(coords, w, k = 10, longlat = FALSE,
                       cl = NULL, loop = FALSE,
                       verbose = FALSE, pfreq = 1) {
  nn = knn(coords = coords, longlat = longlat, k = k)
  N = nrow(coords)
  idx = seq_along(nn)
  lprimes = log(randtoolbox::get.primes(N))

  if (!loop) {
    # get list of list of logical vectors
    czones = scsg2_cpp(nn, w, idx = idx, nlevel = k, lprimes = lprimes, verbose = verbose)
    # convert to zone indices
    czones = logical2zones(czones, nn, idx)
    # return distinct zones
    return(czones[distinct(czones)])
  } else {
    czones = list()
    pri = randtoolbox::get.primes(N)
    czones_id = numeric(0) # unique identifier of each zone
    for (i in seq_len(N)) {
      if (verbose) {
        if ((i %% pfreq) == 0) {
          message(i, "/", N, ". Starting region ", i,
                  " at ", Sys.time(), ".")
        }
      }
      # logical vector zones for idxi
      izones = scsg2_cpp(nn, w, i, k, lprimes, verbose = FALSE)
      # convert to region ids
      izones = logical2zones(izones, nn, idx = i)
      # determine unique ids for izones
      izones_id = sapply(izones, function(xi) sum(lprimes[xi]))
      # determine if some izones are duplicated with czones
      # remove duplicates and then combine with czones
      dup_id = which(izones_id %in% czones_id)
      if (length(dup_id) > 0) {
        czones = combine.zones(czones, izones[-dup_id])
        czones_id = c(czones_id, izones_id[-dup_id])
      } else {
        czones = combine.zones(czones, izones)
        czones_id = c(czones_id, izones_id)
      }
    }
    return(czones)
  }
}

#' Convert logical vector to zone
#'
#' @param czones List of list of logical vectors
#' @param nn List of nearest neighbor indices
#' @param idx Relevant nn indices
#' @return List of list of zones
#' @export
#' @keywords internal
logical2zones = function(czones, nn, idx = seq_along(nn)) {
  # for each element of czones,
  # strip the element (which is a list of logical vectors)
  # for each element of the list of logical vectors
  # get the nn for idx i and then subset with each logical vector
  unlist(lapply(seq_along(czones), function(i) {
    lapply(czones[[i]], function(x) {
      nn[[idx[i]]][x]
    })
  }), recursive = FALSE)
}

