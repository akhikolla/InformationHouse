#' Determine zones for flexibly shaped spatial scan test
#'
#' \code{flex.zones} determines the unique zones to consider
#' for the flexibly shaped spatial scan test of Tango and
#' Takahashi (2005).  The algorithm uses a breadth-first
#' search to find all subgraphs connected to each vertex
#' (region) in the data set of size \eqn{k} or less.
#'
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
#' zones = flex.zones(coords, w = nyw, k = 3)
#' \dontrun{
#' # see what happens when verbose = TRUE
#' zones = flex.zones(coords, w = nyw, k = 3, verbose = TRUE)
#' }
flex.zones = function(coords, w, k = 10, longlat = FALSE,
                      cl = NULL, loop = FALSE,
                      verbose = FALSE, pfreq = 1) {
  nn = knn(coords = coords, longlat = longlat, k = k)
  N = nrow(coords)

  if (!loop) {
    fcall_list = list(X = seq_len(N), function(i, ...) {
      scsg2(nn, w, idx = i)
    }, cl = cl)

    # determine which apply function to use
    if (verbose) {
      message("constructing connected subgraphs:")
      fcall = pbapply::pblapply
      fcall_list = list(X = seq_len(N), function(i, ...) {
        scsg2(nn, w, idx = i)
      }, cl = cl)
      # determine zones
      czones = unlist(do.call(fcall, fcall_list),
                      use.names = FALSE, recursive = FALSE)
      return(czones[distinct(czones)])
    } else {
      return(scsg2(nn, w))
    }

    # determine zones
    czones = unlist(do.call(fcall, fcall_list),
                    use.names = FALSE, recursive = FALSE)
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
      # zones for idxi
      izones = scsg2(nn, w, idx = i)
      # determine unique ids for izones
      izones_id = sapply(izones, function(xi) sum(log(pri[xi])))
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
  # determine distinct zones
  # czones[distinct(czones)]
}
