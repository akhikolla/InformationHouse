#' Determine zones for flexibly shaped spatial scan test
#'
#' \code{rflex_zones} determines the unique zones to
#' consider for the flexibly shaped spatial scan test of
#' Tango and Takahashi (2012).  The algorithm uses a
#' breadth-first search to find all subgraphs connected to
#' each vertex (region) in the data set of size \eqn{k} or
#' less with the constraint that the middle p-value of each
#' region must be less than \code{alpha1}.
#'
#' @inheritParams rflex.test
#' @param nn An n by k matrix providing the k nearest
#'   neighbors of each region, presumably produced by the
#'   \code{\link{knn}} function.
#' @param pop The population size associated with each
#'   region.  The default is \code{NULL} since this argument
#'   is only needed for \code{type = "binomial"}.
#' @param loop A logical value indicating whether a loop
#'   should be used to implement the function instead of
#'   \code{\link[pbapply]{pbapply}}.  The default is
#'   \code{FALSE}. If \code{TRUE}, then memory-saving steps
#'   are also taken.
#' @param verbose A logical value indicating whether
#'   progress messages should be provided.
#'   The default is \code{FALSE}.  If both \code{loop} and
#'   \code{verbose} are \code{TRUE}, informative messages
#'   are displayed that can be useful for diagnosing where
#'   the sequences of connected subgraphs are slowing down
#'   or having problems.
#' @param pfreq The frequency that messages are reported
#'   from the loop (if \code{verbose = TRUE}). The default
#'   is \code{pfreq = 1}, meaning a message is returned for
#'   each index of the loop.
#' @return Returns a list of zones to consider for
#'   clustering.  Each element of the list contains a vector
#'   with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Tango, T. and Takahashi, K. (2012), A
#'   flexible spatial scan statistic with a restricted
#'   likelihood ratio for detecting disease clusters.
#'   Statist. Med., 31: 4207-4218. <doi:10.1002/sim.5478>
#' @seealso rflex.midp
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$x, nydf$y)
#' nn = knn(coords, k = 5)
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ex = pop * sum(cases)/sum(pop)
#' # zones for poisson model
#' pzones = rflex_zones(nn, w = nyw, cases = cases, ex = ex)
#' \dontrun{
#' pzones = rflex_zones(nn, w = nyw, cases = cases,
#'                       ex = ex, verbose = TRUE)
#' # zones for binomial model
#' bzones = rflex_zones(nn, w = nyw, cases = cases, ex = ex,
#'                      type = "binomial", pop = pop)
#' }
rflex_zones = function(nn, w, cases, ex, alpha1 = 0.2,
                       type = "poisson", pop = NULL,
                       cl = NULL, loop = FALSE,
                       verbose = FALSE, pfreq = 1) {
  arg_check_rflex_zones(nn, w, cases, ex, alpha1, type, pop,
                        loop, pfreq, verbose)

  # compute mid p-value
  p = rflex.midp(cases, ex, type = type, pop = pop)
  # determine which regions are "hot" (keep) or "cold" (remove)
  keep = which(p < alpha1)
  nkeep = length(keep)
  k = max(sapply(nn, length))

  if (length(keep) > 0) {
    remove = setdiff(seq_along(ex), keep)

    # remove connections when p >= alpha1
    w[, remove] = 0

    N = nrow(w)
    idx = keep
    lprimes = log(randtoolbox::get.primes(N))

    if (!loop) {
      # get list of list of logical vectors
      czones = scsg2_cpp(nn, w, idx = idx, nlevel = k, verbose = verbose, lprimes)
      # convert to zone indices
      czones = logical2zones(czones, nn, idx)
      # return distinct zones
      return(czones[distinct(czones)])
    } else {
      czones = list()
      pri = randtoolbox::get.primes(N)
      czones_id = numeric(0) # unique identifier of each zone
      for (i in keep) {
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
  } else {
    czones = vector("list", 1)
    czones[[1]] = numeric(0)
    return(czones)
  }
}

#' Check arguments of rflex.zones
#'
#' Check the arguments of rflex.zones
#'
#' @return NULL
#' @noRd
arg_check_rflex_zones = function(nn, w, cases, ex,
                                 alpha1, type, pop,
                                 loop, pfreq, verbose) {
  if (is.list(dim(nn))) stop("nn must be a list of nn vectors")
  N = length(nn)
  if (nrow(w) != N) stop("nrow(w) must match length(nn)")
  if (length(cases) != N) stop("length(cases) must match length(nn)")
  if (length(alpha1) != 1 | alpha1 <= 0) {
    stop("alpha1 must be in (0, 1]")
  }
  if (!is.element(type, c("poisson", "binomial"))) {
    stop("type must be 'poisson' or 'binomial'")
  }
  if (!is.null(pop)) {
    if (length(pop) != N) {
      stop("length(pop) != length(nn)")
    }
  }
  if (length(loop) != 1 | !is.logical(loop)) {
    stop("loop must be a logical value")
  }
  if (length(pfreq) != 1 | pfreq < 1) {
    stop("pfreq must be a positive integer")
  }
  if (length(verbose) != 1 | !is.logical(verbose)) {
    stop("verbose must be a logical value")
  }
}
