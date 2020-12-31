#' Minimum spanning tree sequence
#'
#' \code{mst.seq} finds the sequence of connected regions
#' that maximize the spatial scan statistic (the likelihood
#' ratio test statistic) from a starting region. The set of
#' connected regions at each step is a candidate zone. The
#' zone continues to grow until no region should be added to
#' the zone due to relevant constraints (size, connectivity,
#' or other stopping criteria).  This function is not
#' intended to be used by users directly, but it can be
#' quite educational for seeing the spread of the cluster.
#' Consequently, it prioritizes efficiency over user
#' friendliness.
#'
#' The function can be used to construct candidate zones for
#' the dynamic minimum spanning tree (dmst), early stopping
#' dynamic minimum spanning tree (edmst), double connection
#' spatial scan test (dc), and maximum linkage spatial scan
#' test (mlink).
#'
#' \code{type} is a character vector indicating what should
#' be returned by the function.  If \code{type = "maxonly"},
#' then only the maximum of the log likelihood ratio test
#' statistic across all candidate zones is returned.  If
#' \code{type = "pruned"},, the function returns a list that
#' includes the location ids, test statistic, total cases,
#' expected cases, and total population for the zone with
#' the maximum test statistic.  It \code{type = "all"}, the
#' same information the same information is returned for the
#' entire sequence of zones.
#'
#' If \code{nlinks = "one"}, then a region only needs to be
#' connected to one other region in the current zone to be
#' considered for inclusion in the next zone.  If
#' \code{nlinks = "two"}, then the region must be connected
#' to at least two other regions in the current zone.  If
#' \code{nlinks = "max"}, then only regions with the maximum
#' number of connections to the current zone are considered
#' for inclusion in the next zone.
#'
#' @param start The initial region to start the candidate
#'   zone.
#' @param neighbors A vector containing the neighbors for
#'   the starting region (in ascending order of distance
#'   from the region).  The staring region itself is
#'   included among the neighbors.
#' @inheritParams scan.test
#' @inheritParams scan.stat
#' @inheritParams uls.test
#' @param max_pop The population upperbound (in total
#'   population) for a candidate zone.
#' @param type One of \code{"maxonly"}, \code{"pruned"}, or
#'   \code{"all"}.  The default is \code{"maxonly"}.  See
#'   Details.
#' @param nlinks A character vector.  The options are
#'   \code{"one"}, \code{"two"}, or \code{"max"}.  See
#'   Details.
#' @param early A logical value indicating whether the
#'   "early" stopping criterion should be used.  If
#'   \code{TRUE}, the sequence is stopped when the next
#'   potential zone doesn't produce a test statistic larger
#'   than the current zone. The default is \code{FALSE}.
#' @return Returns a list of relevant information.  See
#'   Details.
#' @author Joshua French
#' @export
#' @examples
#' # load data
#' data(nydf)
#' data(nyw)
#'
#' # create relevant data
#' coords = nydf[,c("longitude", "latitude")]
#' cases = floor(nydf$cases)
#' pop = nydf$population
#' w = nyw
#' ex = sum(cases)/sum(pop)*pop
#' ubpop = 0.5
#' ubd = 0.5
#' ty = sum(cases)   # total number of cases
#' # intercentroid distances
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' # upperbound for population in zone
#' max_pop = ubpop * sum(pop)
#' # upperbound for distance between centroids in zone
#' max_dist = ubd * max(d)
#' # create list of neighbors for each region (inclusive of region itself)
#' all_neighbors = nndist(d, ubd)
#' # find the dmst max zone
#' mst.seq(start = 1, all_neighbors[[1]], cases, pop, w, ex,
#'         ty, max_pop)
#' mst.seq(start = 1, all_neighbors[[1]], cases, pop, w, ex,
#'         ty, max_pop, "pruned")
#' bigout = mst.seq(start = 1, all_neighbors[[1]], cases, pop,
#'                  w, ex, ty, max_pop, "all")
#' head(bigout)
mst.seq = function(start, neighbors, cases, pop, w, ex, ty,
                   max_pop, type = "maxonly", nlinks = "one",
                   early = FALSE) {
  loglikrat = yin = ein = popin = numeric(length(neighbors))
  region = start
  yin[1] = cases[region]
  ein[1] = ex[region]
  popin[1] = pop[region]
  loglikrat[1] = scan.stat(yin[1], ein[1], ty - ein[1], ty)
  uz = max_neighbors = vector("list", length(neighbors))
  uz[[1]] = region
  max_neighbors[[1]] = setdiff(neighbors, region)

  stop = FALSE
  # double check that pop constraint is satisfied
  # if not, stop and set loglikrat to 0 (so that it's not
  # the maximum in the search)
  if (popin[1] > max_pop) {
    stop = TRUE
    loglikrat[1] = 0
  }

  i = 1
  while (!stop & i < length(neighbors)) {
    # get current zone to extend
    c_zone = uz[[i]]

    # current max set of neighbors for current zone
    cmn = max_neighbors[[i]]

    # which regions are connected to c_zone
    # and satisfy the population constraints
    sw = w[c_zone, cmn, drop = FALSE]

    # count number of connections between c_zone and
    # neighbors
    if (i == 1 | nlinks == "one") {
      # can only be single connected region first iteration
      connected = cmn[colSums(sw) >= 1]
    } else if (nlinks == "two") {
      connected = cmn[colSums(sw) >= 2]
    } else if (nlinks == "max") {
      nconnect = colSums(sw)
      connected = cmn[nconnect == max(nconnect)]
    }

    # new popin when adding each potential neighbor to c_zone
    p_popin = popin[i] + pop[connected]

    # candidate regions that satisfy constraints
    in_size = which(p_popin <= max_pop)
    cand_regions = connected[in_size]

    if (length(cand_regions) > 0) {
      # yin and ein for candidate zones (cur zone plus cand_regions)
      yin_cand = yin[i] + cases[cand_regions]
      ein_cand = ein[i] + ex[cand_regions]

      # test statistics for candidate locations
      stat_cand = scan.stat(yin = yin_cand, ein = ein_cand,
                            eout = ty - ein_cand, ty = ty)

      # index of max stat_cand
      max_idx = which.max(stat_cand)

      # only update if new candidate statistics
      # are greater than old statistics
      if (early & (stat_cand[max_idx] <= loglikrat[i])) {
        stop = TRUE
      } else {
        # update for next best zone
        loglikrat[i + 1] = stat_cand[max_idx]
        yin[i + 1] = yin_cand[max_idx]
        ein[i + 1] = ein_cand[max_idx]
        popin[i + 1] = p_popin[in_size][max_idx]
        max_neighbors[[i + 1]] = setdiff(max_neighbors[[i]],
                                         cand_regions[max_idx])
        uz[[i + 1]] = c(c_zone, cand_regions[max_idx])
        i = i + 1 # update counter
      }
    } else {
      # end algorithm
      stop = TRUE
    }
  }
  # return desired output
  if (type == "maxonly") {
    return(max(loglikrat))
  } else if (type == "pruned") {
    which_max = which.max(loglikrat)

    return(list(locids = uz[[which_max]],
                loglikrat = loglikrat[which_max],
                cases = yin[which_max],
                expected = ein[which_max],
                population = popin[which_max]))
  } else if (type == "all") {
    # return only non-zero elements
    return(list(locids = uz[[i]],
                loglikrat = loglikrat[1:i],
                cases = yin[1:i],
                expected = ein[1:i],
                population = popin[1:i]))
  }
}
