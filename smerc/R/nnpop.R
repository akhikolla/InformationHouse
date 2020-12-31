#' Determine nearest neighbors with population constraint
#'
#' \code{scan.nn} determines the nearest
#' neighbors for a set of observations based on the
#' distance matrix according to a population-based
#' upperbound.
#'
#' This function determines the nearest neighbors of each
#' centroid based on the intercentroid distance.  The number
#' of nearest neighbors is limited by the sum of the
#' population values among the nearest neighbors.  The set
#' of nearest neighbors can contain no more than \code{ubpop
#' * sum(pop)} members of the population.  The nearest
#' neighbors are ordered from nearest to farthest.
#'
#' @param d An \eqn{n\times n} square distance matrix
#'   containing the intercentroid distance between the
#'   \eqn{n} region centroids.
#' @inheritParams scan.test
#'
#' @return Returns the indices of the nearest neighbors as a
#'   list.  For each element of the list, the indices are
#'   ordered from nearest to farthest from each centroid.
#' @author Joshua French
#' @export
#' @aliases nnpop
#' @rdname scan.nn
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' d = as.matrix(dist(coords))
#' nn = scan.nn(d, pop = nydf$pop, ubpop = 0.1)
nnpop = function(d, pop, ubpop) {
  tpop = sum(pop) # total population
  # order distances for each region
  # results has each column showing order of indices
  # from shortest to largest distance
  od = apply(d, 1, order)

  # for each row of ordered distance matrix
  # sum the cumulative population size for
  # the expanding collection of regions

  # return the largest collection of regions for which the
  # cumulative population is less than the desired
  # proportion of the total popuation
  return(apply(od, 2,
               FUN = function(x) {
                 csum = cumsum(pop[x])
                 x[which(csum <= tpop * ubpop)]
               }))
}

#' @export
#' @rdname scan.nn
scan.nn = nnpop
