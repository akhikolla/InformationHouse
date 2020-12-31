#' Determine zones for the maxima likelihood
#' first algorithm.
#'
#' \code{mlf.zones} determines the most likely cluster zone
#' obtained by implementing the maxima likelihood first
#' scann method of Yao et al. (2011).  Note that this is
#' really just a special case of the dynamic minimum
#' spanning tree (DMST) algorithm of Assuncao et al. (2006)
#'
#' Each step of the mlf scan test seeks to maximize the
#' likelihood ratio test statistic used in the original
#' spatial scan test (Kulldorff 1997).  The first zone
#' considered is the region that maximizes this likelihood
#' ration test statistic, providing that no more than
#' \code{ubpop} proportion of the total population is in the
#' zone.  The second zone is the first zone and the
#' connected region that maximizes the scan statistic,
#' subject to the population and distance constraints.  This
#' pattern continues until no additional zones can be added
#' due to population or distance constraints.
#'
#' Every zone considered must have a total population less
#' than \code{ubpop * sum(pop)} in the study area.
#' Additionally, the maximum intercentroid distance for the
#' regions within a zone must be no more than \code{ubd *
#' the maximum intercentroid distance across all regions}.
#'
#' @inheritParams mlf.test
#' @inheritParams dmst.zones
#' @inherit dc.zones return
#' @author Joshua French
#' @references Yao, Z., Tang, J., & Zhan, F. B. (2011).
#'   Detection of arbitrarily-shaped clusters using a
#'   neighbor-expanding approach: A case study on murine
#'   typhus in South Texas. International Journal of Health
#'   Geographics, 10(1), 1.
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("x", "y")])
#' mlf.zones(coords, cases = floor(nydf$cases),
#'           pop = nydf$pop, w = nyw, longlat = TRUE)
mlf.zones = function(coords, cases, pop, w,
                     ex = sum(cases) / sum(pop) * pop,
                     ubpop = 0.5, ubd = 1, longlat = FALSE) {
  # sanity checking
  arg_check_dmst_zones(coords = coords, cases = cases,
                       pop = pop, w = w, ex = ex,
                       ubpop = ubpop, ubd = ubd,
                       longlat = longlat, type = "all",
                       progress = FALSE)
  # total number of cases
  ty = sum(cases)

  # calculate test statistics for each individual region
  # yin, ein, eout, ty
  tobs = scan.stat(cases, ex, ty - ex, ty)

  # determine starting region for maxima likelihood first algorithm
  start = which.max(tobs)

  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = longlat)

  # upperbound for population in zone
  max_pop = ubpop * sum(pop)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)

  # find neighbors of starting region
  start_neighbors = which(d[start, ] <= max_dist)

  # return sequence of candidate zones (or a subset depending on type)
  mst.seq(start, start_neighbors, cases, pop, w, ex, ty,
          max_pop, type = "pruned")
}
