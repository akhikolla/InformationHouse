#' Determine sequence of ULS zones.
#'
#' \code{uls.zones} determines the unique zones obtained by
#' implementing the ULS (Upper Level Set) test of Patil and
#' Taillie (2004).
#'
#' The zones returned must have a total population less than
#' \code{ubpop * sum(pop)} of all regions in the study area.
#'
#' @inheritParams uls.test
#' @return Returns a list of zones to consider for
#'   clustering.  Each element of the list contains a vector
#'   with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Patil, G.P. & Taillie, C. Upper level set
#'   scan statistic for detecting arbitrarily shaped
#'   hotspots. Environmental and Ecological Statistics
#'   (2004) 11(2):183-197.
#'   <doi:10.1023/B:EEST.0000027208.48919.7e>
#' @examples
#' data(nydf)
#' data(nyw)
#' uls.zones(cases = nydf$cases, pop = nydf$population, w = nyw)
uls.zones = function(cases, pop, w, ubpop = 0.5,
                     check.unique = FALSE) {
  arg_check_uls_zones(cases = cases, pop = pop, w = w,
                      ubpop = ubpop,
                      check.unique = check.unique)

  # order rates from largest to smallest
  or = order(cases / pop, decreasing = TRUE);
  # reorder rows and columns by order of r
  w = w[or, ]
  w = w[, or]

  #unique zones, first zone always the starting point
  uz = vector("list", nrow(w))
  uz[[1]] = 1

  # current zones
  cz = 1
  for (i in 2:nrow(w)) {
    #are regions adjacent
    which_adjacent = which(w[1:(i - 1), i] == 1)
    # if there are no neighbors for the new vertex
    if (length(which_adjacent) == 0) {
      # add new zone to list
      uz[[i]] = i
      # update current zones
      cz = c(cz, i)
    } else {
      # which zones intersect
      wzi = which(unlist(lapply(uz[cz], function(x) {
        length(intersect(x, which_adjacent)) > 0
      }
      ), use.names = FALSE))
      # add new zone to list
      uz[[i]] = c(unlist(uz[cz[wzi]], use.names = FALSE), i)
      # update current zones
      cz = c(cz[-wzi], i)
    }
  }
  # convert back to original location ids
  uz = lapply(uz, function(x) or[x])
  # return only the zones that meet constraint for population upper bound
  popin = zones.sum(uz, pop)
  uz = uz[which(popin <= sum(pop) * ubpop)]

  if (check.unique) {
    # special case when rates are the same in some regions
    g = cases / pop
    ug = unique(g)
    if (length(g) > length(ug)) {
      uz2 = vector("list", length(g))
      og = g[or]
      counter = 0
      for (j in sort(ug, decreasing = TRUE)) {
        wl = which(og == j)
        if (length(wl) == 1) {
          counter = counter + 1
          uz2[counter] = uz[wl]
        } else {
          on = order(sapply(uz[wl], length), decreasing = TRUE)
          uzwlon = uz[wl[on]]
          uzwlon = uzwlon[noz(uzwlon)]
          for (k in seq_along(uzwlon)) {
            counter = counter + 1
            uz2[counter] = uzwlon[k]
          }
        }
      }
      uz = uz2[seq_len(counter)]
    }
  }
  return(uz)
}

#' Check uls.zones arguments
#'
#' @param cases A vector of cases
#' @param pop A vector of populations for each region
#' @param w A spatial adjacency matrix
#' @param ubpop A population upperbound
#' @param check.unique A logical value. TRUE means check for
#' unique rates between regions.
#' @return NULL
#' @noRd
arg_check_uls_zones = function(cases, pop, w, ubpop,
                               check.unique) {
  N = length(cases)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_w(w, N)
  arg_check_ubpop(ubpop)
  arg_check_check_unique(check.unique)
}
