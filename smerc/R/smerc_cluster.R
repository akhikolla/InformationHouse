#' Prepare \code{smerc_cluster}
#'
#' \code{smerc_cluster} prepares a \code{smerc_cluster}.
#'
#' @param tobs The vector of observed test statistics for each zone
#' @param zones A list of zones
#' @param pvalue The p-value associated with each test statistic
#' @inheritParams flex.test
#' @param d A precomputed distance matrix based on \code{coords}
#' @param method A character string indicating the method
#' used to construct the \code{smerc_cluster}.
#' @param rel_param A names list with the relevant parameters
#' associated with \code{method}.
#' @param a A single value >= 0 indicating the penalty to use
#' for \code{\link{elliptic.test}}.
#' @param shape_all A vector of shape parameters associated
#' with \code{zones}.
#' @param angle_all A vector of angle parameter associated with
#' \code{zones}.
#' @param alpha The significance level of the test.
#' @inheritParams scan.test
#' @return A \code{smerc_cluster} object. The object
#' generally has the following components:
#' \item{clusters}{A list containing information about the significant
#' clusters.  See further details below.}
#' \item{coords}{The matrix of centroid coordinates.}
#' \item{number_of_regions}{The number of regions considered.}
#' \item{total_population}{The total population in the regions.}
#' \item{total_cases}{The total number of cases in the regions.}
#' \item{cases_per_100k}{The rate of cases per 100,000 persons.}
#' \item{method}{The name of the method applied.}
#' \item{rel_param}{A list of relevant method parameters.}
#' \item{alpha}{The significance level.}
#' \item{longlat}{A logical value indicating which type
#' of distance was used.}
#'
#' Each element of the \code{clusters} component has:
#' \item{locids}{The ids of the regions in the cluster.}
#' \item{centroid}{The cluster centroid.}
#' \item{r}{The radius of the region (from the starting
#' region to last region of the cluster).}
#' \item{max_dist}{The maximum intercentroid distance between
#' all the regions in the cluster.}
#' \item{population}{The total population in the cluster.}
#' \item{cases}{The number of cases in the cluster.}
#' \item{expected}{The expected number of cases in the cluster.}
#' \item{smr}{Standardized mortality ratio
#' (\code{cases/expected}) in the cluster.}
#' \item{rr}{Relative risk in the cluster window. This is
#' \code{(cases/pop)/((total_cases - cases)/
#' (total_population - population))}.}
#' \item{loglikrat}{The log of the likelihood ratio test
#' statistic for the cluster. Only valid for the scan-type
#' tests.}
#' \item{test_statistic}{The test statistic for the cluster.}
#' \item{pvalue}{The p-value of the test statistic
#' associated with the cluster.}
#' \item{w}{The adjacency information for the cluster.}
#'
#' For \code{\link{elliptic.test}}, \code{clusters}
#' additionally has:
#' \item{semiminor_axis}{The semi-minor axis length for the
#' ellipse.}
#' \item{semimajor_axis}{The
#' semi-major axis length for the ellipse.}
#' \item{angle}{The rotation angle of the ellipse.}
#' \item{shape}{The shape of the ellipse.}
#' @export
smerc_cluster = function(tobs, zones, pvalue,
                         coords, cases, pop, ex, longlat,
                         method, rel_param,
                         alpha,
                         w = NULL, d = NULL,
                         a = NULL, shape_all = NULL,
                         angle_all = NULL) {
  arg_check_smerc_cluster(tobs = tobs, zones = zones,
                          pvalue = pvalue, coords = coords,
                          cases = cases, pop = pop, ex = ex,
                          longlat = longlat,
                          method = method,
                          rel_param = rel_param,
                          w = w, d = d, a = a,
                          shape_all = shape_all,
                          angle_all = angle_all,
                          alpha)
  new_smerc_cluster(tobs = tobs, zones = zones,
                    pvalue = pvalue, coords = coords,
                    cases = cases, pop = pop, ex = ex,
                    longlat = longlat,
                    method = method,
                    rel_param = rel_param, alpha = alpha,
                    w = w, d = d, a = a,
                    shape_all = shape_all,
                    angle_all = angle_all)
}

#' Construct \code{smerc_cluster}
#'
#' Doesn't check arguments, which is done in
#' \code{smerc_cluster}
#' @return A \code{smerc_cluster}
#' @noRd
new_smerc_cluster = function(tobs, zones, pvalue, coords,
                             cases, pop, ex, longlat,
                             method, rel_param, alpha, w, d,
                             a, shape_all, angle_all) {
  # total cases and population
  ty = sum(cases)
  tpop = sum(pop)

  # for the unique, non-overlapping clusters in order of significance,
  # find the associated centroid,
  # zone radius, cases in window, expected cases in window,
  # population in window, standarized mortality ratio, relative risk
  centroid_id = sapply(zones, utils::head, n = 1)
  boundary_id = sapply(zones, utils::tail, n = 1)

  if (!is.null(d)) {
    zone_r = d[cbind(centroid_id, boundary_id)]
  } else {
    zone_r = rep(NA, length(centroid_id))
  }
  # maximum centroid distance
  max_dist = unname(sapply(zones, function(x) {
    max(sp::spDists(coords[x,, drop = FALSE], longlat = longlat))
  }))
  # zone centroid
  # cases, expected, population in zone
  centroid = coords[centroid_id,, drop = FALSE]
  yin = zones.sum(zones, cases)
  ein = zones.sum(zones, ex)
  popin = zones.sum(zones, pop)
  smr = yin / ein
  rr_num = yin / popin
  rr_den = (ty - yin) / (tpop - popin)
  rr = rr_num / rr_den
  if (is.null(w)) {
    sig_w = lapply(zones, function(x) {
      matrix(c(0, rep(1, length(x) - 1)), nrow = 1)
    })
  } else {
    sig_w = sapply(zones, function(x) {
      w[x, x, drop = FALSE]
    }, simplify = FALSE)
  }
  if (!is.null(a)) {
    minor = unname(sapply(seq_along(zones), function(i) {
      first = zones[[i]][1]
      last = utils::tail(zones[[i]], 1)
      dist.ellipse(coords[c(first, last),, drop = FALSE],
                   shape = shape_all[i],
                   angle = angle_all[i])[1, 2]
    }))
    major = minor * shape_all
    loglikrat = stat.poisson(yin, ty - yin, ein, ty - ein)
  } else if (method == "Besag-Newell") {
    loglikrat = NA
  } else {
    loglikrat = tobs
  }

  # reformat output for return
  clusters = vector("list", length(zones))
  for (i in seq_along(clusters)) {
    clusters[[i]]$locids = zones[[i]]
    clusters[[i]]$centroid = centroid[i,, drop = FALSE]
    clusters[[i]]$r = zone_r[i]
    clusters[[i]]$max_dist = max_dist[i]
    if (!is.null(a)) {
      clusters[[i]]$semiminor_axis = minor[i]
      clusters[[i]]$semimajor_axis = major[i]
      clusters[[i]]$angle = angle_all[i]
      clusters[[i]]$shape = shape_all[i]
    }
    clusters[[i]]$population = popin[i]
    clusters[[i]]$cases = yin[i]
    clusters[[i]]$expected = ein[i]
    clusters[[i]]$smr = smr[i]
    clusters[[i]]$rr = rr[i]
    clusters[[i]]$loglikrat = loglikrat[i]
    clusters[[i]]$test_statistic = tobs[i]
    clusters[[i]]$pvalue = pvalue[i]
    clusters[[i]]$w = sig_w[[i]]
  }
  structure(list(clusters = clusters,
                 coords = coords,
                 number_of_regions = length(cases),
                 total_population = tpop,
                 total_cases = ty,
                 cases_per_100k = ty / tpop * 1e5,
                 method = method,
                 rel_param = rel_param,
                 alpha = alpha,
                 longlat = longlat),
            class = "smerc_cluster")
}

#' Check arguments for \code{smerc_cluster}
#'
#' \code{smerc_cluster} prepares and returns a
#' \code{smerc_cluster}.
#'
#' @param tobs The vector of observed test statistics for each zone
#' @param zones A list of zones
#' @param pvalue The p-value associated with each test statistic
#' @inheritParams flex.test
#' @param d A precomputed distance matrix based on \code{coords}
#' @param method A character string indicating the method
#' used to construct the \code{smerc_cluster}.
#' @param rel_param A names list with the relevant parameters
#' associated with \code{method}.
#' @param a A single value >= 0 indicating the penalty to use
#' for \code{\link{elliptic.test}}.
#' @param shape_all A vector of shape parameters associated
#' with \code{zones}. Relevant for \code{\link{elliptic.test}}.
#' @param angle_all A vector of angle parameter associated with
#' \code{zones}. Relevant for \code{\link{elliptic.test}}.
#' @param alpha The significance level.
#' @return A \code{smerc_cluster} object
#' @noRd
arg_check_smerc_cluster = function(tobs, zones, pvalue,
                                   coords, cases, pop, ex,
                                   longlat, method,
                                   rel_param, w, d, a,
                                   shape_all, angle_all,
                                   alpha) {
  arg_check_tobs(tobs)
  Ntobs = length(tobs)
  arg_check_zones(zones, Ntobs)
  arg_check_pvalue(pvalue, Ntobs)
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_ex(ex, N)
  arg_check_longlat(longlat)
  arg_check_method(method)
  arg_check_rel_param(rel_param)
  if (!is.null(w)) {
    arg_check_w(w, N)
  }
  if (!is.null(d)) {
    arg_check_d(d, N)
  }
  if (!is.null(a)) {
    arg_check_a(a)
  }
  if (!is.null(shape_all)) {
    arg_check_shape_all(shape_all, Ntobs)
  }
  if (!is.null(angle_all)) {
    arg_check_angle_all(angle_all, Ntobs)
  }
  arg_check_alpha(alpha)
}
