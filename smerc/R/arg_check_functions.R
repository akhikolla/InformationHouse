#' Check coords argument
#'
#' @param coords An matrix-like object with N rows
#' @noRd
arg_check_coords = function(coords) {
  if (!(is.matrix(coords) | is.data.frame(coords))) {
    stop("coords should be a matrix or a data frame")
  }
  if (ncol(coords) != 2) {
    stop("coords must have two columns")
  }
}

#' Check cases argumnent
#'
#' @param cases A numeric vector of cases
#' @param N The number of rows in coords
#' @noRd
arg_check_cases = function(cases, N) {
  if (length(cases) != N) {
    stop("length(cases) != nrow(coords)")
  }
  if (!is.numeric(cases)) {
    stop("cases should be a numeric values")
  }
  if (!is.vector(cases)) {
    stop("cases should be a vector")
  }
  if (min(cases) < 0) {
    stop("cases must have non-negative values")
  }
}

#' Check population argument
#'
#' @param pop A vector of population values
#' @param N The number of rows in coords
#' @noRd
arg_check_pop = function(pop, N) {
  if (length(pop) != N) {
    stop("length(pop) != nrow(coords)")
  }
  if (!is.numeric(pop)) {
    stop("pop should be numeric values")
  }
  if (!is.vector(pop)) {
    stop("pop should be a vector")
  }
  if (min(pop) < 1) {
    stop("pop values must be >= 1")
  }
}

#' Check cstar argument
#'
#' @param cstar Case radius
#' @param cases Numeric vector of cases
#' @noRd
arg_check_cstar = function(cstar, cases) {
  if (length(cstar) != 1 || !is.numeric(cstar)) {
    stop("cstar should be a numeric vector of length 1")
  }
  if (cstar < 1 || cstar > sum(cases)) {
    stop("cstar should be at least 1 and less than or equal to the sum(cases)")
  }
}

#' Check longlat argument
#'
#' @param longlat A logical value indicating whether longlat
#' distance should be used (TRUE)
#' @noRd
arg_check_longlat = function(longlat) {
  if (length(longlat) != 1) {
    stop("length(longlat) != 1")
  }
  if (!is.logical(longlat)) {
    stop("longlat should be a logical value")
  }
}

#' Check alpha argument
#'
#' @param alpha Signifcance level (single value > 0 and <= 1)
#' @noRd
arg_check_alpha = function(alpha) {
  if (length(alpha) != 1 || !is.numeric(alpha)) {
    stop("alpha should be a numeric vector of length 1")
  }
  if (alpha < 0 || alpha > 1) {
    stop("alpha should be a value [0, 1)")
  }
}

#' Check noc argument
#'
#' @param noc Logical value. Should only non-overlapping clusters be returned.
#' @noRd
arg_check_noc = function(noc) {
  if (length(noc) != 1) {
    stop("length(noc) != 1")
  }
  if (!is.logical(noc)) {
    stop("noc should be a logical value")
  }
}

#' Check ex argument
#'
#' @param ex A vector of expected counts
#' @param N nrow(coords)
#' @noRd
arg_check_ex = function(ex, N) {
  if (length(ex) != N) {
    stop("length(ex) != nrow(coords)")
  }
  if (!is.numeric(ex)) {
    stop("ex should be numeric values")
  }
  if (!is.vector(ex)) {
    stop("ex should be a vector")
  }
}

#' Check modified argument
#'
#' @param modified A logicical value. should bn.test be modified (TRUE)
#' @return NULL
#' @noRd
arg_check_modified = function(modified) {
  if (length(modified) != 1) {
    stop("length(modified) != 1")
  }
  if (!is.logical(modified)) {
    stop("modified should be a logical value")
  }
}

#' Check tobs argument
#'
#' For smerc_cluster function
#'
#' @param tobs Vector of observed test statistics
#' @return NULL
#' @noRd
arg_check_tobs = function(tobs) {
  if (!is.numeric(tobs)) {
    stop("tobs should be numeric values")
  }
  if (!is.vector(tobs)) {
    stop("tobs should be a vector")
  }
  if (min(tobs) < 0) {
    stop("tobs values must be >= 0")
  }
}

#' Check zones argument
#'
#' For smerc_cluster function
#'
#' @param zones A list of zones
#' @param N Number of tobs
#' @return NULL
#' @noRd
arg_check_zones = function(zones, N) {
  if (length(zones) != N) {
    stop("length(zones) != length(tobs)")
  }
  if (!is.list(zones)) {
    stop("zones should be a list")
  }
}

#' Check pvalue argument
#'
#' For smerc_cluster function
#'
#' @param pvalue Vector of p-values
#' @param N length(tobs)
#' @return NULL
#' @noRd
arg_check_pvalue = function(pvalue, N) {
  if (length(pvalue) != N) {
    stop("length(pvalue) != nrow(coords)")
  }
  if (!is.numeric(pvalue)) {
    stop("pvalue should be numeric values")
  }
  if (!is.vector(pvalue)) {
    stop("pvalue should be a vector")
  }
  if (min(pvalue) < 0 | max(pvalue) > 1) {
    stop("pvalue must have values in [0, 1]")
  }
}

arg_check_d = function(d, N) {
  if (!is.matrix(d)) {
    stop("d must be a matrix")
  }
  if (nrow(d) != N | ncol(d) != N) {
    stop("d must be square nrow(d) = nrow(coords)")
  }
}

#' Check method argument
#'
#' @param method A single character vector specifying
#' the method name.
#' @return NULL
#' @noRd
arg_check_method = function(method) {
  if (length(method) != 1) {
    stop("method must be a vector of length 1")
  }
  if (!is.vector(method)) {
    stop("method must be a vector of length 1")
  }
  if (!is.character(method)) {
    stop("method must be a character vector")
  }
}

#' Check method argument
#'
#' @param rel_param A list with relevant paramameters for a method
#' the method name.
#' @return NULL
#' @noRd
arg_check_rel_param = function(rel_param) {
  if (!is.list(rel_param)) {
    stop("rel_param must be a list")
  }
}

#' Check shape_all argument
#'
#' For smerc_cluster function
#'
#' @param shape_all Vector of shapes
#' @param N length(tobs)
#' @return NULL
#' @noRd
arg_check_shape_all = function(shape_all, N) {
  if (length(shape_all) != N) {
    stop("length(shapes_all) != length(tobs)")
  }
  if (!is.numeric(shape_all)) {
    stop("shape_all should be numeric values")
  }
  if (!is.vector(shape_all)) {
    stop("shape_all should be a vector")
  }
  if (min(shape_all) < 1) {
    stop("All shapes must be >= 1")
  }
}

#' Check angle_all argument
#'
#' For smerc_cluster function
#'
#' @param angle_all Vector of shapes
#' @param N length(tobs)
#' @return NULL
#' @noRd
arg_check_angle_all = function(angle_all, N) {
  if (length(angle_all) != N) {
    stop("length(shapes_all) != length(tobs)")
  }
  if (!is.numeric(angle_all)) {
    stop("angle_all should be numeric values")
  }
  if (!is.vector(angle_all)) {
    stop("angle_all should be a vector")
  }
  if (min(angle_all) < 0 | max(angle_all) >= 360) {
    stop("All angles must be in [0, 360)")
  }
}

#' Check a argument
#'
#' @param a Penalty parameter for elliptic.test
#' @return NULL
#' @noRd
arg_check_a = function(a) {
  if (length(a) != 1) {
    stop("a must be a single value")
  }
  if (!is.numeric(a)) {
    stop("a must be a numeric value")
  }
  if (!is.vector(a)) {
    stop("a must be a vector")
  }
  if (a < 0) {
    stop("a must be >= 0")
  }
}

#' Check w argument
#'
#' @param w Spatial connectivity matrix
#' @param N nrow(coords)
#' @return NULL
#' @noRd
arg_check_w = function(w, N) {
  if (!is.matrix(w) & !is.data.frame(w)) {
    stop("w must be a matrix or data.frame")
  }
  if (nrow(w) != N | ncol(w) != N) {
    stop("w must be a square matrix with nrow(w) = nrow(coords)")
  }
  if (any(w != 0 & w != 1)) {
    stop("w must be 0s and 1s")
  }
}

#' Check nsim argument
#'
#' @param nsim A non-negative integer
#'
#' @return NULL
#' @noRd
arg_check_nsim = function(nsim) {
  if (length(nsim) != 1) {
    stop("nsim must be a single value")
  }
  if (!is.numeric(nsim)) {
    stop("nsim must be a numeric value")
  }
  if (!is.vector(nsim)) {
    stop("nsim must be a vector (of length 1)")
  }
  if (min(nsim) < 0) {
    stop("nsim must be a non-negative integer")
  }
}

#' Check ubpop argument
#'
#' @param ubpop A positive value
#' @return NULL
#' @noRd
arg_check_ubpop = function(ubpop) {
  if (length(ubpop) != 1) {
    stop("ubpop must be a single value")
  }
  if (!is.numeric(ubpop)) {
    stop("ubpop must be a numeric value")
  }
  if (!is.vector(ubpop)) {
    stop("ubpop must be a vector (of length 1)")
  }
  if (ubpop <= 0 || ubpop > 1) {
    stop("ubpop should be a value between 0 and 1")
  }
}

arg_check_k = function(k, N) {
  if (length(k) != 1) {
    stop("k must have length 1")
  }
  if (!is.numeric(k)) {
    stop("k must be a numeric value")
  }
  if (!is.vector(k)) {
    stop("k must be a vector (of length 1)")
  }
  if (k < 1) {
    stop("k must be an integer >= 1")
  }
  if (floor(k) > N) {
    stop("k cannot be more than the number of regions")
  }
}

#' Check simdist argument
#'
#' @param simdist Distribution of simulation, single character value
#' @return NULL
#' @noRd
arg_check_simdist = function(simdist) {
  if (!is.null(simdist)) {
    if (length(simdist) != 1) {
      stop("simdist must be of length 1")
    }
    if (!is.element(simdist, c("multinomial", "poisson", "binomial"))) {
      stop("simdist must be 'multinomial', 'poisson', or 'binomial'")
    }
  }
}

arg_check_min_cases = function(min.cases) {
  if (length(min.cases) != 1) {
    stop("min.cases must be a single value")
  }
  if (!is.numeric(min.cases)) {
    stop("min.cases must be a numeric value")
  }
  if (!is.vector(min.cases)) {
    stop("min.cases must be a vector (of length 1)")
  }
  if (min.cases < 1) {
    stop("min.cases must be be >= 1")
  }
}

arg_check_type = function(type) {
  if (length(type) != 1) {
    stop("type must be a single value")
  }
  if (!is.character(type)) {
    stop("type must be a character")
  }
  if (!is.element(type, c("poisson", "binomial"))) {
    stop("type must be 'poisson' or 'binomial'")
  }
}

#' Check shape argument
#'
#' @param shape A vector of shape values >= 1
#' @return NULL
#' @noRd
arg_check_shape = function(shape) {
  if (!is.numeric(shape)) {
    stop("shape must be a numeric vector")
  }
  if (!is.vector(shape)) {
    stop("shape must be a numeric vector")
  }
  if (min(shape) < 1) {
    stop("shape must be >= 1")
  }
}

#' Check nangle argument
#'
#' @param nangle A vector of nangle values >= 1
#' @return NULL
#' @noRd
arg_check_nangle = function(nangle) {
  if (!is.numeric(nangle)) {
    stop("nangle must be a numeric vector")
  }
  if (!is.vector(nangle)) {
    stop("nangle must be a numeric vector")
  }
  if (min(nangle) < 1) {
    stop("nangle must be >= 1")
  }
}

#' Title
#'
#' @param nstar A numeric value indicating the window radius
#' for cepp.test.
#' @return NULL
#' @noRd
arg_check_nstar = function(nstar, pop) {
  if (length(nstar) != 1) {
    stop("nstar should be a single value")
  }
  if (!is.numeric(nstar)) {
    stop("nstar must be a numeric value")
  }
  if (!is.vector(nstar)) {
    stop("nstar must be a vector (of length 1)")
  }
  if (nstar < 1) {
    stop("nstar should be at least 1")
  }
  if (nstar > sum(pop)) {
    stop("nstar should be no more than sum(pop)")
  }
}

#' Check tango.weights/dweights type argument
#'
#' @param type A character vector: basic, rogerson, tango
#' @return NULL
#' @noRd
arg_check_dweights_type = function(type) {
  if (length(type) != 1) {
    stop("type must be a single name")
  }
  if (!is.character(type)) {
    stop("type must be a character")
  }
  if (!is.vector(type)) {
    stop("type must be a vector")
  }
  if (!is.element(type, c("basic", "rogerson", "tango"))) {
    stop("invalid type")
  }
}

#' Check dweights kappa argument
#'
#' @param kappa A positive value
#' @return NULL
#' @noRd
arg_check_dweights_kappa = function(kappa) {
  if (length(kappa) != 1) {
    stop("kappa should be a single value")
  }
  if (!is.numeric(kappa)) {
    stop("kappa must be a numeric value")
  }
  if (!is.vector(kappa)) {
    stop("kappa must be a vector (of length 1)")
  }
  if (kappa <= 0) {
    stop("kappa must be positive")
  }
}

#' Check w argument of tango.test
#'
#' @param w Spatial weights matrix
#' @param N nrow(coords)
#' @return NULL
#' @noRd
arg_check_tango_w = function(w, N) {
  if (!is.matrix(w) & !is.data.frame(w)) {
    stop("w must be a matrix or data.frame")
  }
  if (nrow(w) != N | ncol(w) != N) {
    stop("w must be a square matrix with nrow(w) = nrow(coords)")
  }
  if (!is.numeric(w)) {
    stop("w must be numeric")
  }
}

#' Check check.unique argument (of uls.zones)
#'
#' @param check.unique A single logical value
#' @return NULL
#' @noRd
arg_check_check_unique = function(check.unique) {
  if (length(check.unique) != 1) {
    stop("check.unique must be a single value")
  }
  if (!is.logical(check.unique)) {
    stop("check.unique must be a logical value")
  }
  if (!is.vector(check.unique)) {
    stop("check.unique must be a vector (of length 1)")
  }
}

#' Check total population argument
#'
#' @param tpop Total population
#' @return NULL
#' @noRd
arg_check_tpop = function(tpop) {
  if (length(tpop) != 1) {
    stop("tpop must be a single value")
  }
  if (!is.numeric(tpop)) {
    stop("tpop must be numeric")
  }
  if (!is.vector(tpop)) {
    stop("tpop must be a vector (of length 1)")
  }
  if (tpop <= 0) {
    stop("tpop must be >= 1")
  }
}

#' Check total cases argument
#'
#' @param ty Total cases
#' @return NULL
#' @noRd
arg_check_ty = function(ty) {
  if (length(ty) != 1) {
    stop("ty must be a single value")
  }
  if (!is.numeric(ty)) {
    stop("ty must be numeric")
  }
  if (!is.vector(ty)) {
    stop("ty must be a vector (of length 1)")
  }
  if (ty <= 0) {
    stop("ty must be >= 1")
  }
}

#' Check .sim arguments for type = "poisson"
#'
#' @param ein A vector of expected cases in each zone
#' @param eout A vector of expected cases outside of each zone
#' @param nz The number of zones
#' @return NULL
#' @noRd
arg_check_sim_poisson_type = function(ein, eout, nz) {
  if (is.null(ein) | is.null(eout)) {
    stop("ein and eout must be provided when type='poisson'")
  }
  if (nz != length(ein)) {
    stop("ein has improper length")
  }
  if (!is.vector(ein)) {
    stop("ein must be a vector")
  }
  if (!is.numeric(ein)) {
    stop("ein must be numeric")
  }
  if (nz != length(eout)) {
    stop("eout has improper length")
  }
  if (!is.vector(eout)) {
    stop("eout must be a vector")
  }
  if (!is.numeric(eout)) {
    stop("eout must be numeric")
  }
}

arg_check_sim_binomial_type = function(popin, popout, tpop, nz) {
  if (is.null(popin) | is.null(popout) | is.null(tpop)) {
    stop("popin, popout, and tpop must be provided when type='binomial'")
  }
  if (nz != length(popin)) {
    stop("popin has improper length")
  }
  if (!is.vector(popin)) {
    stop("popin must be a vector")
  }
  if (!is.numeric(popin)) {
    stop("popin must be numeric")
  }
  if (nz != length(popout)) {
    stop("popout has improper length")
  }
  if (!is.vector(popout)) {
    stop("popout must be a vector")
  }
  if (!is.numeric(popout)) {
    stop("popout must be numeric")
  }
  arg_check_tpop(tpop)
}

#' Check ubd argument
#'
#' @param ubd Distance upperbound (in terms of proportion)
#' @return NULL
#' @noRd
arg_check_ubd = function(ubd) {
  if (length(ubd) != 1) {
    stop("ubd must be a single number")
  }
  if (!is.numeric(ubd)) {
    stop("ubd must be numeric")
  }
  if (!is.vector(ubd)) {
    stop("ubd must be a vector (of length 1)")
  }
  if (ubd <= 0 | ubd > 1) {
    stop("ubd must be in (0, 1]")
  }
}

#' Check nclusters argument
#'
#' @param nclusters Number of clusters to plot
#' @param N length(x$clusters) from a smerc_cluster
#'
#' @return NULL
#' @noRd
arg_check_nclusters = function(nclusters, N) {
  if (length(nclusters) != 1) {
    stop("nclusters must have length 1")
  }
  if (!is.numeric(nclusters)) {
    stop("nclusters must be a numeric value")
  }
  if (!is.vector(nclusters)) {
    stop("nclusters must be a vector (of length 1)")
  }
  if (nclusters < 1) {
    stop("nclusters must be >= 1")
  }
  if (nclusters > N) {
    stop("nclusters must be <= length(x$clusters)")
  }
}
