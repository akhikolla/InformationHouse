#' Convert a Pattern to Distances & Ratios
#'
#' Creates `n_sim` null models by permutation of the original pattern and
#' calculates distances between all object of a pattern closer than `max_dist`
#' and determines the fractions of the perimeter of buffers with distance
#' `dist` inside the study area (needed for edge correction).
#'
#' Null models are created by randomly rotating and randomly placing all
#' objects within the study area without overlap. They are used for correcting
#' the biased pcf and constructing a pointwise critical envelope (cf. Nuske
#' et al. 2009).
#'
#' Measuring distances between objects and permutation of the pattern is done
#' with [GEOS](https://trac.osgeo.org/geos) and spatial data are converted to
#' GEOS geometries by [GDAL/OGR](http://www.gdal.org).
#'
#' @param area,pattern Data source name of study area and pattern
#'        (interpretation varies by driver - for some drivers, dsn is a file
#'        name, but may also be a folder, or contain the name and access
#'        credentials of a database)
#' @param max_dist Maximum distance measured in the pattern. Usually smaller
#'        than half the diameter of the study area.
#' @param n_sim Number of simulated patterns (randomizations) to be generated
#'        for computing the envelope and correcting the biased emperical pcf.
#'        Determines together with `n_rank` in [dists2pcf()] the alpha level of
#'        the envelope. If `alpha` and `n_rank` are fix, n_sim can be
#'        calculated by `(n_rank*2/alpha)-1` eg. `(5*2/0.05)-1 = 199`.
#' @param max_tries How often shall a relocation of an object be tried
#'        while randomizing the pattern.
#' @param save_patterns Shall the simulated patterns be saved as Shapefiles
#'        for debugging/later inspections. Might be a large number of files
#'        (4 * n_sim). Can be `NULL` (no export) or a character string
#'        providing a basename optionally including a valid/existing path.
#' @param verbose Provide progress information in the console. Roman numerals
#'        (x: 10, C: 100, D: 500, M: 1000) indicate the progress of the
#'        simulation and 'e' the emperical PCF.
#'
#' @return An object of class [dists].
#'
#' @references
#'  Nuske, R.S., Sprauer, S. and Saborowski J. (2009)
#'  Adapting the pair-correlation function for analysing the spatial
#'  distribution of canopy gaps.
#'  \emph{Forest Ecology and Management}, \bold{259}(1), 107-–116.
#'  \doi{10.1016/j.foreco.2009.09.050}
#'
#' @seealso [dists2pcf()], [plot.fv_pcf()]
#'
#' @examples
#' # it's advised against setting n_sim < 199
#' ds <- pat2dists(area=system.file("shapes/sim_area.shp", package="apcf"),
#'                 pattern=system.file("shapes/sim_pat_reg.shp", package="apcf"),
#'                 max_dist=25, n_sim=3, verbose=TRUE)
#'
#' @export
pat2dists <- function(area, pattern, max_dist, n_sim=199,
                      max_tries=100000, save_patterns=NULL, verbose=FALSE){

  if(missing(area) || missing(pattern))
    stop("area and pattern should specify a data source or filename")

  if(missing(max_dist) || !is.numeric(max_dist))
    stop("max_dist must be given and must be numeric")

  if(length(area) > 1 || length(pattern) > 1)
    warning("using only the first element of area and pattern, respectively")

  if(is.null(save_patterns)){
    save_basename <- ' '
    save_patterns <- FALSE
  } else {
    save_dir <- dirname(save_patterns)
    if(!dir.exists(save_dir) || file.access(save_dir, 2) != 0){
      stop(paste0('Can not write in ', save_dir))
    } else {
      test_file <- paste0(save_patterns, '1.shp')
      if(file.exists(test_file))
        warning(paste0(save_patterns, '* exists and will be overwritten'))
      save_basename <- save_patterns
      save_patterns <- TRUE
    }
  }

  if(file.exists(area))
    area <- normalizePath(area)

  if(file.exists(pattern))
    pattern <- normalizePath(pattern)

  rand_dists_ratios(pattern[1], area[1], max_dist, as.integer(n_sim),
                    as.integer(max_tries), save_patterns, save_basename, verbose)
}
