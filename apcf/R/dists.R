#' Class dists: Collection of Distances and Ratios
#'
#' Advanced Use Only. This low-level function creates an object of class
#' "dists" from raw numerical data.
#'
#' @param df A data frame with the columns `sim`, `dist`, and `ratio`
#'        containing an indicator of the model run (0:n_sim), distances between
#'        the objects of the patterns, and the ratios of a buffer with distance
#'        `dist` inside the study area (needed for Ripley's edge correction).
#' @param area size of the study area in square units
#' @param n_obj number of objects in the pattern
#' @param max_dist maximum distance to be measured in pattern
#' @param obj an R object, preferably of class `dists`
#'
#' @return An object of class `dists`.
#'
#' @export
dists <- function(df, area, n_obj, max_dist){
  stopifnot(is.data.frame(df))
  stopifnot(is.integer(n_obj))
  stopifnot(is.numeric(max_dist))

  if(area <= 0)
    stop("area must be >= 0")
  if(n_obj < 2)
    stop("n_obj must be >> 1")

  if(length(df) != 3 || !all(c("sim", "dist", "ratio") %in% names(df)))
    stop("'df' must have at least columns named sim, dist, ratio")

  attr(df, "area")     <- area
  attr(df, "n_obj")    <- n_obj
  attr(df, "max_dist") <- max_dist
  class(df) <- c("dists", class(df))

  return(df)
}


#' @rdname dists
#' @export
is.dists <- function(obj){
  inherits(obj, "dists")
}
