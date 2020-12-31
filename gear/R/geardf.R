#' Construct a \code{geardf}
#'
#' \code{geardf} constructs a \code{geardf}, which is simply
#' a \code{data.frame} with an attribute indicating which
#' columns refer to the coordinates. The function either
#' combines \code{x} and \code{coords} if \code{coordnames}
#' isn't provided, or identifies the columns of \code{x} to
#' which \code{coordnames} refers. The \code{geardf} class
#' is only provided to make plotting results of the
#' \code{predict.geolm*} functions simple without depending
#' on the \code{sp} or \code{sf} packages. See
#' \code{\link[gear]{plot.geardf}} for easy plotting.
#'
#' @param data A data.frame
#' @param coords A data.frame or matrix with column names
#' @inheritParams evgram
#' @return A \code{geardf}
#' @seealso \code{\link[gear]{plot.geardf}}
#' @export
#' @examples
#' dtf = data.frame(a = 1:2, b = 3:4)
#'
#' # create geardf with matrix coords (note column names)
#' coords = matrix(rnorm(4), ncol = 2)
#' colnames(coords) = c("u", "v")
#' geardf(dtf, coords)
#'
#' # create geardf with data.frame coords
#' coords = as.data.frame(coords)
#' geardf(dtf, coords)
#'
#' # create geardf using coordnames
#' dtf2 = cbind(dtf, u = rnorm(2), v = rnorm(2))
#' # vector form of coordnames
#' geardf(dtf2, coordnames = c("u", "v"))
#' # formula form of coordnames
#' geardf(dtf2, coordnames = ~ u + v)
#' # column index forum of coordnames
#' geardf(dtf2, coordnames = 3:4)
#'
#' gdf = geardf(dtf2, coordnames = 3:4)
#' # looks like a data.frame
#' gdf
#' # but slightly more complicated
#' class(gdf)
#' attr(gdf, "coordnames")
geardf = function(data, coords, coordnames) {
  if (missing(coords)) coords = NULL
  if (missing(coordnames)) coordnames = NULL
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  if (!is.null(coords)) {
    coordnames = arg_check_geardf(data, coords, coordnames)
    data = cbind(data, coords)
  } else {
    coordnames = arg_check_coordnames(coordnames, names(data))
  }
  class(data) = c("geardf", "data.frame")
  attr(data, "coordnames") = coordnames
  return(data)
}

#' Check arguments of geardf
#'
#' @param data A data.frame
#' @param coords NULL or a matrix or data.frame with 2 columns
#' @param coordnames A vector of coordinate names
#' @return A vector of coordinate names
#' @noRd
arg_check_geardf = function(data, coords, coordnames) {
  if (!is.null(coords) & !is.data.frame(coords) & !is.matrix(coords)) {
    stop("coords must be a data.frame or matrix or NULL")
  }
  if (is.null(coords) & is.null(coordnames)) {
    stop("either coords and coordnames must be provided")
  }
  if (!is.null(coords) & !is.null(coordnames)) {
    stop("only one of coords, coordnames can be provided")
  }
  nr = ifelse(is.null(coords), nrow(data), nrow(coords))
  if (nrow(data) != nr) {
    stop("nrow(data) != nrow(coords)")
  }
  if (is.data.frame(coords) | is.matrix(coords)) {
    if (ncol(coords) != 2) {
      stop("coords can only have two columns")
    }
  }
  if (is.matrix(coords)) {
    if (is.null(colnames(coords))) {
      stop("coords doesn't appear to have column names")
    }
    coordnames = colnames(coords)
  }
  if (is.data.frame(coords)) {
    coordnames = names(coords)
  }
  return(coordnames)
}
