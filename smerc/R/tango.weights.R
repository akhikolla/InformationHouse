#' Distance-based weights for \code{tango.test}
#'
#' \code{tango.weights} constructs a distance-based weights
#' matrix. The \code{tango.weights} function can be used to
#' construct a weights matrix \code{w} using the method of
#' Tango (1995), Rogerson (1999), or a basic style.
#'
#' \code{coords} is used to construct an \eqn{n \times n}
#' distance matrix \code{d}.
#'
#' If \code{type = "basic"}, then \eqn{w_{ij} =
#' exp(-d_{ij}/\kappa)}.
#'
#' If \code{type = "rogerson"}, then \eqn{w_{ij} =
#' exp(-d_{ij}/\kappa)/\sqrt(pop_i/pop * pop_j/pop)}.
#'
#' If \code{type = "tango"}, then \eqn{w_{ij} = exp(-4 *
#' d_{ij}^2/\kappa^2)}.
#'
#' @inheritParams scan.test
#' @param kappa A positive constant related to strength of
#'   spatial autocorrelation.
#' @param type The type of weights matrix to construct.
#'   Current options are \code{"basic"}, \code{"tango"}, and
#'   \code{"rogerson"}.  Default is \code{"basic"}.  See
#'   Details.
#'
#' @return Returns an \eqn{n \times n} matrix of weights.
#' @references Tango, T.  (1995) A class of tests for
#' detecting "general" and "focused" clustering of rare
#' diseases.  Statistics in Medicine.  14:2323-2334.
#'
#' Rogerson, P. (1999) The Detection of Clusters Using A
#' Spatial Version of the Chi-Square Goodness-of-fit Test.
#' Geographical Analysis. 31:130-147
#' @author Joshua French
#' @seealso \code{\link{tango.test}}
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' w = tango.weights(coords, kappa = 1, longlat = TRUE)
tango.weights = function(coords, kappa = 1, longlat = FALSE,
                         type = "basic", pop = NULL) {
  type = match.arg(type, c("basic", "rogerson", "tango"))
  arg_check_dweights(coords, kappa, longlat, type, pop)
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  if (type == "basic") {
    w = exp(-d / kappa)
  } else if (type == "rogerson") {
    popp = pop / sum(pop)
    w = exp(-d / kappa) / sqrt(outer(popp, popp))
  } else if (type == "tango") {
    w = exp(-4 * (d / kappa) ^ 2)
  }
  return(w)
}

#' @rdname tango.weights
#' @export
dweights = tango.weights

#' Argument checking for dweights/tango.weights
#'
#' @param coords Matrix of coordinates
#' @param kappa Scale parameter for weight function
#' @param longlat Logical value indicating whether distance
#' should be Euclidean (FALSE) or great circle (TRUE)
#' @param type Type of weights (basic, rogerson, or tango)
#' @param pop vector of population values
#'
#' @return NULL
#' @noRd
arg_check_dweights = function(coords, kappa, longlat, type, pop) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_dweights_kappa(kappa)
  arg_check_longlat(longlat)
  if (!is.null(pop)) {
    arg_check_pop(pop, N)
  }
  arg_check_dweights_type(type)
  if (type == "rogerson" & is.null(pop)) {
    stop("pop must be provided when type = 'rogerson'")
  }
}

