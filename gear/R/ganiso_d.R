#' Anisotropic distance-related characteristics
#'
#' Computes necessary distance-related characteristics when
#' there is geometric anisotropy. This is essentially an
#' internal function to \code{\link{evaluate}} a
#' \code{cmodStd} object produced by \code{\link{cmod_std}}
#' when the anisotropy ratio differs from 1.
#'
#' @param coords2 An \eqn{M \times 2} matrix of spatial
#'   coordinates. Is missing, then \code{coords2 = coords1}.
#' @inheritParams angle2d
#'
#' @return A \code{ganisoD} object with components \code{d}
#'   and \code{angles}, which is the distance matrix between
#'   the coordinates and the angles between the coordinates.
#'   The angles are returned in radians.
#' @export
#' @examples
#' ganiso_d(cbind(0, 0), cbind(1, 1))
ganiso_d = function(coords1, coords2, radians = TRUE, invert = TRUE) {
  if (missing(coords2)) {
    coords2 = coords1
  }
  if (!is.matrix(coords1) | !is.matrix(coords2)) {
    stop("coords1 and coords2 must be matrices")
  }
  if (ncol(coords1) != 2) {
    stop("coords1 should have only 2 columns")
  }
  if (ncol(coords2) != 2) {
    stop("coords2 should have only 2 columns")
  }

  out = vector("list")
  out$d = geodist(coords1, coords2)
  out$angles = t(apply(coords1, 1, function(x) {
    angle2d(matrix(x, nrow = nrow(coords2), ncol = 2, byrow = TRUE), coords2, radians = TRUE, invert = invert)
    }))
  out$radians = radians
  out$invert = invert
  class(out) = "ganisoD"
  return(out)
}
