#' Apply getElement over a list
#'
#' \code{lget} or \code{lgetElement} applies
#' \code{getElement} to a
#' list using \code{\link[base]{lapply}}. \code{sget} and
#' \code{sgetElement} do the same thing with
#' \code{sapply}.
#'
#' @param X A list.
#' @inheritParams base::Extract
#' @inheritParams base::lapply
#' @return A list (\code{lget}) or vector (\code{sget})
#' of the same length as \code{X} with the
#' \code{name} parts of each element of \code{X}.
#' @export
#' @examples
#' e1 = list(x = rnorm(5),
#'           y = letters[c(1:2, 2:1, 3)],
#'           z = c(TRUE, TRUE, FALSE, TRUE, TRUE)
#' )
#' e2 = list(x = rnorm(5),
#'           y = letters[c(1:4, 1)],
#'           z = c(FALSE, TRUE, FALSE, TRUE, FALSE))
#' X = list(e1, e2)
#' lget(X, name = "x")
#' sget(X, name = "y")
lget = function(X, name) {
  base::lapply(X, FUN = getElement, name = name)
}

#' @rdname lget
#' @export
lgetElement = lget

#' @rdname lget
#' @export
sget = function(X, name, simplify = TRUE, USE.NAMES = TRUE) {
  base::sapply(X, FUN = getElement, name = name,
               simplify = simplify, USE.NAMES = USE.NAMES)
}

#' @rdname lget
#' @export
sgetElement = sget
