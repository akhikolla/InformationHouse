#' Invert a symmetric positive-definite matrix
#'
#' Inverts a symmetric positive-definite matrix without requiring the Matrix package.
#'
#' @param obj The symmetric positive-definite matrix
#' @return A matrix of the same size as the input object
#' @details This function does no error checking and it is up to the user to ensure that the input is indeed symmetric, positive-definite, and a matrix.
#'
#' @examples
#' m <- matrix(c(1, 0, .5, .5, 0, 1, .5, .5, .5, .5, 1, .5, .5, .5, .5, 1), 4)
#' sinv(m)
#'
#' @author Claus Ekstrom, \email{claus@@rprimer.dk}.
#' @keywords file
#' @export

sinv <- function(obj) {
    return(chol2inv(chol(obj)))
}
