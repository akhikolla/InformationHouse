#' Solve using cholesky decomposition
#'
#' \code{solve_chol} solves a system of equations using the
#' cholesky decomposition of a positive definite matrix
#' \code{A}, i.e., using \code{a = chol(A)}.
#'
#' Note: Unless you have good reason to suspect that the
#' cholesky decomposition of your matrix will be stable, it
#' is recommended that you use \code{\link[base]{solve}} or
#' perform the solve using the \code{\link[base]{qr}}
#' decomposition.
#'
#' @param a The cholesky decomposition of a positive
#'   definite matrix.
#' @param b A numeric or complex vector or matrix giving the
#'   right-hand side(s) of the linear system. If missing,
#'   \code{b} is taken to be an identity matrix and solve
#'   will return the inverse of \code{a}.
#' @param ... Currently unused.
#' @author Joshua French
#' @export
#' @examples
#' set.seed(10)
#' # create positive definite matrix a
#' a = crossprod(matrix(rnorm(25^2), nrow = 25))
#' # create vector x and matrix b
#' # x can be used to check the stability of the solution
#' x = matrix(rnorm(25))
#' b = a %*% x
#'
#' # standard solve
#' x1 = solve(a, b)
#' all.equal(x, x1)
#'
#' # solve using cholesky decomposition
#' chola = chol(a)
#' x2 = solve_chol(chola, b)
#' all.equal(x, x2)
#'
#' # solve using qr decomposition
#' qra = qr(a)
#' x3 = solve.qr(qra, b)
#' all.equal(x, x3)
#'
#' # compare direct inversion
#' ai1 = solve(a)
#' ai2 = solve_chol(chola) #using cholesky decomposition
#' all.equal(ai1, ai2) # should be TRUE
#' @seealso \code{\link[base]{solve}},
#'   \code{\link[base]{chol}}, \code{\link[base]{qr}},
#'   \code{\link[base]{solve.qr}}
#' @export
solve_chol = function(a, b, ...){
  if (missing(b)) {
    chol2inv(a)
  } else {
    forwardsolve(a, backsolve(a, b, transpose = TRUE), upper.tri = TRUE)
  }
}

#' @export
#' @rdname solve_chol
solve.chol = solve_chol
