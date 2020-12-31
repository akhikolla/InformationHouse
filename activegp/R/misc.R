#' Get the distance between subspaces defined as the ranges of A and B
#'
#' @param A A matrix or const_C object.
#' @param B Another matrix with the same number of rows as A, or const_C object of the same dimension.
#' @param r A scalar integer, the dimension of the subspace to compare (only necessary if either A or B is a const_C object).
#' @return A nonnegative scalar giving the cosine of the first principle angle between the two subspaces.
#' @export
#' @importFrom stats rnorm
#' @importFrom methods is
subspace_dist <- function(A, B, r) {
    if (is(A, "numeric")) {
        A <- matrix(A, ncol = 1)
    } else if (is(A, "const_C")) {
      if (missing(r)) {
        stop("Need to specify r when passing a const_C object")
      }
      A <- matrix(eigen(A$mat)$vectors[,1:r], ncol = r)
    }
    if (is(B, "numeric")) {
        B <- matrix(B, ncol = 1)
    } else if (is(B, "const_C")) {
      if (missing(r)) {
        stop("Need to specify r when passing a const_C object")
      }
      B <- matrix(eigen(B$mat)$vectors[,1:r], ncol = r)
    }
    if(!all(dim(A) == dim(B))) stop("Matrices not of same dimension")
  
    m <- nrow(A)
    #Normalize columns
    # Complete A's columns randomly to form a basis for Rn
    A_b <- cbind(A, matrix(rnorm(m*(m-ncol(A))), ncol = m-ncol(A)))
    # An Orthonormal basis for A's complement.
    A_c <- qr.Q(qr(A_b))[,(ncol(A)+1):m]
    # An orthonormal basis for B
    B_b <- qr.Q(qr(B))

    return(norm(t(A_c) %*% B_b, '2'))
}

#' Change a function's inputs to live in [-1, 1]
#' 
#' Given an m dimensional function whose inputs live in bounded intervals [a1, b1], ..., [am, bm], return a wrapped version of the function whose inputs live in [-1, 1], ..., [-1, 1].
#' @param f The function to wrap, should have a single vector-valued input.
#' @param domain A list of real tuples, indicating the original domain of the function.
#' @return A function wrapping f.
#' @export
domain_to_unit <- function(f, domain) {
    n11_2_ab <- function(x, a, b) {
        a + (b - a) * (x + 1) / 2
    }
    function(x) {
        xt <- sapply(1:length(x), function(i) n11_2_ab(x[i], domain[[i]][1], domain[[i]][2]))
        f(xt)
    }
}

#' f:[-1, 1] -> R Becomes f:[0,1] -> R
#' @param f initial function
#' @return The same function with domain shifted.
#' @export
n11_2_01 <- function(f) {
    function(x) f(2*(x - 0.5))
}

#' Rectangular Domain -> Unbounded Domain
#' 
#' Given an m dimensional function whose inputs live in bounded intervals [a1, b1], ..., [am, bm], return a wrapped version of the function whose inputs live in R^m. Transformed using the logit function.
#' @param f The function to wrap, should have a single vector-valued input.
#' @param domain A list of real tuples, indicating the original domain of the function.
#' @return A function wrapping f.
#' @export
#' @importFrom stats plogis
domain_to_R <- function(f, domain) {
    R_2_ui <- plogis
    ui_2_ab <- function(x, a, b) (b-a) * x + a
    trans <- function(x, a, b) ui_2_ab(R_2_ui(x), a, b)
    function(x) {
        xt <- sapply(1:length(x), function(i) trans(x[i], domain[[i]][1], domain[[i]][2]))
        f(xt)
    }
}

#' Estimate the Active Subspace of a Cheap Function using Gradients
#'
#' Looks between [-1, 1]
#'
#' @param f The function to eval
#' @param r The max dim of the active subspace
#' @param m The dimension of the underlying/embedding space.
#' @param M optional budget of evaluations, default to \code{2 * r * log(m)}
#' @param scale Scale all gradients to have norm 1?
#' @return A list with sub, the active subspace, sv, the singular values (all m of them), fs, which gives function values, gs, function grads, and X, which gives sampled locations.
#' @export
#' @importFrom stats runif
#' @importFrom numDeriv grad
grad_est_subspace <- function(f, r, m, M = NULL, scale = FALSE) {
    if (is.null(M)) {
        M <- ceiling(2 * r * log(m))
    }

    X <- matrix(runif(m*M, -1, 1), ncol = M)#Each column is a point we're going to evaluate the func at.
    grads <- sapply(1:M, function(i) grad(f, X[,i], method = 'simple'))
    if (scale) {
        grads <- apply(grads, 2, function(i) i / sqrt(sum((i)^2)))
    }
    decomp <- svd(grads)

    return(list(sub = matrix(decomp$u[,1:r], ncol = r), sv = decomp$d, gs = grads, X = X))
}
