#' An adjacency matrix of a dynamic Bayesian network 
#'
#' An adjacency matrix representing the ground truth DBN used to generate a synthetic dataset \code{\link{DBNdata}}. The matrix is a compact representation of a 2-step DBN, such that initial structure is stored in the first 15 columns of the matrix and transitional structure is stored in the last 12 columns of the matrix.
#'
#' @format A binary matrix with 27 rows and 27 columns representing an adjacency matrix of a DBN. Rows and columns of the matrix correspond to 15 variables of a DBN across 2 time slices.
#'
"DBNmat"