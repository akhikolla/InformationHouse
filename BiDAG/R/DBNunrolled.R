#' An unrolled adjacency matrix of a dynamic Bayesian network 
#'
#' An adjacency matrix representing the ground truth DBN used to generate a synthetic dataset \code{\link{DBNdata}}. The matrix is an unrolled representation of a 2-step DBN, such that the static variables are represented in the first 3 columns/rows of the matrix.
#'
#' @format A binary matrix with 63 rows and 63 columns representing an adjacency matrix of a DBN. Rows and columns of the matrix correspond to 15 variables (s1, s2, s3, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12) of a DBN across 5 time slices.
#'
"DBNunrolled"