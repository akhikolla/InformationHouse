#' Double-Centre a Three-way Array
#' 
#' Double-centre the matrix slices of a three-way array.
#' 
#' @param array A three-way array
#' @param margin The way of the array over which the centring must be done
#' @param rows Logical indicating whether to centre the rows of the matrix slices
#' @param columns Logical indicating whether to centre the columns of the matrix slices
#' 
#' @export
carray <- function(array, margin = 3L, rows = TRUE, columns = TRUE) {
  
  ## Checks
  if (length(dim(array)) != 3L)
    stop("Argument 'array' must be a three-way array.")
  if (length(margin) != 1L) {
    warning("Only the first element of argument 'margin' is used.")
    margin <- margin[1L]
  }
  
  ## Indices and dimensions
  ind <- seq_len(3L)
  ind_rm <- ind[-margin]
  row_ind <- ind_rm[1L]
  col_ind <- ind_rm[2L]
  dims <- dim(array)
  dims_rm <- dims[-margin]
  row_dim <- dims_rm[1L]
  col_dim <- dims_rm[2L]
  
  ## Centre the rows
  if (rows) {
    rmns <- apply(array, sort(c(row_ind, margin)), mean)
    rmns <- rep(rmns, each = col_dim)
    rmns <- array(rmns, dim = c(col_dim, dims[-col_ind]))
    perm <- switch(margin, c(2:3, 1), c(2, 3, 1), c(2:1, 3))
    rmns <- aperm(rmns, perm)
    array <- array - rmns
  }
  
  ## Centre the columns
  if (columns) {
    cmns <- apply(array, sort(c(col_ind, margin)), mean)
    cmns <- rep(cmns, each = row_dim)
    cmns <- array(cmns, dim = c(row_dim, dims[-row_ind]))
    perm <- switch(margin, c(2, 1, 3), 1:3, 1:3)
    cmns <- aperm(cmns, perm)
    array <- array - cmns
  }
  
  return(array)
}