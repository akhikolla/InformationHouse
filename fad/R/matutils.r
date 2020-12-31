
# A funciton to check if a matrix is from package Matrix
is.Matrix <- function(z) inherits(z,"sparseMatrix")


# SD of each column for a sparse matrix. Currently only CSC matrices are supported.
colSD <- function(X,means=NULL)
{
  stopifnot(class(X) == "dgCMatrix")
  if(is.null(means)) means <- colMeans(X)

  v <- colMeans(X^2) - means^2

  if(any(v <= 0)) stop("Some variables have zero variance or are badly scaled.")

  return(sqrt(v))
}

