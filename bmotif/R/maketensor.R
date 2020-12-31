maketensor <- function(A, B){
  x <- dim(A)[1]
  y <- dim(A)[2]
  z <- dim(B)[2]
  C <- array(rep(1, x * y * z), c(x, y, z))
  for (j in 1 : y){
    for (k in 1 : z){
      C[,j,k] <- A[,j] * B[,k]
    }
  }
  C
}