#' Updates the estimates of the connection signs by running mini OLS models as described in Striaukas and Weber (2018).
#' 
#' @param BX Bx element
#' @param By{By element}
#' @param beta \eqn{\hat{\beta}} estimated value
#' @param xi \eqn{\hat{\xi}} matrix estimated at the previous step
#' @param M Penalty matrix
#' @return xi Updated \eqn{\hat{\xi}} matrix
get.xi <- function(Bx ,By ,beta ,xi ,M){
  
  # check only lower triangular values, since M is symmetric
  M[upper.tri(M)] <- 0
  non.zero <- which(M!=0,arr.ind = T)
  # in case diagonal has a non zero entry, we exlude that
  non.diag <- non.zero[,1]!=non.zero[,2]
  # store indicies of lower triagular signs being not zero
  e.con    <- non.zero[non.diag,]
  # number of iterations needed 
  dd       <- nrow(e.con)
  if (dd!=0){
    for (i in 1:dd){
      beta.exc   <- beta[-c(e.con[i,1],e.con[i,2])]
      beta.tilde <- By[,i] - Bx[,,i]%*%beta.exc
      # updating only lower triangular part
      xi[e.con[i,1],e.con[i,2]] <- sign(beta.tilde[1])*sign(beta.tilde[2])
      
    }
    # updating upper triagular part
    xi[upper.tri(xi)] <- t(xi)[upper.tri(xi)]
    diag(xi) <- 1
  } else {
    xi <- M
  }
  return(xi)
}