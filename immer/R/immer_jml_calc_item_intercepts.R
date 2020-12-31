## File Name: immer_jml_calc_item_intercepts.R
## File Version: 0.05

immer_jml_calc_item_intercepts <- function(A, xsi, b_fixed=0)
{
    K <- dim(A)[2]
    b <- matrix( 0, nrow=dim(A)[1], ncol=K )
    b <- b + b_fixed
    for (kk in 1:K){
        b[,kk] <- b[,kk] + A[,kk,] %*% xsi
    }
    return(b)
}
