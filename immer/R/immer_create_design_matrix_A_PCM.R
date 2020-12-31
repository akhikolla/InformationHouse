## File Name: immer_create_design_matrix_A_PCM.R
## File Version: 0.04

immer_create_design_matrix_A_PCM <- function(maxK, I, K)
{
    NX <- sum(maxK)
    A <- array(0, dim=c(I,K,NX) )
    dimnames(A)[[1]] <- names(maxK)
    dimnames(A)[[2]] <- paste0("Cat", 1:K)
    dimnames(A)[[3]] <- paste0("parm", 1:NX)

    pp <- 1
    for (ii in 1:I){
        start_ii <- pp
        for (kk in 1:maxK[ii] ){
            A[ii,kk, seq(start_ii,pp)] <- 1
            dimnames(A)[[3]][pp] <- paste0( names(maxK)[ii], "_Cat", kk )
            pp <- pp + 1
        }
    }
    return(A)
}
