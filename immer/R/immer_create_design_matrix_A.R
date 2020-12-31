## File Name: immer_create_design_matrix_A.R
## File Version: 0.16

immer_create_design_matrix_A <- function( maxK, A, b_fixed, irtmodel, NA_val=99)
{
    I <- length(maxK)
    K <- max(maxK)
    if ( ! is.null(A) ){
        irtmodel <- "user"
    }

    #--- create design matrix A
    if ( is.null(A) ){
        if (irtmodel=="PCM"){
            A <- immer_create_design_matrix_A_PCM(maxK=maxK, I=I, K=K)
        }
        if (irtmodel=="PCM2"){
            A <- immer_create_design_matrix_A_PCM2(maxK=maxK, I=I, K=K)
        }

    }

    if ( is.null( dimnames(A) ) ){
        dimnames(A)[[1]] <- names(maxK)
        dimnames(A)[[2]] <- paste0("Cat", 1:K)
        dimnames(A)[[3]] <- paste0("parm", seq(1,(dim(A))[3]) )
    }

    #--- fixed b values
    if (is.null(b_fixed)){
        b_fixed <- matrix( 0, nrow=I, ncol=K)
        rownames(b_fixed) <- dimnames(A)[[1]]
        colnames(b_fixed) <- dimnames(A)[[2]]
        for (ii in 1:I){
            if ( maxK[ii] < K){
                b_fixed[ii, seq(maxK[ii]+1,K) ] <- NA_val
            }
        }

    }

    #--- output
    res <- list( A=A, b_fixed=b_fixed, irtmodel=irtmodel )
    return(res)
}
