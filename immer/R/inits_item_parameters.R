## File Name: inits_item_parameters.R
## File Version: 0.12

################################################
# inits item parameters
inits_itempars <- function( dat, prior, sd_init=1 )
{
    maxK <- apply( dat, 2, max, na.rm=TRUE )
    K <- max(maxK)
    I <- ncol(dat)
    b <- matrix( NA, nrow=I, ncol=K)
    b1 <- colMeans( dat, na.rm=TRUE ) / maxK
    for (ii in 1:I){
        b[ii, seq(1,maxK[ii],1) ] <- b1[ii] + seq( -2, 2, length=maxK[ii]  )
        # b <- sd_init * b
    }

    if ( ! is.null( prior$b$M ) ){
        b <- prior$b$M
    }
    if ( ! is.null( prior$a$M ) ){
        a <- prior$a$M
    }
    a <- rep(1,I)
    res <- list( b=b, maxK=maxK, K=K, a=a, I=I)
    return(res)
}
#################################################
