## File Name: sampling_hrm_xi_item.R
## File Version: 0.16


##############################################################
# sampling xi HRM
sampling_hrm_xi_item <- function( x, theta, b, a, phi, psi, K, pid, rater,
        ND=NULL, x_ind=NULL, useRcpp, xi_ind  )
{
    if ( is.null(ND) ){
        ND <- length(x)
    }
    eps <- 1E-20
    N <- length(theta)
    # matrix with probabilities
    probs <- matrix( NA, nrow=N, ncol=K+1)
    # calculate probability GPCM
    for (kk in 0:K){
        # kk <- 0
        p1 <- probs_gpcm( x=rep(kk,N), theta, b=b, a=a, K=K, x_ind=xi_ind,
                    useRcpp=useRcpp )
        p2 <- probs_hrm( x=x, xi=rep( kk, ND ), phi=phi[ rater ], psi=psi[rater ], K=K,
                    x_ind=x_ind, useRcpp=useRcpp)
        p2 <- log( p2 + eps )
        p3 <- rowsum( p2, pid )        # sum over logarithms
        probs[, kk+1] <- p1 * exp( p3[,1] )
    }
    probs <- probs / rowSums( probs )
    rn <- stats::runif(N)
    if ( ! useRcpp ){
        probs1 <- sirt::rowCumsums.sirt(probs)
        xi <- sirt::rowIntervalIndex.sirt(matr=probs1,rn=rn) - 1
    } else {
        # xi <- .Call("sample_prob_index", as.matrix(probs), rn, PACKAGE="immer")
        xi <- sample_prob_index( as.matrix(probs), rn )
    }
    return(xi)
}
################################################################
