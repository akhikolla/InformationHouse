## File Name: probs_hrm_R.R
## File Version: 0.17

######################################################################
probs_hrm_R <- function( x, xi, phi, psi, K, x_ind=NULL )
{
    N <- length(xi)
    KM <- matrix( 0:K, nrow=N, ncol=K+1, byrow=TRUE )
    p1 <- exp( - ( KM - ( xi + phi ) )^2 / ( 2 * psi ) )
    probs <- p1 / rowSums(p1, na.rm=TRUE)
    if ( ! is.null(x) ){
        ind <- cbind( 1:N, x+1 )
        probs <- probs[ind ]
    }
    if ( ! is.null( x_ind) ){
        probs <- ifelse( x_ind==0, 1, probs )
    }
    return(probs)
}
####################################################
