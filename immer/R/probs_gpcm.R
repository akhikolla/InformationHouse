## File Name: probs_gpcm.R
## File Version: 2.29

###****************
# probs GPCM
# INPUT:
#  x ... vector of categories
#  theta ... vector of abilities
#  b ... vector of item parameters
#  a ... item discrimination
#  K ... maximum category (1,2,...)
probs_gpcm <- function( x, theta, b, a, K, x_ind=NULL, useRcpp=FALSE, eps=1E-20, use_log=FALSE)
{
    if ( ! useRcpp ){
        probs <- probs_gpcm_R( x, theta, b, a, K, x_ind  )
    }
    if ( useRcpp ){
        if ( is.null(x_ind) ){
            x_ind <- rep(1, length(theta))
        }
        probs <- probs_gpcm_rcpp( x, theta, b, a, K, x_ind )
    }
    if (use_log){
        probs <- log( probs + eps )
    }
    return(probs)
}
##################################################################

