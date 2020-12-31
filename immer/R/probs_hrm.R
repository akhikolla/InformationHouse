## File Name: probs_hrm.R
## File Version: 0.17


#######################################################
# Probabilities in hierarchical rater model
probs_hrm <- function( x, xi, phi, psi, K, x_ind=NULL, useRcpp=FALSE)
{
    if ( ! useRcpp ){
        probs <- probs_hrm_R( x, xi, phi, psi, K, x_ind )
    } else {
        probs <- probs_hrm_rcpp( x, xi, phi, psi, K, x_ind )
    }
    return(probs)
}
####################################################
