## File Name: mnlfa_parameter_regularization.R
## File Version: 0.15


mnlfa_parameter_regularization <- function( parms, parms_indices, L, regular_type,
    regular_lam, regular_tau, regular_alpha, center_group_parms=FALSE )
{
    if (regular_type !="none"){
        x <- parms[ parms_indices ]
        NP <- length(parms_indices)
        if ( NP==1 ){
            # single parameter regularization
            x <- x*L
            x2 <- CDM::cdm_parameter_regularization(x=x, regular_type=regular_type, regular_lam=regular_lam,
                            regular_tau=regular_tau, regular_alpha=regular_alpha)
            x2 <- x2 / L
        } else {
            # group parameter regularization
            if (center_group_parms){
                x <- x - mean(x)
            }
            x_norm <- sqrt(sum(x^2))
            if (center_group_parms){
                lam_fac <- sqrt(NP-1)
            } else {
                lam_fac <- sqrt(NP)
            }
            regular_lam <- lam_fac*regular_lam
            x2 <- x_norm*L
            x2 <- CDM::cdm_parameter_regularization(x=x2, regular_type=regular_type, regular_lam=regular_lam,
                            regular_tau=regular_tau, regular_alpha=regular_alpha)
            x2 <- x2 / L
            x2 <- x2 * x / x_norm
        }
        parms[ parms_indices ] <- x2
    }
    #-- output
    return(parms)
}
