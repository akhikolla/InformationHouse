## File Name: mnlfa_penalty_values.R
## File Version: 0.16


mnlfa_penalty_values <- function(parms, parms_indices, regular_type, regular_lam, N_item,
    center_group_parms )
{
    x <- parms[ parms_indices ]
    NI <- length(x)
    val <- rep(0,NI)
    regularized <- rep(0,NI)
    estimated <- rep(1,NI)
    if (regular_type !="none"){
        if ( NI > 1 ){
            x_norm <- sqrt(sum(x^2))
            if (center_group_parms){
                lam_fac <- sqrt(NI-1)
            } else {
                lam_fac <- sqrt(NI)
            }
            regular_lam <- lam_fac*regular_lam
            x <- x_norm
        }
        val <- cdm_penalty_values(x=x, regular_type=regular_type, regular_lam=regular_lam,
                    regular_tau=NULL, regular_alpha=NULL)
        val <- N_item * val
        #** number of regularized parameters
        eps <- 1E-5
        regularized <- sum( ( val < eps ) * ( regular_lam > 0 ) )
        estimated <- sum( val >=eps )
        if (NI > 1 ){
            val <- rep( val / NI, NI )
            if ( sum(regularized) > 0 ){
                regularized <- c( rep(regularized, NI-1), 0)
            } else {
                regularized <- rep(0, NI)
            }
            if ( sum(estimated) > 0 ){
                estimated <- c( rep(1, NI-1), 0)
            } else {
                estimated <- rep(0, NI)
            }
        }
    }
    #--- output
    res <- list( val=val, regularized=regularized, estimated=estimated )
    return(res)
}
