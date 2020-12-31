## File Name: mnlfa_mstep_item_2pl.R
## File Version: 0.510


mnlfa_mstep_item_2pl <- function(y, y_resp, theta, parms, Xdes_int,
        Xdes_slo, post,    b_index, a_index, parms_iterations, h, N_item,
        parms_regular_types, parms_regular_lam, center_group_parms, msteps,
        conv_mstep, eps=1E-15, L_max=.25 )
{
    NH <- length(parms_iterations)
    iterate <- TRUE
    iter <- 1
    fct_mstep <- 'mnlfa_mstep_item_loglike_2pl_update_parameters'

    #--- begin iterations M-steps
    while(iterate){
        parms0 <- parms
        for (hh in seq_len(NH) ){
            parms_indices <- parms_iterations[[hh]]
            regular_type <- parms_regular_types[parms_indices][1]
            regular_lam <- parms_regular_lam[parms_indices][1]
            args <- list( y=y, y_resp=y_resp, theta=theta, parms=parms, Xdes_int=Xdes_int,
                        Xdes_slo=Xdes_slo, post=post, b_index=b_index,
                        a_index=a_index, parms_indices=parms_indices, h=h, N_item=N_item,
                        regular_type=regular_type, regular_lam=regular_lam,
                        center_group_parms=center_group_parms, eps=eps, L_max=L_max )
            res <- do.call( what=fct_mstep, args=args )
            parms <- res$parms
        }
        change <- max( abs( parms - parms0 ) )
        iter <- iter + 1
        iterate <- ! ( ( iter > msteps ) | ( change < conv_mstep ) )
    }
    #--- end iterations M-steps

    #*** penalty values
    res <- mnlfa_penalty_values_item( parms=parms, parms_iterations=parms_iterations,
                parms_regular_types=parms_regular_types, parms_regular_lam=parms_regular_lam,
                N_item=N_item, center_group_parms=center_group_parms )
    parms_values <- res$parms_values
    parms_regularized <- res$parms_regularized
    parms_estimated <- res$parms_estimated

    #--- output
    res <- list( parms=parms, iter_mstep=iter - 1, parms_change=change, parms_values=parms_values,
                    parms_regularized=parms_regularized, parms_estimated=parms_estimated )
    return(res)
}
