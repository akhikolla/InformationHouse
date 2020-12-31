## File Name: mnlfa_mstep_item_loglike_2pl_update_parameters.R
## File Version: 0.275


mnlfa_mstep_item_loglike_2pl_update_parameters <- function(y, y_resp, theta, parms, Xdes_int,
        Xdes_slo, post,    b_index, a_index, parms_indices, h, N_item,
        regular_type, regular_lam, center_group_parms, eps=1E-15, L_max=.25 )
{
    res <- mnlfa_mstep_item_loglike_2pl_deriv( y=y, y_resp=y_resp, theta=theta,
                parms=parms, Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, post=post, b_index=b_index,
                a_index=a_index, parms_indices=parms_indices, h=h, N_item=N_item, eps=eps )
    incr <- res$incr
    L <- max( res$D2_max, L_max)
    ll <- res$ll
    # Newton-Raphson update
    parms[ parms_indices ] <- parms[ parms_indices ] + incr
    # apply threshold operator
    parms <- mnlfa_parameter_regularization( parms=parms, parms_indices=parms_indices,
                L=L, regular_type=regular_type, regular_lam=regular_lam, regular_tau=NULL,
                regular_alpha=NULL, center_group_parms=center_group_parms )
    #--- output
    res <- list( parms=parms, ll=ll)
    return(res)
}
