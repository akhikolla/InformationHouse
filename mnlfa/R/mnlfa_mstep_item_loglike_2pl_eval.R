## File Name: mnlfa_mstep_item_loglike_2pl_eval.R
## File Version: 0.09


mnlfa_mstep_item_loglike_2pl_eval <- function(y, y_resp, theta, parms, Xdes_int, Xdes_slo, post,
        b_index, a_index, N_item, eps=1E-15 )
{
    b_ii <- parms[ b_index ]
    a_ii <- parms[ a_index ]
    res <- mnlfa_mstep_item_loglike_2pl_comp( y=y, y_resp=y_resp, theta=theta,
                b_ii=b_ii, a_ii=a_ii, Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, post=post, N_item=N_item,
                eps=eps )
    return(res)
}
