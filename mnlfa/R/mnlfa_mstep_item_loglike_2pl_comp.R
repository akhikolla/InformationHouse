## File Name: mnlfa_mstep_item_loglike_2pl_comp.R
## File Version: 0.07

mnlfa_mstep_item_loglike_2pl_comp <- function(y, y_resp, theta, b_ii, a_ii, Xdes_int, Xdes_slo,
        post, N_item, eps=1E-15 )
{
    b <- Xdes_int %*% b_ii
    a <- Xdes_slo %*% a_ii
    like_ii <- mnlfa_rcpp_calc_probs_2pl(a=a, b=b, theta=theta, y=y, y_resp=y_resp )
    log_probs <- log(like_ii+eps)
    ll <- sum( post * log_probs )
    ll <- ll / N_item
    return(ll)
}
