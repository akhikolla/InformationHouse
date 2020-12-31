## File Name: mnlfa_mstep_item_loglike_2pl_deriv.R
## File Version: 0.186


mnlfa_mstep_item_loglike_2pl_deriv <- function(y, y_resp, theta, parms, Xdes_int, Xdes_slo, post,
        b_index, a_index, parms_indices, h, N_item, eps=1E-15 )
{
    args <- list( y=y, y_resp=y_resp, theta=theta, parms=parms, Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, post=post,
                b_index=b_index, a_index=a_index, N_item=N_item, eps=eps )
    fct <- 'mnlfa_mstep_item_loglike_2pl_eval'
    ll0 <- do.call( what=fct, args)
    NP <- length(parms_indices)
    ll1 <- rep(NA,NP)
    ll2 <- rep(NA,NP)
    parms0 <- parms
    for (pp in 1:NP){
        parms1 <- parms0
        parms1[ parms_indices[pp] ] <- parms0[ parms_indices[pp] ] + h
        args$parms <- parms1
        ll1[pp] <- do.call( what=fct, args)
        parms1[ parms_indices[pp] ] <- parms0[ parms_indices[pp] ] - h
        args$parms <- parms1
        ll2[pp] <- do.call( what=fct, args)
    }
    #-- compute derivatives
    res <- mnlfa_differences_derivative(ll0=ll0, ll1=ll1, ll2=ll2, h=h)
    D1 <- res$D1
    D2 <- res$D2
    D2_max <- res$D2_max
    #-- updates
    incr <- D1 / D2_max
    #-- output
    res <- list( ll=ll0, D1=D1, D2=D2, D2_max=D2_max, incr=incr)
    return(res)
}
