## File Name: immer_jml_facets_update_main_effects.R
## File Version: 0.03

immer_jml_facets_update_main_effects <- function(parm, des_names1, parm_sign,
    suff_stat, design, N, K, maxcat, is_step, weights, max_incr, center)
{
    parm0 <- parm
    se <- list()
    parm_change <- list()
    for (dd in des_names1){
        #dd <- des_names1[2]
        parm0_dd <- parm[[dd]]
        parm_sign_dd <- parm_sign[[dd]]
        score_dd <- suff_stat[[dd]]$score
        design_dd <- design[,dd]
        res <- immer_jml_facets_calc_probs( maxcat=maxcat, N=N, K=K, design=design,
                    des_names1=des_names1, parm=parm, parm_sign=parm_sign, is_step=is_step )
        probs <- res$probs
        M <- res$M
        Var <- res$Var
        # first derivative
        d1 <- parm_sign_dd * ( score_dd - rowsum( weights * M, design_dd )[,1] )
        # second derivative
        d2 <- abs( rowsum( weights * Var, design_dd )[,1] )
        incr <- d1 / d2
        incr <- immer_trim_increment(incr=incr, max_incr=max_incr)
        parm[[ dd ]] <- parm[[ dd ]] + incr
        if (center[[dd]]){
            parm[[dd]] <- immer_center(x=parm[[dd]])
        }
        se[[dd]] <- sqrt( 1 / d2 )
        parm_change[[dd]] <- max( abs( parm0[[dd]] - parm[[dd]] ) )
    }

    #--- output
    res <- list(parm=parm, parm_change=parm_change, se=se, probs=probs,
                    M=M, Var=Var)
    return(res)
}
