## File Name: immer_jml_facets_update_step_effects.R
## File Version: 0.11

immer_jml_facets_update_step_effects <- function(parm, des_names1, parm_sign,
    suff_stat, design, N, K, maxcat, weights, max_incr, center, max_step,
    se, parm_change)
{
    parm0 <- parm
    dd <- "step"
    parm0_dd <- parm[[dd]]
    parm_sign_dd <- parm_sign[[dd]]
    score_dd <- suff_stat[[dd]]
    design_dd <- design[,dd]
    for (mm in 1:maxcat){
        parm_dd <- parm[[dd]]
        res <- immer_jml_facets_calc_probs( maxcat=maxcat, N=N, K=K, design=design,
                        des_names1=des_names1, parm=parm, parm_sign=parm_sign, is_step=TRUE )
        probs <- res$probs
        M <- res$M
        Var <- res$Var
        se_step <- 0*parm0_dd
        parm_dd <- parm[[dd]]
        p_mm <- probs[,seq(mm+1,maxcat+1),drop=FALSE]
        d1 <- parm_sign_dd * ( score_dd[,mm] - rowsum( weights * rowSums(p_mm), design_dd )[,1] )
        d2 <- abs( rowsum( weights * rowSums(p_mm - ( rowSums(p_mm)^2 ) ), design_dd )[,1] )
        incr <- d1 / d2
        incr <- immer_trim_increment(incr=incr, max_incr=max_incr)
        parm_dd[,mm] <- parm_dd[,mm] + incr
        se_step[,mm] <- sqrt( 1 / d2 )
        parm[[dd]] <- parm_dd
    }
    if (center[[dd]]){
        parm_dd <- parm[[dd]]
        parm_dd <- parm_dd - rowMeans( parm_dd)
        parm[[dd]] <- parm_dd
    }
    se[[dd]] <- se_step
    parm_change[[dd]] <- max( abs( parm0[[dd]] - parm[[dd]] ) )
    #--- output
    res <- list(parm=parm, parm_change=parm_change, se=se, probs=probs,
                    M=M, Var=Var)
    return(res)
}
