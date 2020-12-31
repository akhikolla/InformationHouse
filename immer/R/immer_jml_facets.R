## File Name: immer_jml_facets.R
## File Version: 0.26


immer_jml_facets <- function(y, design, center=NULL, weights=NULL, max_equal=FALSE,
        eps=.3,    bc=TRUE, max_incr=1, max_iter=500)
{
    CALL <- match.call()
    s1 <- Sys.time()

    #-- data processing for design
    res <- immer_jml_facets_proc_data( y=y, design=design, weights=weights,
                max_equal=max_equal, eps=eps, bc=bc, center=center )
    design <- res$design
    design_labels <- res$design_labels
    N <- res$N
    weights <- res$weights
    data <- res$data
    parm <- res$parm
    suff_stat <- res$suff_stat
    parm_sign <- res$parm_sign
    des_names <- res$des_names
    center <- res$center
    maxcat <- res$maxcat
    K <- res$K
    des_names1 <- res$des_names1
    is_step <- res$is_step
    max_item <- res$max_item
    max_step <- res$max_step
    y_design <- res$y_design


    #** algorithmic defaults
    iter <- 0

    #** start algorithm
    while( iter < max_iter){

        #* update main effects
        res <- immer_jml_facets_update_main_effects( parm=parm, des_names1=des_names1,
                    parm_sign=parm_sign, suff_stat=suff_stat, design=design, N=N, K=K, maxcat=maxcat,
                    is_step=is_step, weights=weights, max_incr=max_incr, center=center )
        parm <- res$parm
        parm_change <- res$parm_change
        se <- res$se
        probs <- res$probs
        M <- res$M
        Var <- res$Var

        #* update step parameters
        if (is_step){
            res <- immer_jml_facets_update_step_effects( parm=parm, des_names1=des_names1,
                        parm_sign=parm_sign, suff_stat=suff_stat, design=design, N=N, K=K, maxcat=maxcat,
                        weights=weights, max_incr=max_incr, center=center, max_step=max_step, se=se,
                        parm_change=parm_change )
            parm <- res$parm
            parm_change <- res$parm_change
            se <- res$se
            probs <- res$probs
            M <- res$M
            Var <- res$Var
        }

        # compute log-likelihood
        eps <- 1e-50
        ll <- sum( weights * log( probs[ y_design ] + eps ) )
# Revalpr("ll")
        parm_change <- abs( unlist( parm_change ) )
# Revalpr("parm_change")
        iter <- iter + 1

        utils::flush.console()


    }



}
