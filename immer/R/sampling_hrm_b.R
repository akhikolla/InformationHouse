## File Name: sampling_hrm_b.R
## File Version: 0.26

#####################################################################
# sampling of item parameters
sampling_hrm_b <- function( xi, xi_ind, b, a, maxK, prior, MHprop, I,
        theta, useRcpp, eps=1E-20 )
{
    # refresh count
    MHprop$refresh_count$b <- MHprop$refresh_count$b + 1
    for (ii in 1:I ){
        for (kk in seq(1,maxK[ii]) ){
            b_new <- b_old <- b
            b_new[ii,kk] <- stats::rnorm( 1, mean=b_old[ii,kk], sd=MHprop$SD$b[ii,kk] )
            p_new <- stats::dnorm( b_new[ii,kk], mean=prior$b$M[ii,kk], sd=prior$b$SD[ii,kk] )
            p_old <- stats::dnorm( b_old[ii,kk], mean=prior$b$M[ii,kk], sd=prior$b$SD[ii,kk] )
            ll_new <- sum( probs_gpcm( x=xi[,ii], theta, b=b_new[ii,],
                            a=a[ii], K=maxK[ii], x_ind=xi_ind[,ii],
                            useRcpp=useRcpp, use_log=TRUE  ) )
            ll_old  <- sum( probs_gpcm( x=xi[,ii], theta, b=b_old[ii,],
                            a=a[ii], K=maxK[ii], x_ind=xi_ind[,ii],
                            useRcpp=useRcpp, use_log=TRUE ) )
            ratio <- immer_sampling_calc_ratio( ll_old=ll_old, ll_new=ll_new,
                            p_old=p_old, p_new=p_new, eps=eps )
            if ( ratio > stats::runif(1) ){
                MHprop$accept$b[ii,kk] <- MHprop$accept$b[ii,kk] + 1
                b <- b_new
            }
        }  # end category kk
    }  # end item ii
    res <- list( b=b, MHprop=MHprop )
    return(res)
}
#####################################################################
