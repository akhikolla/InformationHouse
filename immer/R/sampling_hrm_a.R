## File Name: sampling_hrm_a.R
## File Version: 0.28

#####################################################################
# sampling of item parameters
sampling_hrm_a <- function( xi, xi_ind, b, a, maxK, prior, MHprop, I,
        theta, eps=1E-20, useRcpp=TRUE )
{
    # refresh count
    MHprop$refresh_count$a <- MHprop$refresh_count$a + 1

    for (ii in 1:I ){
        a_new <- a_old <- a
        a_new[ii] <- stats::rlnorm( 1, meanlog=log(a_old[ii]), sdlog=MHprop$SD$a[ii] )
        p_new <- stats::dlnorm( a_new[ii], meanlog=prior$a$M[ii], sdlog=prior$a$SD[ii] )
        p_old <- stats::dlnorm( a_old[ii], meanlog=prior$a$M[ii], sdlog=prior$a$SD[ii] )
        ll_new <- sum( probs_gpcm( x=xi[,ii], theta, b=b[ii,],
                        a=a_new[ii], K=maxK[ii], x_ind=xi_ind[,ii], useRcpp=useRcpp, use_log=TRUE ) )
        ll_old  <- sum( probs_gpcm( x=xi[,ii], theta, b=b[ii,],
                        a=a_old[ii], K=maxK[ii], x_ind=xi_ind[,ii], useRcpp=useRcpp, use_log=TRUE ) )
        #--- distinguish case of fixed and estimated a parameters
        if (prior$a$SD[ii] > 1E-6){
            ratio <- immer_sampling_calc_ratio( ll_old=ll_old, ll_new=ll_new,
                            p_old=p_old, p_new=p_new, eps=eps )
        } else {
            ratio <- 0
        }
        if ( ratio > stats::runif(1) ){
            MHprop$accept$a[ii] <- MHprop$accept$a[ii] + 1
            a <- a_new
        }
    }  # end item ii
    res <- list( a=a, MHprop=MHprop )
    return(res)
}
#####################################################################
