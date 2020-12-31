## File Name: immer_mcmc_draw_variance.R
## File Version: 0.09


###################################################################
# draw variances from inverse chi square distribution
immer_mcmc_draw_variance <- function( N, w0, sig02, n, sig2 )
{
    # INPUT:
    # N ... number of random draws
    # w0 ... sample size prior
    # sig02 ... prior variance
    # n ... empirical sample size
    # sig2 ... empirical variance
    res <- 1/ stats::rgamma( N, 0.5 * (w0+n), 0.5 * ( w0*sig02 + n*sig2 ) )
    return(res)
}
#####################################################################
