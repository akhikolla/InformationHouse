## File Name: IRT.likelihood_immer.R
## File Version: 0.03


IRT.likelihood.immer_jml <- function(object, theta=seq(-9,9,len=41), ...)
{
    b <- cbind(0,object$b)
    dat <- object$dat
    dat_resp <- object$dat_resp
    I <- ncol(dat)
    a <- rep(1,I)
    #-- compute individual likelihood
    like <- immer_IRT_likelihood_gpcm_helper(theta=theta, b=b, a=a, dat=dat,
                dat_resp=dat_resp)
    return(like)
}
