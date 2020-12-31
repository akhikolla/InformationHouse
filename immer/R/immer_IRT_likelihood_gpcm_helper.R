## File Name: immer_IRT_likelihood_gpcm_helper.R
## File Version: 0.06

immer_IRT_likelihood_gpcm_helper <- function(theta, b, a, dat, dat_resp=NULL)
{
    dat <- as.matrix(dat)
    if ( is.null(dat_resp) ){
        dat_resp <- 1 - is.na(dat)
    }
    dat_resp <- as.matrix(dat_resp)
    dat[ is.na(dat) ] <- 0
    b <- as.matrix(b)
    K <- ncol(b)
    TP <- length(theta)
    #-- compute probabilities
    probs <- immer_gpcm_calc_probs( theta=theta, b=b, a=a )
    #-- compute individual likelihood
    like <- immer_irt_likelihood_gpcm( probs=probs, dat=dat, dat_resp=dat_resp, TP=TP, K=K )
    attr(like, "theta") <- theta
    attr(like, "prob.theta") <- NA
    attr(like, "skillspace") <- NA
    attr(like, "G") <- NA
    return(like)
}
