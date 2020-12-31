## File Name: mnlfa_proc_inits_prior.R
## File Version: 0.02

mnlfa_proc_inits_prior <- function(theta, N)
{
    TP <- nrow(theta)
    w1 <- stats::dnorm( theta[,1] )
    pi.k <- w1 / sum(w1)
    prior <- matrix( pi.k, nrow=N, ncol=TP, byrow=TRUE)
    return(prior)
}
