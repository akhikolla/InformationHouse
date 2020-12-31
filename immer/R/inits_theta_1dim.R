## File Name: inits_theta_1dim.R
## File Version: 0.13

#####################################################################
inits_theta_1dim <- function( dat, pid, eps=.05, sd_init=1 )
{
    N <- length(unique(pid))
    # initial values
    maxK <- apply( dat, 2, max, na.rm=TRUE )
    I <- ncol(dat)
    K <- max(maxK)
    # ability inits
    theta <- stats::aggregate( dat, list(pid), mean, na.rm=TRUE )[,-1]
    theta <- theta / matrix( maxK, nrow=N, ncol=I, byrow=TRUE )
    theta <- rowMeans(theta, na.rm=TRUE)
    theta <- ( theta + eps ) / ( 1 + 2*eps )
    theta <- stats::qlogis(theta)
    theta <- theta - mean(theta)
    theta <- sd_init / stats::sd(theta) * theta
    return(theta)
}
#####################################################################
