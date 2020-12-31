## File Name: mnlfa_proc_theta_design.R
## File Version: 0.02

mnlfa_proc_theta_design <- function(theta)
{
    if ( is.null(theta) ){
        theta <- matrix( seq(-5,5,len=15), ncol=1 )
        colnames(theta) <- "theta1"
    }
    TP <- nrow(theta)
    if ( ncol(theta) > 1 ){
        stop("Multidimensional models not yet implemented.\n")
    }
    #-- outut
    res <- list(theta=theta, TP=TP)
    return(res)
}
