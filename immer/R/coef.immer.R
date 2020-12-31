## File Name: coef.immer.R
## File Version: 0.13

coef.immer_cml <- function(object, ... )
{
    return( object$coefficients )
}

coef.immer_ccml <- function(object, ... )
{
    return( object$coef )
}


coef.immer_latent_regression <- function(object, ... )
{
    return( object$coef )
}
