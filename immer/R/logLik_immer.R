## File Name: logLik_immer.R
## File Version: 0.14

###############################################################
# log-likelihood immer_hrm
logLik.immer_hrm <- function (object, ...)
{
    out <- object$like
    attr(out, "df") <- sum(object$ic$Npars)
    attr(out, "nobs") <- object$ic$N
    class(out) <- "logLik"
    return(out)
}
#################################################################


###############################################################
# log-likelihood lc2_agreement
logLik.lc2_agreement <- function (object, ...)
{
    out <- object$loglike
    attr(out, "df") <- object$model_output$npars
    attr(out, "nobs") <- object$nobs
    class(out) <- "logLik"
    return(out)
}
#################################################################


###############################################################
#--- logLik immer_cml
logLik.immer_cml <- function (object, ...)
{
    out <- object$loglike
    attr(out, "df") <- sum(object$npars)
    attr(out, "nobs") <- object$N
    class(out) <- "logLik"
    return(out)
}
#################################################################


#--- logLik immer_latent_regression
logLik.immer_latent_regression <- function (object, ...)
{
    out <- -object$deviance / 2
    attr(out, "df") <- object$ic$np
    attr(out, "nobs") <- object$N
    class(out) <- "logLik"
    return(out)
}



#--- logLik immer_jml
logLik.immer_jml <- function (object, ...)
{
    out <- -object$ic$dev / 2
    attr(out, "df") <- object$ic$np
    attr(out, "nobs") <- object$n
    class(out) <- "logLik"
    return(out)
}
