## File Name: immer_jml_ic.R
## File Version: 0.05

#######################################################
# information criteria immer_jml
immer_jml_ic <- function( loglike, N, center_theta, xsi, I )
{
    ic <- list( dev=-2*loglike)
    ic$n <- N
    ic$I <- I
    ic$ntheta <- N - center_theta
    ic$nxsi <- length(xsi)
    Npars <- ic$ntheta + ic$nxsi
    ic$Npars <- Npars
    ic$np <- ic$Npars
    ic$ND <- NA
    ic$R <- NA
    ic <- immer_IC_calc(ic=ic)
    return(ic)
}
############################################################
