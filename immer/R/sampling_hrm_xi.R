## File Name: sampling_hrm_xi.R
## File Version: 3.48

######################################################
# sample complete xi data frame
sampling_hrm_xi <- function( dat, theta, b, a, phi, psi, K, pid, rater,  ND,
            dat_ind, N, I, maxK, useRcpp, xi_ind )
{
    xi <- matrix( NA, nrow=N, ncol=I)
    eps <- 1E-20
    for (ii in 1:I){
        x <- dat[,ii]
        x_ind <- dat_ind[,ii]
        xi[,ii] <- sampling_hrm_xi_item( x=x, theta=theta, b=b[ii,], a=a[ii],
                        phi=phi[ii,], psi=psi[ii,], K=maxK[ii], pid=pid,
                        rater=rater, ND=ND, x_ind=x_ind, useRcpp, xi_ind=xi_ind[,ii])
    }
    return(xi)
}
##############################################################

