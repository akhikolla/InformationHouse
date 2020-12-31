## File Name: immer_jml_update_theta_R.R
## File Version: 0.30

immer_jml_update_theta_R <- function(score_pers, I, K, N, theta, b, dat_resp, maxiter_update,
    conv_update, center_theta, max_incr, shortcut_index, weights )
{
    KM <- matrix( 0:K, nrow=N, ncol=K+1, byrow=TRUE )
    iterate <- TRUE
    eps <- 1E-7
    iter <- 0
    probs <- array( 0, dim=c(N,I,K+1) )
    while(iterate){
        theta0 <- theta
        der1 <- score_pers
        der2 <- rep(0,N)
        for (ii in 1:I){
            b_ii <- c(0,b[ii,])
            probs_ii <- probs_pcm_one_item(theta=theta, b_ii=b_ii)
            probs[,ii,] <- probs_ii
            M_ii <- rowSums( KM * probs_ii * dat_resp[,ii] )
            Var_ii <- rowSums( KM^2 * probs_ii * dat_resp[,ii] ) - M_ii^2
            der1 <- der1 - M_ii
            der2 <- der2 - Var_ii
        }
        #-- update theta
        incr <- der1 / ( abs(der2) + eps )
        incr <- immer_trim_increment(incr=incr, max_incr=3)
        theta <- theta + incr
        iter <- iter + 1
        theta_change <- max( abs(theta - theta0) )
        if (iter > maxiter_update){ iterate <- FALSE }
        if (theta_change < conv_update){ iterate <- FALSE }
    }
    theta <- immer_jml_center_theta( theta=theta, center_theta=center_theta )
    #--- output
    res <- list(theta=theta, theta_der2=der2, probs=probs)
    return(res)
}
