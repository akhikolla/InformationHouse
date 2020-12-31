## File Name: immer_jml_update_item_R.R
## File Version: 0.79


immer_jml_update_item_R <- function( score_items, ItemScore, I, K, b, A, xsi, theta,
        N, dat_resp, max_incr,     maxiter_update, conv_update, b_fixed, shortcut_index,
        weights )
{
    iterate <- TRUE
    iter <- 0
    eps <- 1E-7
    NX <- length(xsi)

    r <- matrix(0, nrow=I, ncol=K)
    rr <- array(0, dim=c(I,K,K) )
    probs <- array( 0, dim=c(N,I,K+1) )
    max_incr <- rep(5,NX)

    while(iterate){
        b0 <- b
        for (ii in 1:I){
            sc_ii <- score_items[ii,]
            b_ii <- c(0, b[ii,] )
            probs_ii <- probs_pcm_one_item(theta=theta, b_ii=b_ii)
            probs[,ii,] <- probs_ii
            r[ii,] <- colSums( weights * probs_ii * dat_resp[,ii] )[-1]
            for (kk1 in 1:K){
                for (kk2 in kk1:K){
                    rr[ii,kk1,kk2] <- sum( weights * probs_ii[,kk1+1] * probs_ii[,kk2+1] * dat_resp[,ii] )
                    if (kk1 < kk2){
                        rr[ii,kk2,kk1] <- rr[ii,kk1,kk2]
                    }
                }
            }

        }  # end ii
        A_Sq <- AA_bari <- A_bari <- matrix( 0, nrow=NX, ncol=I )
        for (kk in 1:K){
            A_bari <- A_bari + t( A[, kk, ] * r[, kk ] )
            AA_bari <- AA_bari + t( A[, kk, ]^2 * r[, kk ] )
        }
        for (kk1 in 1:K){
            for (kk2 in 1:K){
                A_Sq <- A_Sq + t( A[,kk1,] * A[,kk2,] * rr[, kk1, kk2 ] )
            }
        }
        expected <- rowSums(A_bari, na.rm=TRUE) # sum over items
        err <- rowSums(AA_bari - A_Sq, na.rm=TRUE)   #sum over the items
        err_inv <- abs(1/( abs(err) + eps ))
        scores <- ItemScore - expected
        incr <- - err_inv*scores
        incr <- immer_trim_increment(incr=incr, max_incr=max_incr)
        xsi <- xsi + incr
        max_incr <- abs(incr)

        b <- immer_jml_calc_item_intercepts(A=A, xsi=xsi, b_fixed=b_fixed)
        iter <- iter + 1
        b_change <- max( abs( b - b0) )
        if (iter >=maxiter_update){ iterate <- FALSE }
        if (b_change < conv_update){ iterate <- FALSE }
    }
    #--- output
    res <- list(b=b, xsi=xsi, xsi_der2=err, probs=probs )
    return(res)
}
