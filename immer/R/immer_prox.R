## File Name: immer_prox.R
## File Version: 0.06


immer_prox <- function(dat)
{
    dat_resp <- 1 - is.na(dat)
    dat[ is.na(dat) ] <- 0
    I <- ncol(dat)
    # descriptives
    N <- colSums(dat_resp)
    M <- colSums( dat*dat_resp) / N
    SD <- sqrt( colSums( dat^2*dat_resp) / N - M^2 )
    # score matrix
    maxK <- apply( dat, 2, max, na.rm=TRUE )
    K <- max(maxK)
    scores <- matrix( 0, nrow=I, ncol=K+1)
    colnames(scores) <- paste0("Cat", 0:K )
    rownames(scores) <- colnames(dat)
    for (kk in 0:K ){
        scores[,kk+1] <- colSums( ( dat==kk ) * dat_resp )
    }
    #- item difficulties
    b_item <- rep(0, I)
    names(b_item) <- rownames(scores)
    for (ii in 1:I){
        if ( maxK[ii] %% 2==1 ){
            ind1 <- seq(1, maxK[ii], by=2)
            ind0 <- seq(0, maxK[ii], by=2)
            ratio <- log( sum(scores[ii, ind1 + 1]) / sum(scores[ii, ind0 + 1]) )
            b_item[ii] <- M[ii] - sqrt(  1 + SD[ii]^2 / 2.9 ) * ratio
        }
        if ( maxK[ii] %% 2==0 ){
            ind1 <- seq(1, maxK[ii]-1, by=2)
            ind0 <- seq(0, maxK[ii]-1, by=2)
            ratio <- log( sum(scores[ii, ind1 + 1]) / sum(scores[ii, ind0 + 1]) )
            h1 <- M[ii] - sqrt(  1 + SD[ii]^2 / 2.9 ) * ratio
            ind1 <- seq(2, maxK[ii], by=2)
            ind0 <- seq(1, maxK[ii], by=2)
            ratio <- log( sum(scores[ii, ind1 + 1]) / sum(scores[ii, ind0 + 1]) )
            h2 <- M[ii] - sqrt(  1 + SD[ii]^2 / 2.9 ) * ratio
            b_item[ii] <- ( h1 + h2 ) / 2
        }
    }
    #- item-category parameters
    b <- matrix(NA, nrow=I, ncol=K)
    for (ii in 1:I){
        b[ii,1:maxK[ii] ] <- b_item[ii] * seq(1,maxK[ii] )
    }
    #---- output
    res <- list(N=N, M=M, SD=SD, I=I, maxK=maxK, K=K, scores=scores, b_item=b_item, b=b)
    return(res)

}
