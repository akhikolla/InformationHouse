## File Name: probs_pcm_one_item.R
## File Version: 0.06

probs_pcm_one_item <- function(theta, b_ii )
{
    K <- length(b_ii) - 1
    N <- length(theta)
    KM <- immer_matrix2( 0:K, nrow=N )
    bM <- immer_matrix2( b_ii, nrow=N )
    probs_ii <- exp( theta*KM - bM )
    probs_ii <- immer_rownorm( probs_ii )
    return(probs_ii)
}
