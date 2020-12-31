## File Name: immer_jml_facets_calc_probs.R
## File Version: 0.09

immer_jml_facets_calc_probs <- function(maxcat, N, K, design, des_names1,
    parm, parm_sign, is_step)
{
    K <- maxcat + 1
    KM <- matrix(0:maxcat, nrow=N, ncol=K, byrow=TRUE)
    probs <- matrix(0, nrow=N, ncol=K)
    for (uu in des_names1){
        des_uu <- design[,uu]
        probs <- probs + parm[[ uu ]][ des_uu ] * KM * parm_sign[[ uu ]]
    }
    # add step parameters
    if (is_step){
        uu <- "step"
        des_uu <- design[,uu]
        parm_uu <- parm[[uu]][ des_uu, ] * parm_sign[[ uu ]]
        for (mm in 1:maxcat){
            probs[, mm+1] <- probs[, mm+1] + rowSums( parm_uu[, 1:mm, drop=FALSE ] )
        }
    }
    # compute probabilities
    probs <- exp(probs)
    probs_sum <- rowSums(probs)
    probs <- probs / probs_sum

    #- compute case-wise mean and variance
    M <- rowSums( KM * probs )
    Var <- rowSums( KM^2 * probs ) - M^2

    #-- output
    res <- list(probs=probs, M=M, Var=Var, maxcat=maxcat, K=K)
    return(res)
}
