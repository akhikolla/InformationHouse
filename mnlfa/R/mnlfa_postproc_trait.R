## File Name: mnlfa_postproc_trait.R
## File Version: 0.01

mnlfa_postproc_trait <- function(parm_trait)
{
    v1 <- parm_trait$mu
    if ( length(v1) > 0 ){
        trait <- data.frame(type="mu", par=names(v1), est=v1, estim=1)
    } else {
        trait <- data.frame(type="mu", par="mu0", est=0, estim=0)
    }
    v1 <- parm_trait$sigma
    dfr1 <- data.frame(type="sigma", par=names(v1), est=v1, estim=1)
    trait <- rbind( trait, dfr1)
    rownames(trait) <- NULL
    return(trait)
}
