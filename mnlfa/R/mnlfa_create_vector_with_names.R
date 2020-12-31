## File Name: mnlfa_create_vector_with_names.R
## File Version: 0.02

mnlfa_create_vector_with_names <- function( vec, val=0 )
{
    NP <- length(vec)
    res <- rep(val, NP)
    names(res) <- vec
    return(res)
}
