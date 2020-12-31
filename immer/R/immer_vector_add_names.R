## File Name: immer_vector_add_names.R
## File Version: 0.02

immer_vector_add_names <- function(vec, pre)
{
    n1 <- paste0(pre, 1:length(vec) )
    names(vec) <- n1
    return(vec)
}
