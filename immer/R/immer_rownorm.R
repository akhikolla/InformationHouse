## File Name: immer_rownorm.R
## File Version: 0.01

immer_rownorm <- function(x)
{
    y <- x / rowSums(x)
    return(y)
}
