## File Name: immer_center.R
## File Version: 0.02

immer_center <- function(x, na.rm=TRUE)
{
    y <- x - mean(x, na.rm=na.rm)
    return(y)
}
