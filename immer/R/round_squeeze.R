## File Name: round_squeeze.R
## File Version: 0.04


#***********************************************
# round and squeeze
round_squeeze <- function( x, digits=0, min=-Inf, max=Inf)
{
    x1 <- round( x, digits=digits )
    x1 <- ifelse( x1 < min, min, x1 )
    x1 <- ifelse( x1 > max, max, x1 )
    return(x1)
}
#***********************************************
