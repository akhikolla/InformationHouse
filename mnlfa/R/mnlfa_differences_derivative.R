## File Name: mnlfa_differences_derivative.R
## File Version: 0.01

mnlfa_differences_derivative <- function(ll0, ll1, ll2, h)
{
    D1 <- (ll1 - ll2 ) / (2*h)
    D2 <- (ll1 - 2*ll0 + ll2 ) / h^2
    D2_max <- max(abs(D2))*(1+1E-3)
    #-- output
    res <- list(D1=D1, D2=D2, D2_max=D2_max)
    return(res)
}
