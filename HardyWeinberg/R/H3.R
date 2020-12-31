H3 <- function(x.f, x.m, alpha.f=c(0.333,0.333,0.333), alpha.m=c(0.5,0.5)) { 
    if (!identical(length(x.f), length(alpha.f))) 
        stop("x.f and alpha.f differ in length.")
    if (!identical(length(x.m), length(alpha.m))) 
        stop("x.m and alpha.m differ in length.")
    x <- as.vector(x.f)
    alpha <- as.vector(alpha.f)
    dens.f <- lgamma(sum(x)+1) +  lgamma(sum(alpha)) +   sum(lgamma(x+alpha)) -
		sum(lgamma(x+1))  -  sum(lgamma(alpha))  -   lgamma(sum(x+alpha))
    x <- as.vector(x.m)
    alpha <- as.vector(alpha.m)
    dens.m <- lgamma(sum(x)+1) +  lgamma(sum(alpha)) +   sum(lgamma(x+alpha)) -
		sum(lgamma(x+1))  -  sum(lgamma(alpha))  -   lgamma(sum(x+alpha))
    dens <- dens.f + dens.m
    dens <- exp(dens)
    return(dens)
}
