H2 <- function(x.f, x.m, alpha.f=c(0.5,0.5), alpha.m=c(0.5,0.5)) { 
    a <- alpha.f[1] + 2*x.f[1] + x.f[2] 
    b <- alpha.f[2] + 2*x.f[3] + x.f[2] 
    dens.f <- lfactorial(sum(x.f)) + lgamma(sum(alpha.f)) + x.f[2]*log(2) + lgamma(a) + lgamma(b) - 
		sum(lfactorial(x.f)) - sum(lgamma(alpha.f)) - lgamma(a+b)
    dens.m <- lgamma(sum(x.m)+1) +  lgamma(sum(alpha.m)) +   sum(lgamma(x.m+alpha.m)) -
		sum(lgamma(x.m+1))  -  sum(lgamma(alpha.m))  -   lgamma(sum(x.m+alpha.m))
    dens <- dens.f + dens.m
    dens <- exp(dens)
    return(dens)
}
