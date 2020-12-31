H0 <- function(x.f, x.m, alpha=c(0.5,0.5)) { 
    a <- alpha[1] + 2*x.f[1] + x.f[2] + x.m[1]
    b <- alpha[2] + 2*x.f[3] + x.f[2] + x.m[2]
    dens <- lfactorial(sum(x.f)) + lfactorial(sum(x.m)) + lgamma(sum(alpha)) + x.f[2]*log(2) + lgamma(a) + lgamma(b) - 
		sum(lfactorial(x.f)) - sum(lfactorial(x.m))  -  sum(lgamma(alpha)) - lgamma(a+b)
    dens <- exp(dens)
    return(dens)
}
