M11p <-
function(RR, alpha=c(0.5,0.5)) { 
 x.m <- RR[1:3]
 x.f <- RR[4:6]
    a <- alpha[1] + 2*x.f[1] + x.f[2] + 2*x.m[1] + x.m[2]
    b <- alpha[2] + 2*x.f[3] + x.f[2] + 2*x.m[3] + x.m[2]

    dens <- lfactorial(sum(x.f)) + lfactorial(sum(x.m)) - sum(lfactorial(x.f)) - sum(lfactorial(x.m)) + (x.f[2]+x.m[2])*log(2) +
		lgamma(sum(alpha)) - sum(lgamma(alpha)) +
		lgamma(a) + lgamma(b) - lgamma(a+b)

        dens <- exp(dens)

    return(dens)
}
