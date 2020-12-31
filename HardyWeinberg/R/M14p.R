M14p <-
function(RR, alpha=c(0.333,0.333, 0.333)) { 
 x.m <- RR[1:3]
 x.f <- RR[4:6]
    a <- alpha[1] + x.f[1] + x.m[1] 
    b <- alpha[2] + x.f[2] + x.m[2]
    c <- alpha[3] + x.f[3] + x.m[3]

    dens <- lfactorial(sum(x.f)) + lfactorial(sum(x.m)) - sum(lfactorial(x.f)) - sum(lfactorial(x.m)) +
		lgamma(sum(alpha)) -  sum(lgamma(alpha)) +
		lgamma(a) + lgamma(b) + lgamma(c) - lgamma(a+b+c)

        dens <- exp(dens)
    return(dens)

}
