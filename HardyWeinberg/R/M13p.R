M13p <-
function(RR, alpha=c(0.333,0.333, 0.333),r=0.1) { 
 x.m <- RR[1:3]
 x.f <- RR[4:6]
# r <- 0.005


 AA <- expand.grid(seq(r/2,1-r/2,by=r), seq(r/2,1-r/2,by=r))

 paa.m <- AA[,1]
 pab.m <- AA[,2]
 pp <- cbind(paa.m,pab.m)

 matriu <- pp[paa.m < 1 -r/2 -pab.m, ]



ddir <- function (x, alpha, log = T) {  # log = F
    x <- as.vector(x)
    alpha <- as.vector(alpha)
    if (!identical(length(x), length(alpha))) 
        stop("x and alpha differ in length.")
    dens <- lgamma(sum(alpha)) +  sum((alpha-1)*log(x)) -
		sum(lgamma(alpha)) 
    if (log == FALSE) 
        dens <- exp(dens)
    return(dens)
}


 ff3 <- apply(matriu, 1, function(el){
		pAA.m <- el[1]
		pAB.m <- el[2]
		pBB.m <- 1-pAA.m-pAB.m
		pAA.f <- ((2*pAA.m+pAB.m)/2)^2
		pAB.f <- 2*((2*pAA.m+pAB.m)/2)*(1-(2*pAA.m+pAB.m)/2)
		pBB.f <- 1-pAA.f - pAB.f
  	 	aux3 <- dmultinom(x.f, size=sum(x.f), prob=c(pAA.f,pAB.f,pBB.f),log=T) + 
			  dmultinom(x.m, size=sum(x.m), prob=c(pAA.m,pAB.m,pBB.m),log=T) +
			  ddir(x=c(pAA.m,pAB.m,pBB.m),alpha,log=T)
	 	return(exp(aux3))
 	})

 dens <- sum(ff3)*r*r  

 return(dens)

}
