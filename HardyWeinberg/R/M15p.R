M15p <-
function(RR, alpha.f=c(0.333,0.333, 0.333), alpha.m=c(0.333, 0.666), r=0.1) { 


 x.m <- RR[1:3]
 x.f <- RR[4:6]

 AA <- expand.grid(seq(r/2,1-r/2,by=r), seq(r/2,1-r/2,by=r), seq(r/2,1-r/2,by=r))

 paa.f <- AA[,1]
 pab.f <- AA[,2]
 paa.m <- AA[,3]

 pp <- cbind(paa.f,pab.f,paa.m)

 matriu <- pp[paa.f<=(1 -r/2 - pab.f) & 
	    paa.m >= (2*paa.f + pab.f -1)  &  
          paa.m <= (paa.f+ pab.f/2)  , ]

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


 ff3 <- apply(matriu, 1, function(el,r){

		pAA.f <- el[1]
		pAB.f <- el[2]
		pBB.f <- 1-pAA.f - pAB.f
		pAA.m <- el[3]
		pAB.m <- 2*pAA.f + pAB.f - 2*pAA.m
		pBB.m <- 1 - pAA.m - pAB.m


  	 	aux3 <- dmultinom(x.f, size=sum(x.f), prob=c(pAA.f,pAB.f,pBB.f),log=T) + 
			  dmultinom(x.m, size=sum(x.m), prob=c(pAA.m,pAB.m,pBB.m),log=T) +
			  ddir(x=c(pAA.f,pAB.f,pBB.f),alpha.f,log=T) + 
			  ddir(x=c(pAA.m,1-pAA.m),alpha.m,log=T) 


		aux4 <- exp(aux3)/(pbeta(pAA.f+ pAB.f/2,alpha.m[1],alpha.m[2]) -
			   pbeta(2*pAA.f + pAB.f -1,alpha.m[1],alpha.m[2]))

	 	return(aux4)

 })

 dens <- sum(ff3)*r*r*r  

 return(dens)

}
