M24p <-
function(RR, alpha.f=c(0.333,0.333, 0.333), alpha.m=c(0.5, 0.5), r=0.1) { 
 x.m <- RR[1:3]
 x.f <- RR[4:6]

 AA <- expand.grid(seq(r/2,1-r/2,by=r), seq(r/2,1-r/2,by=r), seq(r/2,1-r/2,by=r))


 paa.f <- AA[,1]
 pab.f <- AA[,2]
 pa.m <- AA[,3]

 t1 <- (4*paa.f^2 + 4*paa.f*pab.f - 4*paa.f + pab.f^2)/(2*pab.f)
 t2 <- (-4*paa.f^2 - 4*paa.f*pab.f + 4*paa.f - pab.f^2 + 2*pab.f)/(2*pab.f)


 pp <- cbind(paa.f,pab.f,pa.m)

 matriu <- pp[paa.f<(1 -r/2 - pab.f) &
		   pa.m>(4*paa.f^2 + 4*paa.f*pab.f - 4*paa.f + pab.f^2)/(2*pab.f) &
		   pa.m<(-4*paa.f^2 - 4*paa.f*pab.f + 4*paa.f - pab.f^2 + 2*pab.f)/(2*pab.f) ,]



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

		pAA.f <- el[1]
		pAB.f <- el[2]
		pBB.f <- 1-pAA.f - pAB.f

		pA.m <- el[3]

		pAA.m <- pA.m^2 + pA.m*(1-pA.m)*(pAA.f-(pAA.f+0.5*pAB.f)^2)/((pAA.f+0.5*pAB.f)*(1-pAA.f-0.5*pAB.f))
		pAB.m <- 2*pA.m*(1-pA.m)*(1-(pAA.f-(pAA.f+0.5*pAB.f)^2)/((pAA.f+0.5*pAB.f)*(1-pAA.f-0.5*pAB.f)))
		pBB.m <- 1- pAA.m - pAB.m


  	 	aux3 <- dmultinom(x.f, size=sum(x.f), prob=c(pAA.f,pAB.f,pBB.f),log=T) + 
			  dmultinom(x.m, size=sum(x.m), prob=c(pAA.m,pAB.m,pBB.m),log=T) +
			  ddir(x=c(pAA.f,pAB.f,pBB.f),alpha.f,log=T) + 
			  ddir(x=c(pA.m,1-pA.m),alpha.m,log=T)


   		a <- ( 4*pAA.f^2 + 4*pAA.f*pAB.f - 4*pAA.f + pAB.f^2)/(2*pAB.f)
  		b <- (-4*pAA.f^2 - 4*pAA.f*pAB.f + 4*pAA.f - pAB.f^2 + 2*pAB.f)/(2*pAB.f)


		aux4 <- exp(aux3)/(pbeta(b,alpha.m[1],alpha.m[2]) - pbeta(a,alpha.m[1],alpha.m[2]))


	 	return(aux4)
 })

 dens <- sum(ff3)*r*r*r  

 return(dens)

}
