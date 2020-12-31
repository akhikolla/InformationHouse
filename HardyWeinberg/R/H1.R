H1 <- function(x.f, x.m, alpha=c(0.333,0.333,0.333)) {
 matriu <- support()
 ff3 <- apply(matriu, 1, function(el){t1 <- el[1]; t2 <- el[2]; t3 <- 1-el[1]-el[2]; phi <- (t1+0.5*t2)
                                      aux3 <- dmultinom(x.f, size=sum(x.f), prob=c(t1,t2,t3),log=TRUE) +
                                              dmultinom(x.m, size=sum(x.m), prob=c(phi, 1-phi), log=TRUE) +
                                              ddir(x=c(t1,t2,t3),alpha,log=TRUE)
	 return(exp(aux3))
 	})
 r <- matriu[2,1]-matriu[1,1]
 dens <- sum(ff3)*r*r  
 return(dens)
}
