HWABO <- function(x,p=c(1/3,1/3,1/3),maxiter=50,tol=1e-10,verbose=TRUE) {
   fA <- x[1]
   fB <- x[2]
   nAB <- x[3]
   nOO <- x[4]
   n <- sum(x)
   nt <- 2*n
   
   di <- 10
   i <- 0
   It.hist <- NULL
   while(di > tol & i <= maxiter) {
    
     EuAA <- fA*p[1]/(p[1]+2*p[3])
     EuAO <- 2*fA*p[3]/(p[1]+2*p[3]) 
     EuBB <- fB*p[2]/(p[2]+2*p[3])
     EuBO <- 2*fB*p[3]/(p[2]+2*p[3])
     
     pt <- (2*EuAA + EuAO + nAB)/nt
     qt <- (2*EuBB + EuBO + nAB)/nt
     rt <- 1-(pt+qt)

     nA <- 2*EuAA + EuAO + nAB
     nB <- 2*EuBB + EuBO + nAB
     nhet <- nAB + EuAO + EuBO
     Kn <- lfactorial(n)
     Kd <- lfactorial(nOO) + lfactorial(fA) + 
           lfactorial(fB) + lfactorial(nAB) 
     K <- Kn-Kd
     ll <- K + 2*nOO*log(p[3]) + fA*log(p[1]*p[1]+2*p[1]*p[3]) + 
       fB*log(p[2]*p[2] + 2*p[2]*p[3]) + nAB*log(2*p[1]*p[2])
     res <- c(p,ll)
     pn <- c(pt,qt,rt)
     di <- pn-p
     di <- sum(di*di)  
     i <- i + 1
     p <- pn
     It.hist <- rbind(It.hist,res)
   }
   rownames(It.hist) <- 0:(nrow(It.hist)-1)
   colnames(It.hist) <- c("p","q","r","ll")
   if(verbose) {
     cat("Iteration history:\n")
     print(It.hist)
   }
   
   eAA <- n*p[1]^2
   eBB <- n*p[2]^2
   eOO <- n*p[3]^2
   eAB <- n*2*p[1]*p[2]
   eAO <- n*2*p[1]*p[3]
   eBO <- n*2*p[2]*p[3]
   expected <- c(eAA,eBB,eOO,eAB,eAO,eBO)
   names(expected) <- c("AA","BB","OO","AB","AO","OO")
   
   efA <- eAA+eAO
   efB <- eBB+eBO
   ee <- c(efA,efB,eAB,eOO)
   stats <- rbind(x, ee)
   rownames(stats) <- c("Observed","Expected")
   colnames(stats) <- c("fA","fB","nAB","nOO")
   if(verbose) print(stats)
   chicont <- ((x - ee)^2)/ee
   X2 <- sum(chicont)
   pval <- pchisq(X2,1,lower.tail = FALSE)
   if(verbose) {
     cat("X2 = ",X2,"p-value = ",pval,"\n")
   }
   names(pn) <- c("p","q","r")
   return(list(pn=pn,It.hist=It.hist,stats=stats,X2=X2,pval=pval))
}
