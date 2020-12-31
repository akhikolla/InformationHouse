fobj <-
function(p,nAAf,nABf,nBBf,nA,nm,nAf,n) {
  vQ <- AuxQ(p,nm,nA,nAf)
  vQ2 <- vQ*vQ
  part1 <-  nAAf*(vQ2 - 
                    (nA - 3*p*nA + 2*(p^2)*nA - 2*p*n + 6*(p^2)*n - 4*(p^3)*n)*vQ + 
                    ((nA - 2*p*n)^2)*p*(p-1) - 
                    2*p*vQ2 + (p^2)*vQ2) # ok
  
  part2 <- -nABf*(p*vQ2 - (p^2)*vQ2 + 
                    (2*nA*(p^2) + nA -2*nA*p - 4*(p^3)*n - 2*p*n + 4*(p^2)*n)*vQ + 
                    ((nA - 2*p*n)^2)*p*(1-p)) # ok
  
  part3 <-  nBBf*((p^2)*vQ2 + (nA*p - 2*nA*p*p -2*(p^2)*n +4*(p^3)*n)*vQ + 
                    ((nA - 2*p*n)^2)*p*(p-1)) # ok
  y <- part1 + part2 + part3
  return(y)  
}
