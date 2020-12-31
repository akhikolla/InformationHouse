fder <-
function(p,nAAf,nABf,nBBf,nA,nm,nAf,n) {
  vQ <- AuxQ(p,nm,nA,nAf)
  vQ2 <- vQ*vQ
  part1 <- nAAf*(
    4*nm*(2*p*nm - nA + nAf) - 
      ( (-3*nA + 4*p*nA -2*n + 12*n*p - 12*n*(p^2))*vQ
        + (nA - 3*p*nA + 2*(p^2)*nA - 2*p*n + 6*(p^2)*n - 4*(p^3)*n)*2*nm)
    + 2*(nA - 2*p*n)*(-2*n)*p*(p-1) + 
      ((nA - 2*p*n)^2)*(2*p-1) -
      (2*vQ2 + 8*p*nm*(2*p*nm -nA +nAf)) +
      2*p*vQ2 + (p^2)*4*nm*(2*p*nm - nA + nAf)) # ok
  part2 <- -nABf*(vQ2 + 4*p*nm*(2*p*nm - nA + nAf) - 
                    (2*p*vQ2 + 4*(p^2)*nm*(2*p*nm - nA + nAf)) +
                    (4*nA*p - 2*nA - 12*n*(p^2) - 2*n + 8*p*n)*vQ +
                    (2*nA*(p^2) + nA - 2*nA*p - 4*(p^3)*n - 2*p*n + 4*(p^2)*n)*2*nm +
                    2*(nA - 2*p*n)*(-2*n)*p*(1-p) + ((nA - 2*p*n)^2)*(1-2*p)) # ok
  part3 <- +nBBf*(2*p*vQ2 + p*p*4*nm*(2*p*nm - nA + nAf) + 
                    (nA - 4*nA*p - 4*p*n + 12*n*p*p)*vQ +
                    (nA*p - 2*nA*p*p - 2*n*p*p + 4*p*p*p*n)*2*nm +
                    2*(nA - 2*p*n)*(-2*n)*(p*(p-1)) + 
                    ((nA - 2*p*n)^2)*(2*p - 1)) # ok
  y <- part1 + part2 + part3
  return(y)  
}
