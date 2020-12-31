mloglik0 <-
function(jtms,TT=max(jtms),nu,g,Ig=function(x)sapply(x,function(y)integrate(g,0,y)$value)){
  n.jmp <- length(jtms)
  term1 <- -sum(log(nu(jtms)+
                    sapply(1:n.jmp,
                           function(i)sum(g(jtms[i]-head(jtms,i-1)))
                           )
                    ))
  ## term2 <- integrate(nu,0,TT)$value+sum(Ig(TT-jtms)-Ig(0))
  term2 <- integrate(nu,0,TT)$value+sum(Ig(TT-jtms))  
  term1+term2
}
