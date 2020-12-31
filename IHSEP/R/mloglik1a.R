mloglik1a <-
function(jtms,TT=max(jtms),nu,g,Ig=function(x)sapply(x,function(y)integrate(g,0,y)$value),
         tol.abs=1e-12,tol.rel=1e-12,limit=1000){
  .External("mloglik1a",jtms,length(jtms),TT,nu,g,Ig,new.env(),
            tol.abs,tol.rel,limit,PACKAGE="IHSEP")
}
