HWLRtest <- function(x,y,scene.null="S1",scene.alt="S6",verbose=TRUE,tracing=0) {
  z <- x+y
  pam <- af(x)
  paf <- af(y)
  pa  <- af(z)
#  fm <- ifelse(mac(x)==0,1,HWf(x))
#  ff <- ifelse(mac(y)==0,1,HWf(y))
  fm <- HWf(x)
  ff <- HWf(y)
  f  <- HWf(z)
  out.null <- switch(scene.null,
         S1 = loglik.M11(pa,z), # EoAF and HWP (f's zero)
         S2 = loglik.M14(pa,z,f), # EoAF and EoIC
         S3 = loglik.M15(pa,f,x,y,tracing), # EAF only
         S4 = loglik.M21(x,y,pam,paf), # HWP for both sexes
         S5 = loglik.M24(x,y,tracing), #no EAF but with EIC
         S6 = loglik.M25(x,y,pam,paf,fm,ff))
  out.alt <- switch(scene.alt,
         S1 = loglik.M11(pa,z), # EoAF and HWP (f's zero)
         S2 = loglik.M14(pa,z,f), # EoAF and EoIC
         S3 = loglik.M15(pa,f,x,y,tracing), # EAF only
         S4 = loglik.M21(x,y,pam,paf), # HWP for both sexes
         S5 = loglik.M24(x,y,tracing), #no EAF but with EIC
         S6 = loglik.M25(x,y,pam,paf,fm,ff))
  loglik0 <- out.null[1]
  loglik1 <- out.alt[1]
  df <- out.alt[2] - out.null[2]
  if(df <= 0) stop("This test is not possible (no degrees of freedom)")
  loglambda <- loglik0 - loglik1
  G2 <- -2*loglambda
  pval <- pchisq(G2,df = df,lower.tail = FALSE)
  if(verbose) cat("G2 = ",G2,"df = ",df,"p-value = ",pval,"\n")
  return(list(G2=G2,df=df,pval=pval))
}
