hwe.modelC.obj <- function(params,x) {
  p   <- params[1] # A allele frequency
  f.m <- params[2] # inbreeding coefficient males
  f.f <- params[3] # inbreeding coefficient females
  n <- sum(x)
  pvec <- c(p^2+f.m*(1-p)*p,
            2*p*(1-p)*(1-f.m),
            (1-p)^2 + p*(1-p)*f.m,
            p^2+f.f*(1-p)*p,
            2*p*(1-p)*(1-f.f),
            (1-p)^2 + p*(1-p)*f.f)
  ind <- !(x==0)
  lpvec <- log(pvec[ind])  
  llik <- sum(x[ind]*lpvec)
  -llik
}
