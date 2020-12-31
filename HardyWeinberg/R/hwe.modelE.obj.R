hwe.modelE.obj <-
function(params,x) {
  p.m  <- params[1] # male A allele frequency
  p.f  <- params[2] # female A allele frequency 
  f    <- params[3] # overall inbreeding coefficient
  n <- sum(x)
  pvec <- c(p.m^2 + (1-p.m)*p.m*f,
            2*p.m*(1-p.m)*(1-f),
            (1-p.m)^2 + p.m*(1-p.m)*f,
            p.f^2+f*(1-p.f)*p.f,
            2*p.f*(1-p.f)*(1-f),
            (1-p.f)^2 + p.f*(1-p.f)*f)
  ind <- !(x==0)
  lpvec <- log(pvec[ind])
  llik <- sum(x[ind]*lpvec)  
  -llik
}
