eqnE <-
function(p,i1) {
  ma.m <- min(p[1],1-p[1])
  ma.f <- min(p[2],1-p[2])
  z1=p[3] + ma.m/(1-ma.m)   # f has to satisfy the -pm/(1-pm) limit for both sexes.
  z2=p[3] + ma.f/(1-ma.f)
  z3=p[1] + p[3]/(1 - p[3]) # p males satisfying its restrictions
  z4=p[1] - 1/(1 - p[3])
  z5=p[2] + p[3]/(1 - p[3]) # p females satisfying its restrictions
  z6=p[2] - 1/(1 - p[3])
  return(c(z1,z2,z3,z4,z5,z6))
}
