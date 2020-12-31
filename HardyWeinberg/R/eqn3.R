eqn3=function(p,i1,i2) {
  # p[1]: p
  # p[2]: pm
  # p[3]: pf
  ma <- min(p[1],1-p[1])
  z1=p[2] + ma/(1-ma)
  z2=p[3] + ma/(1-ma)
  z3=p[1] + p[2]/(1 - p[2]) # p overall satisfying its restrictions on f.m
  z4=p[1] - 1/(1 - p[2])
  z5=p[1] + p[3]/(1 - p[3]) # p overall satisfying its restrictions on f.f
  z6=p[1] - 1/(1 - p[3])
  return(c(z1,z2,z3,z4,z5,z6))
}
