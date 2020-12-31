
#Create a vector of factors for nb_cond conditions with nb_sample samples in each condition.
gen.cond=function(nb_cond=2,nb_sample=5){
x=NULL
for (i in 1:nb_cond){
  x=c(x,rep(i,nb_sample))
}
return(as.factor(x))
}

