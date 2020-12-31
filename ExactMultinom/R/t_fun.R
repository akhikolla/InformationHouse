# Functions calculating test statistic values
t_prob = function(x,p){
  n = sum(x)
  2*sum(n*p*log(p) - lgamma(n*p+1) - x*log(p) + lgamma(x+1))
}
t_chisq = function(x,p){
  n = sum(x)
  sum(x^2/n/p) - n 
} 
t_llr = function(x,p){
  n = sum(x)
  2*sum(ifelse(x == 0,0,x*log(x/n/p)))
}
