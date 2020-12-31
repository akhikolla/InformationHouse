

#' Title
#'
#' @param null a vector of likelihood ratio under the null hypothesis. Generated from simulation
#' @param LR the likelihood of the data
#' @keywords internal
getp_au1 = function(null, LR){
  t1 = mean(null)
  t2 = mean(null^2)
  B = max(0,1-3*t1^2/t2) # avoid negative probabiltiy
  A = t1/(1 - B)
  p =(1-B)*pchisq(A*LR, df = 1,  lower.tail = F)
  return(list(p= p, prob_0 = B, a = A))
}




