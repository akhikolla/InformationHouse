
#' Title
#'
#' @param null a vector of likelihood ratio under the null hypothesis. Generated from simulation
#' @param LR the likelihood of the data
#' @keywords internal
getp_aud_estimate_pi_first = function(null, LR){
  prop = mean(null == 0) # probability of being zero
  null2 = null[null > 0]
  # Match two moments
  t1 = mean(null2)
  t2 = mean(null2^2)
  d = 2*t1^2/(t2-t1^2)
  a = t1/d
  p =(1-prop)*pchisq(a*LR, df = d,  lower.tail = F)
  return(list(p= p, prob_0 = prop, a = a, d = d ))
}

