#' Title
#'
#' @param null a vector of likelihood ratio under the null hypothesis. Generated from simulation
#' @param LR the likelihood of the data
#' @param qmax quantile of likelihood ratio vector. Most time put as 1, all the likelihood ratio vector under the null hypothesis will be used
#' @keywords internal
getp_Qreg_logP = function(null, LR, qmax = 1){
  # minimize lse between log theoratical p and empirical p
  # qmax =1, all LR0 are used.
  prop = mean(null == 0) # probability of being zero
  null2= null[null > 0]
  null2= null2[null2 >= stats::quantile(null2, 1-qmax)]
  f    = function(x){    # objective function wrt scale, dof
    scl   = x[1]; dof = x[2]
    the_p = 1 - stats::pchisq(null2/scl, df = dof)
    emp_p = 1 - (rank(null2)-0.5)/length(null2) # bigger null2, smaller p
    rr    = sum(log(the_p) -log(emp_p))^2
    return(rr)
  }
  test= stats::optim(c(0.5, 0.1) , f)
  scl = test$par[1]
  dof = test$par[2]
  final_p =(1-prop)*stats::pchisq(scl*LR, df = dof,  lower.tail = F)
  return(list(p= final_p, prob_0 = prop, a = scl, d = dof))
}

#' Title
#'
#' @param null a vector of likelihood ratio under the null hypothesis. Generated from simulation
#' @param LR the likelihood of the data
#' @param qmax quantile of likelihood ratio vector. Most time put as 1, all the likelihood ratio vector under the null hypothesis will be used
#' @keywords internal
getp_Qreg_P = function(null, LR, qmax = 1){
  # minimize lse between theoratical p and empirical p
  # qmax =1, all LR0 are used.
  prop = mean(null == 0) # probability of being zero
  null2= null[null > 0]
  null2= null2[null2 >= stats::quantile(null2, 1-qmax)]
  f    = function(x){    # objective function wrt scale, dof
    scl   = x[1]; dof = x[2]
    the_p = 1 - stats::pchisq(null2/scl, df = dof)
    emp_p = 1 - (rank(null2)-0.5)/length(null2) # bigger null2, smaller p
    rr    = sum(the_p -emp_p)^2
    return(rr)
  }
  test= stats::optim(c(0.5, 0.1) , f)
  scl = test$par[1]
  dof = test$par[2]
  final_p =(1-prop)*stats::pchisq(scl*LR, df = dof,  lower.tail = F)
  return(list(p= final_p, prob_0 = prop, a = scl, d = dof))
}

#' Title
#'
#' @param null  a vector of likelihood ratio under the null hypothesis. Generated from simulation
#' @param LR the likelihood of the data
#' @param qmax quantile of likelihood ratio vector. Most time put as 1, all the likelihood ratio vector under the null hypothesis will be used
#' @keywords internal
getp_Qreg_dir = function(null, LR, qmax = 1){
  # minimize lse between LR0 and theoratical quantile of aChi_d
  prop = mean(null == 0) # probability of being zero
  null2= null[null > 0]
  null2= null2[null2 >= stats::quantile(null2, 1-qmax)]
  f    = function(x){    # objective function wrt scale, dof
    scl   = x[1]; dof = x[2]
    emp_p = (rank(null2)-0.5)/length(null2)
    the_q = scl *stats::qchisq(p = emp_p , df = dof)
    rr    = sum(the_q -null2)^2
    return(rr)
  }
  test= stats::optim(c(0.5, 0.1) , f)
  scl = test$par[1]
  dof = test$par[2]
  final_p =(1-prop)*stats::pchisq(scl*LR, df = dof,  lower.tail = F)
  return(list(p= final_p, prob_0 = prop, a = scl, d = dof))
}
