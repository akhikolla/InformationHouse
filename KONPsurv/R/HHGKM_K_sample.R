# Create permutation for group vector with K values
# 
# @param group A binary vector indicating group, contains values from 1 to K
# @param imputed_status_matrix A matrix contianing the status each obsevations will recieve for each group it belongs.
# @return A permuted vector of group such that there are at least two events in each group .
# 
# @details This is an inner function inside KONPSURV package, not for the users.
sample_group_K_sample <- function(group,imputed_status_matrix)
{
  p_group <- sample(group)
  n <- length(group)
  status <- imputed_status_matrix[cbind(1:n,p_group)] # This is the status vector for the specific permutation
  tf_vec <- sapply(unique(group),function(k){sum(status[p_group==k])<2}) # A vector incdicating whether one of the groups has less then two events
  while(any(tf_vec)) #making sure all groups has et least two events
  {
    p_group <- sample(group)
    status <- imputed_status_matrix[cbind(1:n,p_group)] # This is the status vector for the specific permutation
    tf_vec <- sapply(unique(group),function(k){sum(status[p_group==k])<2}) # A vector incdicating whether one of the groups has less then two events
  }
  return(p_group)
}

#' KONP tests are \eqn{K}-sample Omnibus Non-Proportional hazards tests for right-censored data.
#' 
#' @param time A vector of the observed follow-up times.
#' @param status A vector of event indicators, 0=right censored, 1= event at \code{time}.
#' @param group A vector denoting the group labels, must contain at least two different values.
#' @param n_perm The number of permutations.
#' @param n_impu The number of imputations, for each imputation n_perm permutations will be executed.
#' 
#' @return Three test statistics and their respective p-values are returned: \cr 
#' 
#' \code{pv_chisq} - returns the p-value based on the KONP test chi-square statistic. \cr 
#' \code{pv_lr} - returns the p-value based on the KONP test likelihood ratio statistic. \cr 
#' \code{pv_cauchy} - returns the p-value based on the KONP-based Cauchy-combination test statistic. \cr 
#' \code{chisq_test_stat} - returns the KONP test chi-squared test statistic. \cr 
#' \code{lr_test_stat} - returns the KONP test likelihood-ratio test statistic. \cr
#' \code{cauchy_test_stat} - returns the KONP-based Cauchy-combination test statistic.
#' @examples
#' # gastric cancer data
#' data(gastric)
#' 
#' konp_test(gastric$time, gastric$status, gastric$group, n_perm=10^3) 
#' 
#' @details The KONP tests are powerful non-parametric tests for comparing \eqn{K} (>=2) hazard functions based on right-censored data.
#'These tests are consistent against any differences between the hazard functions of the groups.
#'The KONP tests are often more powerful than other existing tests, especially under non-proportional hazard functions.
#'  
konp_test<-function(time,status,group,n_perm,n_impu = 1)
{
  
  #Input testing 
  
  if (length(unique(group)) == 2)
  {
    return(konp_2_sample(time,status,group,n_perm,n_impu )) #running the faster function when K=2
  }
  
  if (all(unique(status)!= c(0,1)) & all(unique(status)!= c(1,0)) & all(unique(status)!= 1) & all(unique(status)!= 0))
  {
    stop ("ERROR - status vecotr must contain 0 or 1 only\n")
  }
  
  
  if (class(time) != "numeric" & class(time) != "integer")
  {
    stop ("ERROR - time class sould be numeric or integer\n")
  }
  
  if (length(time) != length(group) | length(time) != length(status))
  {
    stop ("ERROR - Vectors time, group and status must be in the same length\n")
  }
  
  if (sum(is.na(time))+ sum(is.na(status)) +sum(is.na(group))>0)
  {
    stop ("ERROR - time or status or group has NA's in the vector\n")
  }
  
  if (n_perm%%1 != 0 | n_perm<1)
  {
    stop ("ERROR - n_perm must be a natural number\n")
  }
  
  if (n_impu%%1 != 0 | n_impu<1)
  {
    stop ("ERROR - n_impu must be a natural number\n")
  }
  
  #Making the group vector be equal to 1 until K
  group_ex <- rep(NA,length(group))
  group_unq <- unique(group)
  K <- length(group_unq) #Nuber of different groups
  new_group <- seq(1,K,1)
  curr_group <- 1
  for (old_group in group_unq)
  {
    group_ex[group==old_group] <- curr_group
    curr_group <- curr_group+1
  }
  group <- group_ex
  
  
  #Making sure we have at least two events in each group
  for (gr in 1:K)
  {
    if (sum(status[gr==group])<2)
    {
      stop ("ERROR - Data must have at least two events in each groups in order to preform test\n")
    }
  }
  
  
  if (min(time)<=0)
  {
    stop ("ERROR - the time vector has negative or zero values\n")
  }
  
  n<-length(time)
  M <- Inf # represents a really large number
  
  test_prep <- sapply(1:K,function(k){ # this function will return all that is neccessery to run the test
    fit <- survival::survfit(survival::Surv(time[group==k], status[group==k]) ~ 1)
    sk <- c(1,fit$surv)
    time_k <- c(0,fit$time)
    n_k=sum(group==k)
    
    max_ev_k <- max(time[group==k & status==1])  #last event time in group k
    
    max_obs_k <- max(time[group==k])  #last event\censor time in group k
    
    tau_k <- ifelse(max_ev_k==max_obs_k,M,max_ev_k) #last time to estimate km in group 0
    
    return(list(sk=sk,time_k=time_k,n_k=n_k,tau_k=tau_k))
  })
  
  tau_k <- unlist(test_prep["tau_k",])
  tau <- sort(tau_k,decreasing = T)[2]
  n_k <- unlist(test_prep["n_k",])
  
  a <- hhgsurv_test_stat_K_sample(s_group = test_prep["sk",],time_group = test_prep["time_k",],n_vec = n_k,
                                  time=time,delta=status,trt=group,tau_k = tau_k ,tau = tau)
  
  chisq_test_stat <- a$chisq_stat
  
  lr_test_stat <- a$lr_stat
  
  #prepering for time only imputation
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  prob.t <- -diff(c(1,fit$surv)) 
  values.t <- fit$time
  if(sum(prob.t)<1)
  {
    prob.t <- c(prob.t,1-sum(prob.t))
    values.t <- c(values.t,max(values.t)+1) #each observation that will get max(values.t)+1 will be censored for sure
  }
  
  
  #prepearing for imputation for censoring (and some of the time)
  
  
  imputed_time_array <- array(data = NA,dim=c(n,K,n_impu)) #a time mat for each observation if the observation had changed a group
  imputed_status_array <- array(data = NA,dim=c(n,K,n_impu)) #a status mat for each observation if the observation had changed a group
  
  cen <- 1 - status
  
  for (k in 1:K)
  {
    fit <- survival::survfit(survival::Surv(time[group==k], cen[group==k]) ~ 1)
    prob.ck <- -diff(c(1,fit$surv)) # probabilities for treatment group censorship
    time_k <- fit$time
    if (sum(prob.ck)==0) #this claim will be correct if and only if there are no censorship in the group
    {
      prob.ck <- rep(1/2,2)
      time_k <- rep(Inf,2) #giving the censorship value a value that will be bigger then the time event for sure
    }else
    {
      if(sum(prob.ck)<1) # giving the correct value of being later then the last observation time
      {
        prob.ck <- c(prob.ck,1-sum(prob.ck))
        time_k <- c(time_k,max(time_k))
      }
    }
    for (subj in (1:n)[group!=k])
    {
      c <- sample(time_k,n_impu,prob=prob.ck,replace = T) + stats::rnorm(n = n_impu,sd = 10^-4)
      t <- rep(time[subj],n_impu)
      prob <- prob.t[values.t > time[subj]]
      values <- values.t[values.t > time[subj]]
      if ((length(prob)>0) & (sum(prob)>0) & status[subj]==0)
      { t <- sample(c(values,values), n_impu, prob = c(prob,prob),replace = T) + stats::rnorm(n = n_impu,sd = 10^-4)
      }
      imputed_status_array[subj,k,] <- ifelse(t<=c,1,0)
      imputed_time_array[subj,k,] <-  pmax(ifelse(t<=c,t,c),10^-10) #pmax is for not getting negative values
    }
    imputed_status_array[group==k,k,] <- status[group==k] #giving the original data when the observation didn't changed group
    imputed_time_array[group==k,k,] <- time[group==k]   #giving the original data when the observation didn't changed group
  }
  
  chisq_stat_perm <- c()
  lr_stat_perm <- c()
  
  pv_chisq_vec <- rep(NA,n_impu)
  pv_lr_vec <- rep(NA,n_impu)

  
  for (imp in 1:n_impu)
  {
    p_group_mat <- replicate(n_perm,sample_group_K_sample(group,imputed_status_array[,,imp])) #creating the permutation for n_perm permutations in a matrix
    
    perm <- get_perm_stats_K_sample(p_group_mat,imputed_time_array[,,imp],imputed_status_array[,,imp],n_perm = n_perm,n_vec = n_k)
    
    pv_chisq_vec[imp] <- (sum(chisq_test_stat<=perm$chisq_stat)+1)/(n_perm +1)
    pv_lr_vec[imp] <- (sum(lr_test_stat<=perm$lr_stat)+1)/(n_perm +1)
    
    chisq_stat_perm <- c(chisq_stat_perm,perm$chisq_stat)
    lr_stat_perm <- c(lr_stat_perm,perm$lr_stat)
  }
  
  
  pv_chisq <- (sum(chisq_test_stat<=chisq_stat_perm) + 1)/(n_impu*n_perm+1)
  pv_lr <- (sum(lr_test_stat<=lr_stat_perm) + 1)/(n_impu*n_perm+1)
  
  #Calculate Cauchy test statistic:
  #First calculate logrank P-value
  fit_lr <- survival::survdiff(survival::Surv(time, status) ~ group , rho=0)
  Pvalue_logrank <- 1 - stats::pchisq(fit_lr$chisq, 1)
  cau <- mean(c(tan((0.5-pv_chisq)*pi),
                tan((0.5-pv_lr)*pi), tan((0.5-Pvalue_logrank)*pi)))
  pv_cauchy <- 0.5-atan(cau)/pi
  
  return(list(pv_chisq=pv_chisq, pv_lr=pv_lr, pv_cauchy=pv_cauchy,
              chisq_test_stat=chisq_test_stat, lr_test_stat=lr_test_stat,
              cauchy_test_stat = cau))
}
