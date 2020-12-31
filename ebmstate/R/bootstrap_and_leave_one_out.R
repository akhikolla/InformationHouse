
#' Bootstrap confidence intervals for transition probabilities
#' 
#' Generates 95\% highest density bootstrap interval estimates for transition probabilities computed using \code{probtrans_ebmstate} (semi-Markov version).
#' 
#' @param coxrfx_fits_boot The list of CoxRFX objects obtained by running \code{boot_coxrfx}.
#' @param patient_data (Single) patient data in `long format`, possibly with `expanded` covariates
#' (as obtained by running \code{mstate::expand.covs}).
#' @param tmat Transition matrix for the multi-state model, as obtained by running \code{mstate::transMat}
#' @param initial_state The initial state for which transition probability estimates should be computed
#' @param max_time The maximum time for which estimates should be computed
#' @return Interval estimates for transition probabilities. 
#' @author Rui Costa
#' @seealso \code{\link{probtrans_ebmstate}}; \code{\link{boot_coxrfx}}; 
#' \code{\link[mstate:trans]{transMat}}; \code{\link[mstate]{expand.covs}}


boot_probtrans<-function(coxrfx_fits_boot,patient_data,tmat,initial_state,max_time){
  msfit_objects_boot<-vector("list",length(coxrfx_fits_boot))
  probtrans_objects_boot<-vector("list",length(coxrfx_fits_boot))
  for(i in 1:length(coxrfx_fits_boot)){
    print(i)
    covariate_df<-as.data.frame(coxrfx_fits_boot[[i]]$Z)
    covariate_df$strata<-coxrfx_fits_boot[[i]]$strata
    mstate_data_expanded.boot<-list()
    mstate_data_expanded.boot$time<-coxrfx_fits_boot[[i]]$surv[,1]
    mstate_data_expanded.boot$status<-coxrfx_fits_boot[[i]]$surv[,2]
    patient_data2<-patient_data[names(patient_data)%in%names(covariate_df)]
    patient_data2$strata<-patient_data$strata
    
    #environment(coxrfx_fits_boot[[1]]$formula)$covariate_df<-covariate_df
    
    msfit_objects_boot[[i]]<-msfit_generic(coxrfx_fits_boot[[i]],patient_data2,trans=tmat)
    probtrans_objects_boot[[i]]<-probtrans_ebmstate(initial_state,msfit_objects_boot[[i]],"semiMarkov")[[1]]
    probtrans_objects_boot[[i]]<-probtrans_objects_boot[[i]][sapply(seq(from=0,to=max_time,length.out = 400),function(x) which.min(abs(probtrans_objects_boot[[i]]$time-x))),]
    
  }

  probtrans_CIs<-lapply(colnames(tmat),CIs_for_target_state,probtrans_objects_boot=probtrans_objects_boot)
  names(probtrans_CIs)<-colnames(tmat)
  return(list(probtrans_CIs=probtrans_CIs,probtrans_objects_boot=probtrans_objects_boot, msfit_objects_boot=msfit_objects_boot))
}

#' Ancillary function to \code{boot_ebmstate}.
#' 
#' Extracts the bootstrap estimates of transition probabilities for
#' target state `tstate` from a list
#' with bootstrap estimates of transition probabilities into multiple states.
#' This function is not meant to be called by the user.
#' 
#' @param list_object A list in which each individual element is a single
#' bootstrap estimate of the probability of transition
#' into different states.
#' @param tstate The state whose bootstrap estimates of transition probabilities we wish to extract
#' from \code{list_object}. 
#' @return Bootstrap estimates of transition probabilities into target state `tstate`. 
#' @details This function is an ancillary function of \code{CIs_for_target_state}, which
#' in turn is an ancillary function of \code{boot_ebmstate}.
#' @author Rui Costa
#' @seealso \code{\link{CIs_for_target_state}}; \code{\link{boot_ebmstate}} 

extract_function<-function(list_object,tstate){
  as.vector(list_object[tstate])
}

#' Ancillary function of \code{boot_ebmstate}.
#' 
#' Computes 95\% highest density bootstrap confidence 
#' intervals for the transition probabilities into \code{target_state}, 
#' given a list object with boostrap estimates of transition probabilities into multiple states. This 
#' function is not meant to be called by the user.
#' 
#' @param target_state The target state for whose transition probabilties the confidence intervals
#' are computed.
#' @param probtrans_objects_boot A list containing bootstrap estimates of transition probabilities.
#' @return 95\% highest density bootstrap confidence intervals for the transition
#' probabilities into \code{target_state}. 
#' @details Uses function \code{extract_function}.
#' @author Rui Costa
#' @seealso \code{\link{boot_ebmstate}}; \code{\link{extract_function}}.

CIs_for_target_state<-function(target_state,probtrans_objects_boot){
  target_state_boot_samples<-as.data.frame(sapply(probtrans_objects_boot, extract_function,tstate=target_state))
  apply(target_state_boot_samples,1,hdi,credMass=0.95)
}

#' Ancillary function of \code{boot_ebmstate}.
#' 
#' Computes 95\% highest density, non-parametric bootstrap confidence 
#' intervals for the cumulative hazard rate functions, 
#' given a list of \code{msfit} objects with boostrap estimates of cumulative hazard rate functions
#' for multiple transitions. This 
#' function is not meant to be called by the user.
#' 
#' @param transition The transition for which transition confidence intervals
#' are computed.
#' @param  msfit_objects_boot List of \code{msfit} objects with boostrap estimates 
#' of cumulative hazard rate functions
#' for multiple transitions.
#' @return 95\% highest density, non-parametric bootstrap confidence intervals for the cumulative
#' hazard rate functions.
#' @author Rui Costa
#' @seealso \code{\link{boot_ebmstate}}.


cumhazCIs_for_target_transition<-function(transition,msfit_objects_boot){
  unique_time_points<-sort(unique(unlist(sapply(msfit_objects_boot,function(x) unique(x[[1]][,"time"])))))
  cumhaz_fun<-function(msfit_object_boot,unique_time_point,transition){
    msfit_for_target_trans<-msfit_object_boot[[1]][msfit_object_boot[[1]][,"trans"]==transition,]
    msfit_for_target_trans[which.max(msfit_for_target_trans[,"time"]>=unique_time_point),"Haz"]
  }
  obj<-sapply(msfit_objects_boot,function(x) sapply(unique_time_points,cumhaz_fun,msfit_object_boot=x,transition=transition))
  output<-apply(obj,1,hdi,credMass=0.95)
  colnames(output)<-unique_time_points
  output
}


#' Bootstrap confidence intervals for regression coefficients
#' 
#' This function computes 95\% highest density bootstrap confidence intervals (non-parametric) for the regression coefficients estimated by CoxRFX.
#' 
#' @param mstate_data_expanded Data in `long format`, possibly with `expanded` covariates (as obtained by running mstate::expand.covs).
#' @param which_group A character vector with the same meaning as the `groups` argument of the function \code{CoxRFX} but named (with the covariate names).
#' @param min_nr_samples The confidence interval of any coefficient is based on a number of bootstrap samples at least as high as this argument. See details.
#' @param output Determines the sort of output. See value.
#' @param ... Further arguments to the CoxRFX function.
#' @return For each regression coefficient, the confidence intervals and the number of bootstrap samples on which they are based, if the `output` argument is equal to `CIs`; if `output` is equal to `CIs_and_coxrfx_fits`, also the \code{CoxRFX} objects for each bootstrap sample.  
#' @details In a given bootstrap sample there might not be enough information to generate 
#' estimates for all coefficients. If a covariate has little or no variation in a given bootstrap sample, 
#' no estimate of its coefficient will be computed. The present function will
#' keep taking bootstrap samples until every coefficient has been estimated
#' at least \code{min_nr_samples} times.
#' @author Rui Costa

boot_coxrfx<-function(mstate_data_expanded,which_group,min_nr_samples=100,output="CIs",...){
  coxrfx_fits_boot<-vector("list")
  rownames(mstate_data_expanded)<-1:nrow(mstate_data_expanded)
  boot_matrix<-matrix(nrow=0,ncol = sum(!names(mstate_data_expanded)%in%c("id","from","to","trans","Tstart","Tstop","time","status","strata","type")),dimnames = list(NULL,names(mstate_data_expanded)[!names(mstate_data_expanded)%in%c("id","from","to","trans","Tstart","Tstop","time","status","strata","type")]))
  j<-1
  repeat{
    boot_samples_trans_1<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==1,]),replace = TRUE)
    boot_samples_trans_2<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==2,]),replace = TRUE)
    boot_samples_trans_3<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==3,]),replace = TRUE)
    boot_samples<-c(boot_samples_trans_1,boot_samples_trans_2,boot_samples_trans_3) 
    
    mstate_data_expanded.boot<-mstate_data_expanded[boot_samples,]
    
    ##exclude covariates without variance and binary covariates with less than 0.05 cases
    # vars_to_exclude<-vector("list",3)
    # for(i in 1:3){
    #   string_split<-strsplit(names(mstate_data_expanded.boot),"[.]")
    #   var_indices<-sapply(string_split,function(x) x[length(x)]==as.character(i))
    #   dummy_dataset<-mstate_data_expanded.boot[mstate_data_expanded.boot$trans==i,var_indices]
    #   which_have_variance<-apply(dummy_dataset, 2, function(x) var(x)>0)
    #   vars_to_exclude[[i]]<-names(dummy_dataset)[!which_have_variance]
    #   dummy_dataset<-dummy_dataset[which_have_variance]
    #   non_categorical_vars<-paste0(c("age_log","hb","anc_log","plt_log","bm_blasts_logit","ring_sideroblasts_logit","ipss","date"),paste0(".",as.character(i)))
    #   percentage_of_ones<-apply(dummy_dataset[!names(dummy_dataset)%in%non_categorical_vars], 2, function(x) sum(x)/length(x))
    #   which_less_than_five_percent<-which(percentage_of_ones<0.05)
    #   vars_to_exclude[[i]]<-c(vars_to_exclude[[i]],names(percentage_of_ones)[which_less_than_five_percent])
    # }
    # mstate_data_expanded.boot<-mstate_data_expanded.boot[!names(mstate_data_expanded.boot)%in%unlist(vars_to_exclude)]
    # 
    
    covariate_df<-mstate_data_expanded.boot[!names(mstate_data_expanded.boot)%in%c("id","from","to","trans","Tstart","Tstop","time","status","type")]
    groups2<-which_group[names(covariate_df)[names(covariate_df)!="strata"]]
    
    coxrfx_fits_boot[[j]]<-CoxRFX(covariate_df,Surv(mstate_data_expanded.boot$time,mstate_data_expanded.boot$status),groups =groups2,... )
    
    if(coxrfx_fits_boot[[j]]$iter[1]!=as.list(coxrfx_fits_boot[[j]]$call)$max.iter & sum(is.na(coxrfx_fits_boot[[j]]$coefficients))==0){
      boot_matrix<-rbind(boot_matrix,rep(NA,ncol(boot_matrix)))
      boot_matrix[j,names(coxrfx_fits_boot[[j]]$coefficients)]<-coxrfx_fits_boot[[j]]$coefficients
      print(min(apply(boot_matrix, 2, function(x) sum(!is.na(x)))))
      j<-j+1
    } 
    if(min(apply(boot_matrix, 2, function(x) sum(!is.na(x))))>=min_nr_samples) break
    
  }
  
  CIs<-apply(boot_matrix,2,hdi,credMass=0.95)
  CIs<-rbind(CIs,apply(boot_matrix, 2, function(x) sum(!is.na(x))))
  dimnames(CIs)[[1]][3]<-"n_samples"
  if(output=="CIs_and_coxrfx_fits"){
    return(list(CIs=CIs,coxrfx_fits_boot=coxrfx_fits_boot))
  }else if(output=="CIs"){
    return(CIs)
  }
}

#' Bootstrap samples and bootstrap interval estimates
#' 
#' This function computes bootstrap samples of regression coefficients,
#' cumulative hazard functions, and transition probability functions.
#' 
#' @param mstate_data_expanded Data in `long format`, possibly with `expanded` covariates (as obtained by running mstate::expand.covs). See details.
#' @param which_group A character vector with the same meaning as the `groups` argument of the function \code{CoxRFX} but named (with the covariate names).
#' @param min_nr_samples The confidence interval of any coefficient is based on a number of bootstrap samples at least as high as this argument. See details.
#' @param patient_data The covariate data for which the estimates of cumulative hazards and transition probabilities are computed. 
#' Must contain: one row of data for each transition, all the covariate columns in the fitted model, and also the 'strata' column. 
#' @param initial_state The initial state for which transition probability estimates should be computed
#' @param tmat Transition matrix for the multi-state model, as obtained by running \code{mstate::transMat}
#' @param backup_file Path to file. Objects generated while the present function is running are stored in this file. 
#' This avoids losing all estimates if and when the algorithm breaks down. See argument \code{input_file}. 
#' @param input_file Path to \code{backup_file} (see argument \code{backup_file}). If this argument is given, all other arguments should be \code{NULL}.
#' @param time_model The model of time-dependency: either 'Markov' or 'semiMarkov'.
#' @param coxrfx_args Named list with arguments to the \code{CoxRFX} function other than \code{Z},\code{surv} and \code{groups}.
#' @param msfit_args Named list with arguments to the \code{msfit_generic.coxrfx} function other than \code{object},\code{newdata} and \code{trans}.
#' @param probtrans_args Named list with arguments to the \code{probtrans_ebmstate} function other than \code{initia_state},\code{cumhaz} and \code{model}.
#' @return A list with: 95\% bootstrap intervals for each regression coefficient and for transition probabilities; 
#' bootstrap samples of regression coefficients, cumulative hazards and transition probabilities.
#' @details In a given bootstrap sample there might not be enough information to generate 
#' estimates for all coefficients. If a covariate has little or no variation in a given bootstrap sample, 
#' no estimate of its coefficient will be computed. The present function will
#' keep taking bootstrap samples until every coefficient has been estimated
#' at least \code{min_nr_samples} times. After expansion, the original covariates should be
#' excluded from \code{mstate_data_expanded}.
#' @author Rui Costa
#' @export


boot_ebmstate<-function(mstate_data_expanded=NULL,which_group=NULL,min_nr_samples=NULL,
                      patient_data=NULL,initial_state=NULL,tmat=NULL,
                      backup_file=NULL,input_file=NULL,time_model=NULL,coxrfx_args=NULL,
                      msfit_args=NULL,probtrans_args=NULL){
  
  list2env(coxrfx_args,envir = environment())
  if(!is.null(input_file)){
    load(input_file)
  }else{
    coxrfx_fits_boot<-vector("list")
    msfit_objects_boot<-vector("list")
    probtrans_objects_boot<-vector("list")
    rownames(mstate_data_expanded)<-1:nrow(mstate_data_expanded)
    boot_matrix<-matrix(nrow=0,
                        ncol = sum(!names(mstate_data_expanded)%in%c("id","from","to","trans",
                                                                     "Tstart","Tstop","strata","time",
                                                                     "status","type")),
                        dimnames = list(NULL,names(mstate_data_expanded)[!names(mstate_data_expanded)%in%c("id","from","to","trans",
                                                                                                           "Tstart","Tstop","strata","time",
                                                                                                           "status","type")]))
    j<-1  
  }
  tol<-unlist(mget("tol",ifnotfound = list(function(tol) 0.001)))
  max.iter<- unlist(mget("max.iter",ifnotfound = list(function(max.iter) 50)))
  sigma0<- unlist(mget("sigma0",ifnotfound = list(function(sigma0) 0.1)))
  sigma.hat<- unlist(mget("sigma.hat",ifnotfound = list(function(sigma.hat) "df")))
  verbose<- unlist(mget("verbose",ifnotfound = list(function(verbose) FALSE)))
  repeat{
    boot_samples_trans_1<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==1,]),replace = TRUE)
    boot_samples_trans_2<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==2,]),replace = TRUE)
    boot_samples_trans_3<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==3,]),replace = TRUE)
    boot_samples<-c(boot_samples_trans_1,boot_samples_trans_2,boot_samples_trans_3) 
    
    mstate_data_expanded.boot<-mstate_data_expanded[boot_samples,]
    covariate_df<-mstate_data_expanded.boot[!names(mstate_data_expanded.boot)%in%c("id","from","to",
                                                                                   "Tstart","Tstop","time","status","type")]
    groups2<-which_group[names(covariate_df)[!names(covariate_df)%in%c("strata","trans")]]
    if(time_model=="semiMarkov"){
      surv_object<-Surv(mstate_data_expanded.boot$time,mstate_data_expanded.boot$status) 
    }else if(time_model=="Markov"){
      surv_object<-Surv(mstate_data_expanded.boot$Tstart,mstate_data_expanded.boot$Tstop,mstate_data_expanded.boot$status) 
    }
    which.mu<-unlist(mget("which.mu",ifnotfound = list(function(which.mu) unique(groups2))))
    coxrfx_fits_boot[[j]]<-CoxRFX(covariate_df,surv_object,groups2,which.mu =which.mu,
                                  tol = tol,
                                  max.iter = max.iter,
                                  sigma0 = sigma0,
                                  sigma.hat = sigma.hat,
                                  verbose = verbose,coxrfx_args)
    #coxrfx_fits_boot[[j]]<-do.call("CoxRFX",c(list(Z=covariate_df,surv=surv_object,groups=groups2),coxrfx_args))
    
    if(sum(is.na(coxrfx_fits_boot[[j]]$coefficients))==0){
      boot_matrix<-rbind(boot_matrix,rep(NA,ncol(boot_matrix)))
      boot_matrix[j,names(coxrfx_fits_boot[[j]]$coefficients)]<-coxrfx_fits_boot[[j]]$coefficients
      
      msfit_objects_boot[[j]]<-do.call("msfit_generic",c(list(object=coxrfx_fits_boot[[j]],newdata=patient_data,trans=tmat),msfit_args))
      probtrans_objects_boot[[j]]<-do.call("probtrans_ebmstate",c(list(initial_state=initial_state,cumhaz=msfit_objects_boot[[j]],model=time_model),probtrans_args))[[1]]
      print(min(apply(boot_matrix, 2, function(x) sum(!is.na(x)))))
      if(j %%5==0){
        save(coxrfx_fits_boot,probtrans_objects_boot,
             msfit_objects_boot,boot_matrix,j, file =backup_file)
      }
      j<-j+1
    } 
    if(min(apply(boot_matrix, 2, function(x) sum(!is.na(x))))>=min_nr_samples) break
    
  }
  
  CIs<-apply(boot_matrix,2,hdi,credMass=0.95)
  CIs<-rbind(CIs,apply(boot_matrix, 2, function(x) sum(!is.na(x))))
  dimnames(CIs)[[1]][3]<-"n_samples"
  
  probtrans_CIs<-lapply(colnames(tmat),CIs_for_target_state,
                        probtrans_objects_boot=probtrans_objects_boot)
  names(probtrans_CIs)<-colnames(tmat)
  
  cumhaz_CIs<-lapply(sort(unique(mstate_data_expanded$trans)),cumhazCIs_for_target_transition,
           msfit_objects_boot=msfit_objects_boot)
  
  
  
  return(list(coefficients_CIs=CIs,coxrfx_fits_boot=coxrfx_fits_boot,
              probtrans_CIs=probtrans_CIs,
              probtrans_objects_boot=probtrans_objects_boot, 
              msfit_objects_boot=msfit_objects_boot,
              patient_data=patient_data,cumhaz_CIs=cumhaz_CIs))
  
}

#' Leave-one-out estimation
#' 
#' This function computes leave-one-out estimation of regression coefficients,
#' cumulative hazard functions, and transition probability functions.
#' 
#' @param mstate_data Data in `long format`.
#' @param mstate_data_expanded Data in `long format`, possibly with `expanded` covariates (as obtained by running mstate::expand.covs).
#' @param which_group A character vector with the same meaning as the `groups` argument of the function \code{CoxRFX} but named (with the covariate names).
#' @param patient_IDs The IDs of the patients whose cumulative hazards and transition probabilities one wishes to estimate.
#' @param initial_state The initial state for which transition probability estimates should be computed
#' @param tmat Transition matrix for the multi-state model, as obtained by running \code{mstate::transMat}
#' @param backup_file Path to file. Objects generated while the present function is running are stored in this file. 
#' This avoids losing all estimates if and when the algorithm breaks down. See argument \code{input_file}. 
#' @param input_file Path to \code{backup_file} (see argument \code{backup_file}). If this argument is given, all other arguments should be \code{NULL}.
#' @param time_model The model of time-dependency: either 'Markov' or 'semiMarkov'.
#' @param coxrfx_args Named list with arguments to the \code{CoxRFX} function other than \code{Z},\code{surv} and \code{groups}.
#' @param msfit_args Named list with arguments to the \code{msfit_generic.coxrfx} function other than \code{object},\code{newdata} and \code{trans}.
#' @param probtrans_args Named list with arguments to the \code{probtrans_ebmstate} function other than \code{initia_state},\code{cumhaz} and \code{model}.
#' @return A list with: 95\% bootstrap intervals for each regression coefficient and for transition probabilities; 
#' bootstrap samples of regression coefficients, cumulative hazards and transition probabilities.
#' @details In a given bootstrap sample there might not be enough information to generate 
#' estimates for all coefficients. If a covariate has little or no variation in a given bootstrap sample, 
#' no estimate of its coefficient will be computed. The present function will
#' keep taking bootstrap samples until every coefficient has been estimated
#' at least \code{min_nr_samples} times.
#' @author Rui Costa
#' @export


loo_ebmstate<-function(mstate_data,mstate_data_expanded,which_group,
                     patient_IDs,initial_state,tmat,
                     backup_file=NULL,input_file=NULL,time_model=NULL,coxrfx_args=list(),
                     msfit_args=NULL,probtrans_args=NULL){
  list2env(coxrfx_args,envir = environment())
  if(!is.null(input_file)){
    load(input_file)
    indices<-j:length(patient_IDs)
  }else{
    coxrfx_fits_loo<-vector("list")
    msfit_objects_loo<-vector("list")
    probtrans_objects_loo<-vector("list")
    indices<-1:length(patient_IDs)
  }
  tol<-unlist(mget("tol",ifnotfound = list(function(tol) 0.001)))
  max.iter<- unlist(mget("max.iter",ifnotfound = list(function(max.iter) 50)))
  sigma0<- unlist(mget("sigma0",ifnotfound = list(function(sigma0) 0.1)))
  sigma.hat<- unlist(mget("sigma.hat",ifnotfound = list(function(sigma.hat) "df")))
  verbose<- unlist(mget("verbose",ifnotfound = list(function(verbose) FALSE)))
  
  for(j in indices){ 
    mstate_data_expanded_loo<-mstate_data_expanded[mstate_data_expanded$id!=patient_IDs[j],]
    covariate_df<-mstate_data_expanded_loo[!names(mstate_data_expanded_loo)%in%c("id","from","to","Tstart","Tstop","time","status","type")]
    groups2<-which_group[names(covariate_df)[!names(covariate_df)%in%c("strata","trans")]]
    if(time_model=="semiMarkov"){
      surv_object<-Surv(mstate_data_expanded_loo$time,mstate_data_expanded_loo$status) 
    }else if(time_model=="Markov"){
      surv_object<-Surv(mstate_data_expanded_loo$Tstart,mstate_data_expanded_loo$Tstop,mstate_data_expanded_loo$status) 
    }
    which.mu<-unlist(mget("which.mu",ifnotfound = list(function(which.mu) unique(groups2))))
    coxrfx_fits_loo[[j]]<-CoxRFX(covariate_df,surv_object,groups2,which.mu =which.mu,
                                 tol = tol,
                                 max.iter = max.iter,
                                 sigma0 = sigma0,
                                 sigma.hat = sigma.hat,
                                 verbose = verbose,coxrfx_args)
    if(sum(is.na(coxrfx_fits_loo[[j]]$coefficients))==0){
      patient_data<-mstate_data[mstate_data$id==patient_IDs[j],,drop=FALSE][1,][rep(1,length(unique(mstate_data$trans))),]
      patient_data$trans<-1:length(unique(mstate_data$trans))
      patient_data<-expand.covs(patient_data,
                                covs = names(patient_data)[!names(patient_data)%in%c("id","from","to","trans","strata","Tstart","Tstop","time","status","type")])
      patient_data<-patient_data[names(mstate_data_expanded)]
      patient_data$strata<-unique(mstate_data[c("trans","strata")])[,2]
      msfit_objects_loo[[j]]<-do.call("msfit_generic",c(list(object=coxrfx_fits_loo[[j]],newdata=patient_data,trans=tmat),msfit_args))
      probtrans_objects_loo[[j]]<-do.call("probtrans_ebmstate",c(list(initial_state=initial_state,cumhaz=msfit_objects_loo[[j]],model=time_model),probtrans_args))
      if(j %%5==0){
        save(patient_IDs,coxrfx_fits_loo,msfit_objects_loo,probtrans_objects_loo,j,file=backup_file)
      }
      print(j)
    } 
  }
  
  return(list(coxrfx_fits_loo=coxrfx_fits_loo,
              probtrans_objects_loo=probtrans_objects_loo, 
              msfit_objects_loo=msfit_objects_loo,
              patient_IDs=patient_IDs))
}


