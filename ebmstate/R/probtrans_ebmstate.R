
#'Spline approximations of the cumulative hazard functions
#'
#'Creates a spline approximation for the vector of cumulative hazards of each transition.
#'
#'This function is used by the function \code{probtrans_by_convolution}. It is not meant to be called by the user.
#'
#'@param cumhaz An object of class \code{msfit}, created by 
#'\code{\link{msfit_generic}} or \code{\link{msfit}}.
#'@return A list of estimated cumulative hazard functions (one for each transition).
#'@seealso \code{\link{msfit_generic}}; \code{\link{msfit}}; \code{\link{probtrans_by_convolution}}.
#'@author Rui Costa

cumhaz_splines<-function(cumhaz){
  spline_list<-vector("list",length(unique(cumhaz$Haz$trans)))
  for(i in unique(cumhaz$Haz$trans)){
    cumhaz_subset<-cumhaz$Haz[cumhaz$Haz$trans==i,]
    spline_list[[i]]<-splinefun(cumhaz_subset[,"time"], cumhaz_subset[,"Haz"], method="monoH.FC")
  }
  spline_list
}

#'Find all possible paths until absorption from a given starting state
#'
#'\code{unique_paths} finds all possible sequences of states until absorption
#'when the process has a tree-like structure.
#'
#'This function is used by the function \code{\link{probtrans_by_convolution}}. 
#'It is not meant to be called by the user.
#'
#'@param from_state Initial state.
#'@param tmat A transition matrix describing the states and transitions in 
#' the multi-state model, as can be obtained by running
#' \code{\link{transMat}}. 
#' See argument \code{trans} in \code{\link{msprep}} (\code{mstate}
#' package) for more detailed information.
#'@return A matrix where each column is a sequence of states taken by the process until absorption. 
#'There are as many columns as the number of possible paths until absorption.
#'
#'@author Rui Costa
#'@seealso \code{\link{probtrans_by_convolution}};
#' \code{\link{transMat}}.

unique_paths<-function(from_state,tmat){
  M<-matrix(from_state)
  repeat{
    m<-NULL
    for(i in 1:ncol(M)){
      last_state_of_col_i<-M[dim(M)[1],i]
      if(is.na(last_state_of_col_i)){
        m<-cbind(m,c(M[,i],NA))
        next
      }
      target_states_of_last_state_of_col_i<-names(which(!is.na(tmat[last_state_of_col_i,])))
      if(length(target_states_of_last_state_of_col_i)>0){
        for(j in target_states_of_last_state_of_col_i){
          m<-cbind(m,c(M[,i],j))
        }
      }else{
        m<-cbind(m,c(M[,i],NA))
      }
      
    }
    M<-m
    if(sum(is.na(M[dim(M)[1],]))==dim(M)[2]){
      return(M)
    } 
  }
}

#'Find the unique possible path until an 
#'absorbing state
#'
#'From a \code{unique_paths} object that
#'shows all possible paths until 
#'absorption from an initial state, 
#'\code{successful_transitions} picks the path
#'that finishes in \code{to_state}, if there is one. 
#'The initial state is the one defined in the 
#'argument \code{from_state} 
#'to the function \code{unique_paths}. 
#'The process must have a tree-like structure.
#'
#'This function is used by \code{probtrans_by_convolution_Markov} and \code{probtrans_by_convolution_semiMarkov}.
#'It is not meant to be called by the user.
#'
#'@param unique_paths_object An object created by running 
#'\code{\link{unique_paths}}.
#'@param to_state An absorbing state.
#'@param tmat Transition matrix.
#'@return A vector with the unique sequence of states between two states. 
#'
#'@author Rui Costa
#'@seealso \code{\link{unique_paths}};
#' \code{\link{probtrans_by_convolution_Markov}};
#' \code{\link{probtrans_by_convolution_semiMarkov}}.

successful_transitions<-function(unique_paths_object,to_state,tmat){
  row_of_paths_object_with_to_state<-which(apply(unique_paths_object,1,function(x) sum(na.omit(x==to_state))>0))
  reduced_unique_paths_object<-matrix(unique(unique_paths_object[1:row_of_paths_object_with_to_state,],MARGIN = 2),ncol = ncol(unique_paths_object))
  sucessful_path_column<-which(apply(reduced_unique_paths_object,2,function(x) sum(x==to_state)>0))[1]
  successful_path<-reduced_unique_paths_object[,sucessful_path_column]
  if(length(successful_path)==1){
    return(NULL)
  }else{
    transitions<-NULL
    for(i in 1:(length(successful_path)-1)){
      transitions<-c(transitions,tmat[successful_path[i],successful_path[i+1]])
    }
    return(transitions)
  }
}

#'Compute the cumulative hazard of leaving a given state 
#'
#'@description
#'This function is not meant to be called by the user. It is
#'an internal function of \code{probtrans_by_convolution_Markov}
#'and \code{probtrans_by_convolution_semiMarkov}.
#'
#'\code{joint_cum_hazard_function} returns the cumulative
#'hazard of leaving state \code{i} to any state that can be
#'reached directly from \code{i}, at each of the time points in \code{t}.
#'There is no explicit argument \code{i}: this state 
#'is entirely defined by the transitions
#'that can occur when the patient is in it (and these 
#'transitions
#'are given in the argument \code{competing_transitions}).
#'
#'
#'@param t A vector of time points.
#'@param competing_transitions The transitions that can occur when the process
#'is in state \code{i}.
#'@param spline_list A list whose elements are spline functions 
#'approximating the cumulative hazard of making each possible transition in
#'the process. This is normally a list
#'object created by running \code{cumhaz_splines}.
#'@return A vector with the cumulative hazard of leaving a given state evaluated at given time points.
#'
#'@author Rui Costa
#'@seealso \code{\link{probtrans_by_convolution_Markov}};
#'\code{\link{probtrans_by_convolution_semiMarkov}};
#' \code{\link{cumhaz_splines}}.

joint_cum_hazard_function<-function(t,competing_transitions,spline_list){
  if(length(competing_transitions)>0){
    return(sum(sapply(competing_transitions,function(x) spline_list[[x]](t))))
  }else{
    return(rep(0,length(t)))
  }
}

#'Compute subject-specific transition probabilities
#'using convolution.
#'@details
#'The Markov model is a non-homogeneous Markov model 
#'in which the transition hazard rates depend only on
#'time since the initiating event. The semi-Markov model
#' has a single time scale: the sojourn time in the current 
#' state. This is sometimes called homogeneous semi-Markov
#' model. 
#' 
#' The algorithm behind \code{probtrans_ebmstate} is based 
#' on the convolution of density and survival functions and
#' is suitable for processes with a tree-like transition
#' structure only.
#'
#'@param initial_state The present function 
#'estimates transition probabilities from the state given
#'in this argument.
#'@param cumhaz An \code{msfit} object created by running
#'\code{mstate} or \code{mstate_generic}.
#'@param model Either 'Markov' or
#' 'semiMarkov'. See details.
#'@param max_time The maximum time for which transition probabilities
#'are estimated.
#'@param nr_steps The number of steps in the convolution algorithm
#' (larger increases precision but makes it slower)
#'
#'@return An object of class 'probtrans'. See the 'value' 
#'section in the help page of \code{mstate::probtrans}.
#'@example inst/examples/probtrans_ebmstate-example.R
#'@author Rui Costa & Moritz Gerstung
#'@seealso \code{\link{probtrans}};
#'
#'@export
probtrans_ebmstate<-function(initial_state,cumhaz,model,max_time=NULL,nr_steps=10000){
  if(is.null(max_time)) max_time<-max(cumhaz$Haz$time)
  probtrans_object<-lapply(initial_state, probtrans_by_convolution,tmat=cumhaz$trans,cumhaz=cumhaz,model=model,max_time=max_time,nr_steps=nr_steps)
  probtrans_object$trans<-cumhaz$trans
  probtrans_object$direction<-"forward"
  probtrans_object$predt<-0
  class(probtrans_object)<-"probtrans"
  return(probtrans_object) 
}

#'Compute all transition probabilities from a given state
#'using convolution
#'
#'@description
#'\code{probtrans_by_convolution} is an internal function of \code{probtrans_ebmstate} and
#'is not meant to be called directly by the user.
#'It is itself a wrapper for the functions \code{probtrans_by_convolution_Markov}
#'and \code{probtrans_by_convolution_semiMarkov}, which are the workhorses of the 
#'convolution algorithm.
#'
#'@details For more information on the arguments of this function 
#'see \code{\link{probtrans_ebmstate}}.
#'
#'@param tmat A transition matrix extracted from the \code{cumhaz} argument to 
#'\code{probtrans_ebmstate}. 
#'
#'@param cumhaz \code{msfit} object (argument passed on from \code{probtrans_ebmstate}).
#'@param from_state Initial state (argument passed on from \code{probtrans_ebmstate}).
#'@param model 'Markov' or 'semiMarkov' (argument passed on from \code{probtrans_ebmstate}).
#'@param max_time The maximum time for which transition probabilities
#'are estimated.
#'@param nr_steps The number of steps in the convolution algorithm
#' (larger increases precision but makes it slower)
#'
#'@author Rui Costa & Moritz Gerstung
#'@seealso \code{\link{probtrans_ebmstate}};\code{\link{probtrans_by_convolution_Markov}};
#'\code{\link{probtrans_by_convolution_semiMarkov}}.

probtrans_by_convolution<-function(tmat,cumhaz,from_state,model,max_time,nr_steps){
  spline_list<-cumhaz_splines(cumhaz)
  unique_paths_object<-unique_paths(from_state,tmat)
  all_states<-na.omit(unique(as.vector(unique_paths_object)))
  #all_target_states<-all_states[-which(all_states==from_state)]
  time<-seq(0,max_time,length.out=nr_steps)
  #time<-sort(c(0,sort(coxph.detail(fit)$time)-1,sort(coxph.detail(fit)$time)+1,maximum_time))
  if(model=="semiMarkov"){
    transprobs_for_all_states<-sapply(all_states, probtrans_by_convolution_semiMarkov,cumhaz=cumhaz,tmat=tmat,from_state=from_state,spline_list=spline_list,unique_paths_object=unique_paths_object,time=time)
  }else{
    transprobs_for_all_states<-sapply(all_states, probtrans_by_convolution_Markov,cumhaz=cumhaz,tmat=tmat,from_state=from_state,spline_list=spline_list,unique_paths_object=unique_paths_object,time=time)
  }
  non_reacheable_states<-rownames(tmat)[!rownames(tmat)%in%all_states]
  transprobs_for_all_states<-cbind(transprobs_for_all_states,matrix(0,nrow = length(time),ncol = length(non_reacheable_states),dimnames = list(NULL,non_reacheable_states)))
  
  #order_of_states_according_to_tmat<-rownames(tmat)[rownames(tmat)%in%colnames(transprobs_for_all_states)]
  transprobs_for_all_states<-transprobs_for_all_states[,rownames(tmat)]
  
  as.data.frame(cbind(time,transprobs_for_all_states))
  
}

#'Compute transition probabilities under a non-homogeneous Markov model using
#'a convolution algorithm.
#'
#'@description
#'Compute transition probabilities for a given starting state and target state
#'under a non-homogeneous Markov model, using a convolution algorithm.
#'
#'\code{probtrans_by_convolution_Markov} is an internal function of
#'\code{probtrans_by_convolution} and is not meant to be called directly by the user.
#'
#'@param tmat Transition matrix.
#'@param cumhaz \code{msfit} object.
#'@param from_state Initial state.
#'@param to_state Target state.
#'@param spline_list A list whose elements are spline functions 
#'approximating the cumulative hazard of making each possible transition in
#'the process. This is normally a list
#'object created by running \code{cumhaz_splines}.
#'@param unique_paths_object An object created by running \code{unique_paths}.
#'@param time A vector of ordered time points.
#'
#'@seealso \code{\link{probtrans_ebmstate}};
#'\code{\link{probtrans_by_convolution_semiMarkov}}; 
#'\code{\link{probtrans_by_convolution}};
#'\code{\link{unique_paths}};
#'\code{\link{cumhaz_splines}}.
#'@author Rui Costa & Moritz Gerstung

probtrans_by_convolution_Markov<-function(tmat,cumhaz,from_state,to_state,spline_list,unique_paths_object,time){
  row_of_tmat_of_current_state<-which(rownames(tmat)==from_state)
  competing_transitions<-na.omit(tmat[row_of_tmat_of_current_state,])
  probtrans_vector_1<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
  if(from_state==to_state){
    return(probtrans_vector_1)
  }
  successful_transitions_object<-successful_transitions(unique_paths_object,to_state,tmat)
  #lagged_differences_vector<-diff(spline_list[[successful_transitions_object[1]]](time))
  #integrand_1<-survival_function*c(0,lagged_differences_vector)
  
  
  #successful_transitions_object<-successful_transitions(unique_paths_object,to_state)
  #row_of_tmat_of_current_state<-which(tmat==successful_transitions_object[1],arr.ind = TRUE)[1]
  #competing_transitions<-na.omit(tmat[row_of_tmat_of_current_state,])
  #probtrans_vector_1<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
  
  
  for(i in 1:length(successful_transitions_object)){
    lagged_differences_vector<-c(diff(spline_list[[successful_transitions_object[i]]](time)))
    
    column_of_tmat_with_successful_trans<-which(tmat==successful_transitions_object[i],arr.ind = TRUE)[2]
    name_of_next_state<-colnames(tmat)[column_of_tmat_with_successful_trans]
    if(length(na.omit(tmat[name_of_next_state,]))==0) {
      probtrans_vector_2<-rep(1,length(time))
    }else{
      competing_transitions<-na.omit(tmat[name_of_next_state,])
      probtrans_vector_2<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
    }
    
    probtrans_vector_1<-convolute_Markov(time,lagged_differences_vector,probtrans_vector_1,probtrans_vector_2)
    
  }
  
  return(probtrans_vector_1)
}

#'Compute transition probabilities under a semi-Markov model using
#'a convolution algorithm.
#'
#'@description
#'Compute transition probabilities for a given starting state and target state
#'under a semi-Markov model with a single time scale (sojourn time),
#'using a convolution algorithm.
#'
#'\code{probtrans_by_convolution_semiMarkov} is an internal function of
#'\code{probtrans_by_convolution} and is not meant to be called directly by the user.
#'
#'@param tmat Transition matrix.
#'@param cumhaz \code{msfit} object.
#'@param from_state Initial state.
#'@param to_state Target state.
#'@param spline_list A list whose elements are spline functions 
#'approximating the cumulative hazard of making each possible transition in
#'the process. This is normally a list
#'object created by running \code{cumhaz_splines}.
#'@param unique_paths_object An object created by running \code{unique_paths}.
#'@param time A vector of ordered time points.
#'
#'@seealso \code{\link{probtrans_ebmstate}};
#'\code{\link{probtrans_by_convolution_Markov}}; 
#'\code{\link{probtrans_by_convolution}};
#'\code{\link{unique_paths}};
#'\code{\link{cumhaz_splines}}.
#'@author Rui Costa & Moritz Gerstung

probtrans_by_convolution_semiMarkov<-function(tmat,cumhaz,from_state,to_state,spline_list,unique_paths_object,time){
  row_of_tmat_of_current_state<-which(rownames(tmat)==from_state)
  competing_transitions<-na.omit(tmat[row_of_tmat_of_current_state,])
  survival_function<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
  if(from_state==to_state){
    return(survival_function)
  }
  successful_transitions_object<-successful_transitions(unique_paths_object,to_state,tmat)
  lagged_differences_vector<-diff(spline_list[[successful_transitions_object[1]]](time))
  integrand_1<-survival_function*c(lagged_differences_vector,0)
  for(i in 1:length(successful_transitions_object)){
    column_of_tmat_with_successful_trans<-which(tmat==successful_transitions_object[i],arr.ind = TRUE)[2]
    name_of_next_state<-colnames(tmat)[column_of_tmat_with_successful_trans]
    if(length(na.omit(tmat[name_of_next_state,]))==0) {
      integrand_2<-rep(1,length(time))
    }else if(i==length(successful_transitions_object)){
      competing_transitions<-na.omit(tmat[name_of_next_state,])
      integrand_2<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
    }else{
      competing_transitions<-na.omit(tmat[name_of_next_state,])
      survival_function<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
      lagged_differences_vector<-diff(spline_list[[successful_transitions_object[i+1]]](time))
      integrand_2<-survival_function*c(lagged_differences_vector,0)
    }
    
    integrand_1<-convolute_semiMarkov(time,integrand_1,integrand_2)
    
  }
  
  return(integrand_1)
}
