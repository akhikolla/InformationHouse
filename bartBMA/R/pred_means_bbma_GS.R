#' @title Predictions for bart-bma output obtained from a Gibbs sampler
#' 
#' @description This function produces predictions from BART-BMA by post-hoc Gibbs-sampling from the full conditionals of the terminal node parameters and the variance of the error term. See Hernandez et al. (2018) Appendix D for details.
#' @param object bartBMA object obtained from function bartBMA
#' @param num_iter Total number of iterations of the Gibbs sampler (including burn-in).
#' @param burnin Number of burn-on iterations of the Gibbs sampler.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residua;s.
#' @export 
#' @return The output is a vector of predictions.
#' @examples 
#' set.seed(100)
#' #simulate some data
#' N <- 100
#' p<- 100
#' epsilon <- rnorm(N)
#' xcov <- matrix(runif(N*p), nrow=N)
#' y <- sin(pi*xcov[,1]*xcov[,2]) + 20*(xcov[,3]-0.5)^2+10*xcov[,4]+5*xcov[,5]+epsilon
#' epsilontest <- rnorm(N)
#' xcovtest <- matrix(runif(N*p), nrow=N)
#' ytest <- sin(pi*xcovtest[,1]*xcovtest[,2]) + 20*(xcovtest[,3]-0.5)^2+10*xcovtest[,4]+
#'   5*xcovtest[,5]+epsilontest
#' 
#' #Train the object 
#' bart_bma_example <- bartBMA(x.train = xcov,y.train=y,x.test=xcovtest,zero_split = 1, 
#'                             only_max_num_trees = 1,split_rule_node = 0)
#' #Obtain the prediction intervals
#' pred_means_bbma_GS(bart_bma_example,1000,100,newdata=NULL,update_resids=1)


pred_means_bbma_GS<-function(object,num_iter,burnin,newdata=NULL,update_resids=1){
  #object will be bartBMA object.
  if(update_resids==0){
    if(is.null(newdata) && length(object)==16){
      #if test data specified separately
      gs_chains<-gibbs_sampler_no_update_exp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                             nrow(object$test_data),object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals,object$test_data)
    }else if(is.null(newdata) && length(object)==14){
      #else return Pred Ints for training data
      gs_chains<-gibbs_sampler_no_update2_exp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                              object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals)
      
    }else{
      #if test data included in call to object
      gs_chains<-gibbs_sampler_no_update_exp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                             nrow(newdata), object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals,newdata)
    }
  }else{
    if(is.null(newdata) && length(object)==16){
      #if test data specified separately
      gs_chains<-gibbs_sampler_exp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                   nrow(object$test_data),object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals,object$test_data)
    }else if(is.null(newdata) && length(object)==14){
      #else return Pred Ints for training data
      gs_chains<-gibbs_sampler2_exp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                    object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals)
      
    }else{
      #if test data included in call to object
      gs_chains<-gibbs_sampler_exp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                   nrow(newdata), object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals,newdata)
    }
  }
  
  #y_posterior_sum_trees<-gs_chains[[4]]
  #y_orig_post_sum_trees<-gs_chains[[5]]
  #sigma_chains<-gs_chains[[3]]
  if(is.null(newdata) && length(object)==16){
    #y_posterior_sum_trees<-gs_chains[[4]]
    y_orig_post_sum_trees<-gs_chains[[2]] ##[[5]]
    sigma_chains<-gs_chains[[1]] #[[3]]
  }else if(is.null(newdata) && length(object)==14){
    #y_posterior_sum_trees<-gs_chains[[1]]
    y_orig_post_sum_trees<-gs_chains[[2]] #[[2]]
    sigma_chains<-gs_chains[[1]] #[[3]]
    
  }else{
    #y_posterior_sum_trees<-gs_chains[[4]]
    y_orig_post_sum_trees<-gs_chains[[2]] #[[5]]
    sigma_chains<-gs_chains[[1]] #[[3]]
  } 
  
  sum_of_tree_BIC<- -0.5*object$bic
  weights<-exp(sum_of_tree_BIC-(max(sum_of_tree_BIC)+log(sum(exp(sum_of_tree_BIC-max(sum_of_tree_BIC))))))
  #final_length<-num_iter-burnin
  num_its_to_sample<-round(weights*(num_iter-burnin))
  #final_sigma_chain<-numeric(0)
  
  #final_y_chain<-matrix(nrow=0,ncol=ncol(y_posterior_sum_trees[[1]]))
  final_yorig_chain<-matrix(nrow=0,ncol=ncol(y_orig_post_sum_trees[[1]]))
  
  for(i in 1:length(sigma_chains)){
    sample_its<-sample(burnin:num_iter,num_its_to_sample[i])
    #final_sigma_chain<-c(final_sigma_chain,sigma_chains[[i]][sample_its])
    #now do the same for predicted response updates
    #post_y_i<-y_posterior_sum_trees[[i]]
    post_yorig_i<-y_orig_post_sum_trees[[i]]
    #final_y_chain<-rbind(final_y_chain,post_y_i[sample_its,])
    final_yorig_chain<-rbind(final_yorig_chain,post_yorig_i[sample_its,])
    
  }
  meanpreds<-apply(final_yorig_chain,2,mean)
  
  ret<-meanpreds
  #length(ret)<-1
  #ret[[1]]<-meanpreds
  #ret[[2]] <- final_y_chain
  #ret[[3]] <- post_y_i
  
  class(ret)<-"pred_means_bbma.bartBMA"  
  ret
  
}