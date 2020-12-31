#' @title Variable inclusion probabilities as defined by Linero (2018)
#' 
#' @description This measure defines the posterior inclusion probability of a variable as the model-probability weighted sum of indicator variables for whether the variable was used in any splitting rules in any of the trees in the sum-of-tree model.
#' @param object A bartBMA object obtained using the barBMA function.
#' @export 
#' @return A vector of posterior inclusion probabilities. The variables are ordered in the same order that they occur in columns of the input covariate matrix used to obtain the input bartBMA object.
#' @examples 
#' #set the seed
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
#' #Obtain the variable importances
#' varIncProb(bart_bma_example)

varIncProb<-function(object){
  #object will be bartBMA object.
  imp_vars2=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees)
  res<-apply((imp_vars2[[3]]>0)*imp_vars2[[1]],2,sum)
  #create varImpPlot command
  class(res)<-"varIncProb.bartBMA"
  res
}