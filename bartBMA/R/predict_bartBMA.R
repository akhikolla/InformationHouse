#' @title Predictions for a new dataset using an existing bartbma object
#' 
#' @description This function produces predictions for a new dataset using a previously obtained bartBMA object.
#' @param object A bartBMA object obtained using the barBMA function.
#' @param newdata Covariate matrix for new dataset.
#' @export 
#' @return A vector of predictions for the new dataset.
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
#' predict_bartBMA(bart_bma_example,newdata=xcovtest)

predict_bartBMA<-function(object,newdata){
  #preds<-get_BART_BMA_test_predictions(newdata,object$bic,object$sumoftrees,object$y_minmax)
  #orig_preds<-preds[[1]]
  #class(orig_preds)<-"predict.bartBMA"
  #orig_preds
  
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    ret<-preds_bbma_lin_alg_outsamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,0, 0,object$nrowTrain,
                                    nrow(object$test_data),object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,
                                    object$test_data
    )
  }else{if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    ret<-preds_bbma_lin_alg_insamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,0, 0,object$nrowTrain,
                                   object$a,object$sigma,0,object$nu,
                                   object$lambda#diff_inital_resids
    )
    
  }else{
    #if test data included in call to object
    ret<-preds_bbma_lin_alg_outsamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,0, 0,object$nrowTrain,
                                    nrow(newdata), object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,
                                    newdata
    )
  }}
  
  
  class(ret)<-"predict.bartBMA"
  ret
  
}