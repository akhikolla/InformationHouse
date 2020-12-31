#' @title Prediction intervals for bart-bma output obtained using linear algebra to obtain means and variances, and using bisection to find the quantiles of the mixture of t distributions.
#' 
#' @description This function produces prediction intervals for bart-bma output.
#' @param object bartBMA object obtained from function bartBMA
#' @param l_quant Lower quantile of credible intervals for the ITEs, CATT, CATNT.
#' @param u_quant Upper quantile of credible intervals for the ITEs, CATT, CATNT.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param num_cores Number of cores used in parallel.
#' @param root_alg_precision The algorithm should obtain approximate bounds that are within the distance root_alg_precision of the true quantile for the chosen average of models.
#' @export 
#' @return The output is a list of length 2:
#' \item{PI}{A 3 by n matrix, where n is the number of observations. The first row gives the l_quant*100 quantiles. The second row gives the medians. The third row gives the u_quant*100 quantiles.} 
#' \item{meanpreds}{An n by 1 matrix containing the estimated means.} 
#' @examples
#' \dontrun{
#' #load the package
#' library(bartBMA)
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
#' #Obtain the prediction intervals
#' pred_ints_exact_par(bart_bma_example,0.025,0.975,newdata=NULL,num_cores=1)
#' }

pred_ints_exact_par <-function(object,#min_possible,max_possible,
                               l_quant,u_quant,newdata=NULL,
                               num_cores=1,
                               root_alg_precision=0.00001){
  #object will be bartBMA object.
  
  
  
  
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    ret<-pred_ints_exact_outsamp_par(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,#min_possible, max_possible,
                                     object$nrowTrain,
                                     nrow(object$test_data),object$a,object$sigma,0,object$nu,
                                     object$lambda,#diff_inital_resids,
                                     object$test_data,l_quant,u_quant,num_cores,root_alg_precision
    )
    
    
    
  }else{if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    stop("Code not yet written for insample")
    
    #ret<-pred_ints_lin_alg_insamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
    #                              object$a,object$sigma,0,object$nu,
    #                              object$lambda,l_quant,u_quant
    #)
    
  }else{
    #if test data included in call to object
    ret<-pred_ints_exact_outsamp_par(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,#min_possible, max_possible,
                                     object$nrowTrain,
                                     nrow(newdata), object$a,object$sigma,0,object$nu,
                                     object$lambda,#diff_inital_resids,
                                     newdata,l_quant,u_quant,num_cores,root_alg_precision
    )
    
  }}

  names(ret) <- c("PI", "meanpreds")
  class(ret)<-"pred_intervals.bartBMA"  
  ret
}