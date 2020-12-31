#' @title ITE Predictions (in-sample) using bartBMA and the method described by Hill (2011)
#' 
#' @description This function produces ITE Predictions (in-sample) using bartBMA and the method described by Hill (2011).
#' @param x_covariates Covaraite matrix for training bartBMA.
#' @param z_train treatment vector for traiing bartBMA.
#' @param y_train outcome vector for training bartBMA.
#' @param a This is a parameter that influences the variance of terminal node parameter values. Default value a=3.
#' @param nu This is a hyperparameter in the distribution of the variance of the error term. THe inverse of the variance is distributed as Gamma (nu/2, nu*lambda/2). Default value nu=3.
#' @param sigquant Calibration quantile for the inverse chi-squared prior on the variance of the error term.
#' @param c This determines the size of Occam's Window
#' @param pen This is a parameter used by the Pruned Exact Linear Time Algorithm when finding changepoints. Default value pen=12.
#' @param num_cp This is a number between 0 and 100 that determines the proportion of changepoints proposed by the changepoint detection algorithm to keep when growing trees. Default num_cp=20.
#' @param x.test Test data covariate matrix. Default x.test=matrix(0.0,0,0).
#' @param num_rounds Number of trees. (Maximum number of trees in a sum-of-tree model). Default num_rounds=5.
#' @param alpha Parameter in prior probability of tree node splitting. Default alpha=0.95
#' @param beta Parameter in prior probability of tree node splitting. Default beta=1
#' @param split_rule_node Binary variable. If equals 1, then find a new set of potential splitting points via a changepoint algorithm after adding each split to a tree. If equals zero, use the same set of potential split points for all splits in a tree. Default split_rule_node=0.
#' @param gridpoint Binary variable. If equals 1, then a grid search changepoint detection algorithm will be used. If equals 0, then the Pruned Exact Linear Time (PELT) changepoint detection algorithm will be used (Killick et al. 2012). Default gridpoint=0.
#' @param maxOWsize Maximum number of models to keep in Occam's window. Default maxOWsize=100.
#' @param num_splits Maximum number of splits in a tree
#' @param gridsize This integer determines the size of the grid across which to search if gridpoint=1 when finding changepoints for constructing trees.
#' @param zero_split Binary variable. If equals 1, then zero split trees can be included in a sum-of-trees model. If equals zero, then only trees with at least one split can be included in a sum-of-trees model.
#' @param only_max_num_trees Binary variable. If equals 1, then only sum-of-trees models containing the maximum number of trees, num_rounds, are selected. If equals 0, then sum-of-trees models containing less than num_rounds trees can be selected. The default is only_max_num_trees=1.
#' @param min_num_obs_for_split This integer determines the minimum number of observations in a (parent) tree node for the algorithm to consider potential splits of the node.
#' @param min_num_obs_after_split This integer determines the minimum number of observations in a child node resulting from a split in order for a split to occur. If the left or right chikd node has less than this number of observations, then the split can not occur.
#' @export 
#' @return A list of length 2. The first element is A vector of Individual Treatment Effect Estimates. The second element is a bartBMA object (i.e. the trained BART-BMA model).
#' @examples 
#' n <- 250
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n) 
#' x4 <- rbinom(n,1,0.5)
#' x5 <- as.factor(sample( LETTERS[1:3], n, replace=TRUE))
#' 
#' p= 0
#' xnoise = matrix(rnorm(n*p), nrow=n)
#' x5A <- ifelse(x5== 'A',1,0)
#' x5B <- ifelse(x5== 'B',1,0)
#' x5C <- ifelse(x5== 'C',1,0)
#' 
#' x_covs_train <- cbind(x1,x2,x3,x4,x5A,x5B,x5C,xnoise)
#' 
#' #Treatment effect
#' #tautrain <- 3
#' tautrain <- 1+2*x_covs_train[,2]*x_covs_train[,4]
#' 
#' #Prognostic function
#' mutrain <- 1 + 2*x_covs_train[,5] -1*x_covs_train[,6]-4*x_covs_train[,7] +
#' x_covs_train[,1]*x_covs_train[,3]
#' sd_mtrain <- sd(mutrain)
#' utrain <- runif(n)
#' #pitrain <- 0.8*pnorm((3*mutrain/sd_mtrain)-0.5*x_covs_train[,1])+0.05+utrain/10
#' pitrain <- 0.5
#' ztrain <- rbinom(n,1,pitrain)
#' ytrain <- mutrain + tautrain*ztrain
#' #pihattrain <- pbart(x_covs_train,ztrain )$prob.train.mean
#' 
#' #set lower and upper quantiles for intervals
#' lbound <- 0.025
#' ubound <- 0.975
#' 
#' example_output <- ITEs_bartBMA(x_covariates = x_covs_train, 
#'                                z_train = ztrain,
#'                                y_train = ytrain)

ITEs_bartBMA<-function(x_covariates,z_train ,y_train,
                       a=3,nu=3,sigquant=0.9,c=1000,
                       pen=12,num_cp=20,x.test=matrix(0.0,0,0),
                       num_rounds=5,alpha=0.95,beta=2,split_rule_node=0,
                       gridpoint=0,maxOWsize=100,num_splits=5,gridsize=10,zero_split=1,only_max_num_trees=1,
                       min_num_obs_for_split=2, min_num_obs_after_split=2){
  
  x_train <- cbind(z_train,x_covariates)
  
  trained_bart_BMA <- bartBMA(x.train= x_train,y.train =y_train ,
                              a=a,nu=nu,sigquant=sigquant,c=c,
                              pen=pen,num_cp=num_cp,x.test=x.test,
                              num_rounds=num_rounds,alpha=alpha,beta=beta,split_rule_node=split_rule_node,
                              gridpoint=gridpoint,maxOWsize=maxOWsize,num_splits=num_splits,gridsize=gridsize,zero_split=zero_split,only_max_num_trees=only_max_num_trees,
                              min_num_obs_for_split=min_num_obs_for_split, min_num_obs_after_split=min_num_obs_after_split)
  
  if(nrow(x.test)==0){
    all_treated_data <- cbind(rep(1,nrow(x_covariates)), x_covariates)
    all_control_data <- cbind(rep(0,nrow(x_covariates)), x_covariates)
    
    preds_treated<-get_BART_BMA_test_predictions(all_treated_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
    preds_control<-get_BART_BMA_test_predictions(all_control_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
  }else{
    all_treated_data <- cbind(rep(1,nrow(x.test)), x.test)
    all_control_data <- cbind(rep(0,nrow(x.test)), x.test)
    
    preds_treated<-get_BART_BMA_test_predictions(all_treated_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
    preds_control<-get_BART_BMA_test_predictions(all_control_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
    
  }
  
  ret <- list()
  ret$ITE_ests<-preds_treated[[1]]-preds_control[[1]]
  ret$bbma_object <- trained_bart_BMA
  ret
}