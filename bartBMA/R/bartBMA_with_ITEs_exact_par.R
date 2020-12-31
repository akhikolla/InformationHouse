#' @title Prediction intervals for bart-bma output obtained using linear algebra to obtain means and variances, and using bisection to find the quantiles of the mixture of t distributions.
#' 
#' @description This function produces prediction intervals for bart-bma output.
#' @param l_quant Lower quantile of credible intervals for the ITEs, CATT, CATNT.
#' @param u_quant Upper quantile of credible intervals for the ITEs, CATT, CATNT.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residuals.
#' @param num_cores Number of cores used in parallel.
#' @param root_alg_precision The algorithm should obtain approximate bounds that are within the distance root_alg_precision of the true quantile for the chosen average of models.
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
#' @param exact_residuals Binary variable. If equal to 1, then trees are added to sum-of-tree models within each round of the algorithm by detecting changepoints in the exact residuals. If equals zero, then changepoints are detected in residuals that are constructed from approximate predictions.
#' @param spike_tree If equal to 1, then the Spike-and-Tree prior will be used, otherwise the standard BART prior will be used. The number of splitting variables has a beta-binomial prior. The number of terminal nodes has a truncated Poisson prior, and then a uniform prior is placed on the set of valid constructions of trees given the splitting variables and number of terminal nodes.
#' @param s_t_hyperprior If equals 1 and spike_tree equals 1, then a beta distribution hyperprior is placed on the variable inclusion probabilities for the spike and tree prior. The hyperprior parameters are a_s_t and b_s_t.
#' @param p_s_t If spike_tree=1 and s_t_hyperprior=0, then p_s_t is the prior variable inclusion probability.
#' @param a_s_t If spike_tree=1 and s_t_hyperprior=1, then a_s_t is a parameter of a beta distribution hyperprior.
#' @param b_s_t If spike_tree=1 and s_t_hyperprior=1, then b_s_t is a parameter of a beta distribution hyperprior.
#' @param lambda_poisson This is a parameter for the Spike-and-Tree prior. It is the parameter for the (truncated and conditional on the number of splitting variables) Poisson prior on the number of terminal nodes.
#' @param less_greedy If equal to one, then a less greedy model search algorithm is used.
#' @export 
#' @return The output is a list of length 4:
#' \item{ITE_intervals}{A 3 by n matrix, where n is the number of observations. The first row gives the l_quant*100 quantiles of the individual treatment effects. The second row gives the medians of the ITEs. The third row gives the u_quant*100 quantiles of the ITEs.} 
#' \item{ITE_estimates}{An n by 1 matrix containing the Individual Treatment Effect estimates.} 
#' \item{CATE_estimate}{The Conditional Average Treatment Effect Estimates} 
#' \item{CATE_Interval}{A 3 by 1 matrix. The first element is the l_quant*100 quantile of the CATE distribution, the second element is the median of the CATE distribution, and the thied element is the u_quant*100 quantile of the CATE distribution.} 
#' @examples 
#' \dontrun{
#' #Example of BART-BMA for ITE estimation
#' #Applied to data simulations from Hahn et al. (2020, Bayesian Analysis) 
#' #"Bayesian Regression Tree Models for Causal Inference: Regularization, Confounding, 
#' # and Heterogeneous Effects
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
#' example_output <- bartBMA_with_ITEs_exact_par(l_quant = lbound, 
#'                                            u_quant= ubound, 
#'                                            x_covariates = x_covs_train, 
#'                                            z_train = ztrain, 
#'                                            y_train = ytrain)
#'                                            }

bartBMA_with_ITEs_exact_par <-function(l_quant,u_quant,newdata=NULL,update_resids=1,
                                       num_cores=1,root_alg_precision=0.00001,
                                       x_covariates,z_train ,y_train,
                                       a=3,nu=3,sigquant=0.9,c=1000,
                                       pen=12,num_cp=20,x.test=matrix(0.0,0,0),
                                       num_rounds=5,alpha=0.95,beta=2,split_rule_node=0,
                                       gridpoint=0,maxOWsize=100,num_splits=5,
                                       gridsize=10,zero_split=1,only_max_num_trees=1,
                                       min_num_obs_for_split=2, min_num_obs_after_split=2,
                                       exact_residuals=1,
                                       spike_tree=0, s_t_hyperprior=1, p_s_t=0.5, a_s_t=1,b_s_t=3,
                                       lambda_poisson=10,less_greedy=0){
  
  
  x_train <- cbind(z_train,x_covariates)
  
  object <- bartBMA(x.train= x_train,y.train =y_train ,
                    a=a,nu=nu,sigquant=sigquant,c=c,
                    pen=pen,num_cp=num_cp,x.test=x.test,
                    num_rounds=num_rounds,alpha=alpha,beta=beta,split_rule_node=split_rule_node,
                    gridpoint=gridpoint,maxOWsize=maxOWsize,num_splits=num_splits,
                    gridsize=gridsize,zero_split=zero_split,only_max_num_trees=only_max_num_trees,
                    min_num_obs_for_split=min_num_obs_for_split, min_num_obs_after_split=min_num_obs_after_split,
                    exact_residuals=exact_residuals,
                    spike_tree=spike_tree, s_t_hyperprior=s_t_hyperprior, p_s_t=p_s_t, a_s_t=a_s_t,b_s_t=b_s_t,
                    lambda_poisson=lambda_poisson,less_greedy=less_greedy)
  
  
  
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    ret<-pred_ints_ITE_outsamp_par(object$sumoftrees,object$obs_to_termNodesMatrix,
                                   object$response,object$bic,#min_possible, max_possible,
                                   object$nrowTrain,
                                   nrow(object$test_data),object$a,object$sigma,0,object$nu,
                                   object$lambda,#diff_inital_resids,
                                   object$test_data,l_quant,u_quant,num_cores,
                                   root_alg_precision,x_covariates
    )
    
    
    
  }else{if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    ret <- pred_ints_ITE_insamp_par(object$sumoftrees,
                                    object$obs_to_termNodesMatrix,
                                    object$response,object$bic,#min_possible, max_possible,
                                    object$nrowTrain,#nrow(object$test_data),
                                    object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,object$test_data,
                                    l_quant,u_quant,
                                    num_cores,
                                    root_alg_precision,x_covariates
    )
    
  }else{
    #if test data included in call to object
    ret<-pred_ints_ITE_outsamp_par(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,#min_possible, max_possible,
                                   object$nrowTrain,
                                   nrow(newdata), object$a,object$sigma,0,object$nu,
                                   object$lambda,#diff_inital_resids,
                                   newdata,l_quant,u_quant,num_cores,
                                   root_alg_precision,x_covariates
    )
    
  }}
  
  #PI<-apply(draws_from_mixture,2,function(x)quantile(x,probs=c(l_quant,0.5,u_quant)))
  
  
  
  #each row is a vector drawn from the mixture distribution
  
  
  names(ret)<-c("ITE_intervals",
                "ITE_estimates",
                "CATE_estimate",
                "CATE_Interval")
  
  
  ret
}