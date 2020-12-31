
#' @title Probit BART_BMA for classification of a binary variable
#' 
#' @description This is an implementation of Bayesian Additive Regression Trees (Chipman et al. 2018) using Bayesian Model Averaging (Hernandez et al. 2018). 
#' @param x.train Training data covariate matrix
#' @param y.train Training data outcome vector.
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
#' @param a_s_t If spike_tree=1 and s_t_hyperprior=1, then a_s_t is a parameter of a beta distribution hyperprior
#' @param b_s_t If spike_tree=1 and s_t_hyperprior=1, then b_s_t is a parameter of a beta distribution hyperprior
#' @param lambda_poisson This is a parameter for the Spike-and-Tree prior. It is the parameter for the (truncated and conditional on the number of splitting variables) Poisson prior on the number of terminal nodes.
#' @param less_greedy If equal to one, then a less greedy model search algorithm is used.
#' @param ... Further arguments.
#' @rdname probit_bartBMA
#' @export 
#' @return The following objects are returned by bartbma:
#' \item{fitted.values}{The vector of predictions of the outcome for all training observations.} 
#' \item{sumoftrees}{This is a list of lists of matrices. The outer list corresponds to a list of sum-of-tree models, and each element of the outer list is a list of matrices describing the structure of the trees within a sum-of-tree model. See details.} 
#' \item{obs_to_termNodesMatrix}{This is a list of lists of matrices. The outer list corresponds to a list of sum-of-tree models, and each element of the outer list is a list of matrices describing to which node each of the observations is allocated to at all depths of each tree within a sum-of-tree model. See details.} 
#' \item{bic}{This is a vector of BICs for each sum-of-tree model.} 
#' \item{test.preds}{A vector of test data predictions. This output only is given if there is test data in the input.} 
#' \item{sum_residuals}{CURRENTLY INCORRECT OUTPUT. A List (over sum-of-tree models) of lists (over single trees in a model) of vectors of partial residuals. Unless the maximum number of trees in a model is one, in which case the output is a list (over single tree models) of vectors of partial residuals, which are all equal to the outcome vector.} 
#' \item{numvars}{This is the total number of variables in the input training data matrix.} 
#' \item{call}{match.call returns a call in which all of the specified arguments are specified by their full names.} 
#' \item{y_minmax}{Range of the input training data outcome vector.}
#' \item{response}{Input taining data outcome vector.}
#' \item{nrowTrain}{number of observations in the input training data.}
#' \item{sigma}{sd(y.train)/(max(y.train)-min(y.train))}
#' \item{a}{input parameter}
#' \item{nu}{input parameter}
#' \item{lambda}{parameter determined by the inputs sigma, sigquant, and nu}
#' \item{fitted.probs}{In-sample fitted probabilities}
#' \item{fitted.classes}{In-sample fitted classes}
#' @useDynLib bartBMA, .registration = TRUE
#' @examples 
#' #Example from BART package (McCulloch et al. 2019)
#' set.seed(99)
#' n=100
#' x = sort(-2+4*runif(n))
#' X=matrix(x,ncol=1)
#' f = function(x) {return((1/2)*x^3)}
#' FL = function(x) {return(exp(x)/(1+exp(x)))}
#' pv = FL(f(x))
#' y = rbinom(n,1,pv)
#' probit_bartBMA(x.train = X,y.train = y)

probit_bartBMA<-function(x.train,...)UseMethod("probit_bartBMA")

#' @rdname probit_bartBMA
#' @export probit_bartBMA.default
#' @export
probit_bartBMA.default<-function(x.train,y.train,
                                 a=3,nu=3,sigquant=0.9,c=1000,
                                 pen=12,num_cp=20,x.test=matrix(0.0,0,0),
                                 num_rounds=5,alpha=0.95,beta=2,split_rule_node=0,
                                 gridpoint=0,maxOWsize=100,num_splits=5,gridsize=10,zero_split=1,only_max_num_trees=1,
                                 min_num_obs_for_split=2, min_num_obs_after_split=2,
                                 exact_residuals=1,
                                 spike_tree=0, s_t_hyperprior=1, p_s_t=0.5, a_s_t=1,b_s_t=3,
                                 lambda_poisson=10,less_greedy=0,...){
  
  
  num_obs=nrow(x.train)
  num_vars=ncol(x.train)
  
  if(is.factor(y.train)) {
    if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
    #binary = TRUE
    #  y.train = as.numeric(y.train)-1
    #stop("Response must be a numeric vector")
  } else {
    if(!((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0))) {
      stop("y.train is a numeric with number of unique elements != 2")
      #cat('NOTE: assumming numeric response is binary\n')
      #binary = TRUE
      #stop("Response must be a numeric vector")
    }
  }
  
  
  Zlatent.train <- ifelse(y.train==1,3.1,-3.1)
  binary=FALSE
  start_mean=0
  start_sd=1
  mu=0
  sigma_mu=0
  sigma=sd(Zlatent.train)/(max(Zlatent.train)-min(Zlatent.train))
  qchi = qchisq(1.0-sigquant,nu,1,0)
  lambda = (sigma*sigma*qchi)/nu
  if(is.factor(Zlatent.train)) {
    # if(length(levels(Zlatent.train)) != 2) stop("Zlatent.train is a factor with number of levels != 2")
    binary = TRUE
    #  Zlatent.train = as.numeric(Zlatent.train)-1
    stop("Response must be a numeric vector")
  } else {
    if((length(unique(Zlatent.train)) == 2) & (max(Zlatent.train) == 1) & (min(Zlatent.train) == 0)) {
      cat('NOTE: assumming numeric response is binary\n')
      binary = TRUE
      #stop("Response must be a numeric vector")
    }
  }
  if(gridsize<1) stop("gridsize must be a positive integer")
  
  if(is.vector(x.train) | is.factor(x.train)| is.data.frame(x.train)) x.train = as.matrix(x.train)
  if(is.vector(x.test) | is.factor(x.test)| is.data.frame(x.test)) x.test = as.matrix(x.test)
  
  if(is.matrix(x.train)) {
    if(nrow(x.test)>0) {
      if(!is.matrix(x.test)) stop('x.test must be a matrix')
    } 
  }
  # check input arguments:
  # if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
  # if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
  if((!is.matrix(x.train))) stop("argument x.train must be a double matrix")
  if((!is.matrix(x.test)) ) stop("argument x.test must be a double matrix")
  if(!binary) {
    if((!is.vector(Zlatent.train))) stop("argument Zlatent.train must be a double vector")
  }
  if(nrow(x.train) != length(Zlatent.train)) stop("number of rows in x.train must equal length of Zlatent.train")
  if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
  #if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
  #if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
  #if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
  if(c<1)stop("Value of Occam's Window has to be greater than 0."); 
  if(num_cp<0 || num_cp>100)stop("Value of num_cp should be a value between 1 and 100."); 
  
  bartBMA_call=BART_BMA_sumLikelihood(less_greedy,
                                      spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,
                                      num_obs,num_vars,lambda_poisson,
                                      x.train,Zlatent.train,start_mean,start_sd,a,mu,nu,lambda,c,sigma_mu,
                                      pen,num_cp,x.test,num_rounds,alpha,beta,split_rule_node,gridpoint,maxOWsize,
                                      num_splits,gridsize,zero_split,only_max_num_trees,
                                      min_num_obs_for_split, min_num_obs_after_split,
                                      exact_residuals)
  
  if(length(bartBMA_call)==6){
    #length of bartBMA_call is 6 if test data was included in the call
    names(bartBMA_call)<-c("fitted.values","sumoftrees","obs_to_termNodesMatrix","bic","test.preds","sum_residuals")
    bartBMA_call[[6]]<-bartBMA_call[[6]]#[[length(bartBMA_call[[6]])]]
    bartBMA_call$test_data<-x.test
  }else{
    names(bartBMA_call)<-c("fitted.values","sumoftrees","obs_to_termNodesMatrix","bic","sum_residuals")
    bartBMA_call[[5]]<-bartBMA_call[[5]]#[[length(bartBMA_call[[5]])]]
  }
  
  bartBMA_call$numvars<-ncol(x.train)
  bartBMA_call$call<-match.call()
  bartBMA_call[[2]]<-bartBMA_call[[2]]#[[length(bartBMA_call[[2]])]]
  bartBMA_call[[3]]<-bartBMA_call[[3]]#[[length(bartBMA_call[[3]])]]
  bartBMA_call$y_minmax<-range(Zlatent.train)
  bartBMA_call$response<-Zlatent.train
  bartBMA_call$nrowTrain<-nrow(x.train)
  bartBMA_call$sigma<-sigma
  bartBMA_call$a<-a
  bartBMA_call$nu<-nu
  bartBMA_call$lambda<-lambda
  bartBMA_call$fitted.probs<-pnorm(bartBMA_call$fitted.values)
  bartBMA_call$fitted.classes<-ifelse(bartBMA_call$fitted.probs<=0.5,0,1)
  
  class(bartBMA_call)<-"probit_bartBMA"
  bartBMA_call
}