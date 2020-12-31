#' Model of a single tree of Random Forest or Multivariate Random Forest
#' 
#' Build a Univariate Regression Tree (for generation of Random Forest (RF) ) or Multivariate Regression Tree ( for generation of Multivariate Random Forest (MRF) ) using the training samples,
#'  which is used for the prediction of testing samples.  
#'  
#' @param X Input Feature matrix of M x N, M is the number of training samples and N is the number of input features
#' @param Y Output Feature matrix of M x T, M is the number of training samples and T is the number of ouput features
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node, which must be positive integer and less than N (number of input features)
#' @param min_leaf Minimum number of samples in the leaf node, which must be positive integer and less than or equal to M (number of training samples) 
#' @param Inv_Cov_Y Inverse of Covariance matrix of Output Response matrix for MRF(Input [0 0;0 0] for RF)
#' @param Command 1 for univariate Regression Tree (corresponding to RF) and 2 for Multivariate Regression Tree (corresponding to MRF)
#' @return 
#' Model of a single regression tree (Univariate or Multivariate Regression Tree). An example of the list of the non-leaf node:
#' \item{Flag for determining whether the node is leaf node or branch node. 0 means branch node and 1 means leaf node.}{1}
#' \item{Index of samples for the left node}{int [1:34] 1 2 4 5 ...}
#' \item{Index of samples for the right node}{int [1:16] 3 6 9 ...}
#' \item{Feature for split}{int 34}
#' \item{Threshold values for split, average them}{num [1:3] 0.655 0.526 0.785}
#' \item{List number for the left and right nodes}{num [1:2] 2 3}
#' An example of the list of the leaf node:
#' \item{Output responses}{num[1:4,1:5] 0.0724 0.1809 0.0699 ...}
#' @details
#' The regression tree structure is represented as a list of lists. For a non-leaf node, it contains the splitting criteria 
#' (feature for split and threshold) and for a leaf node, it contains the output responses for the samples contained in the leaf node. 
#' @export
#'
build_single_tree <- function(X, Y, m_feature, min_leaf,Inv_Cov_Y,Command){
  NN=round(nrow(X)/min_leaf)*50
  model=rep( list(NULL), NN )
  i=1
  Index=1:nrow(X)
  
  model=split_node(X,Y,m_feature,Index,i,model,min_leaf,Inv_Cov_Y,Command)
  return(model)
}