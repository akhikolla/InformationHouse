\name{variable_importance_measure}
\alias{variable_importance_measure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates variable Importance of a Regression Tree Model
}
\description{
Number of times a variable has been picked in the branch nodes of a (single) regression tree.
}
\usage{
variable_importance_measure(Model_VIM,NumVariable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model_VIM}{ Regression Tree model in which the variable importance is measured}
  \item{NumVariable}{ Number of variables in the training or testing matrix}
}
\details{
In time of calculating node cost of a tree of a random forest, a user defined number of variables are randomly picked. Among this, the best variable is chosen for the node using the node cost. While an important variable for a model will always come out as the best.
This function calculates the number of times a variable has been picked in the regression tree. It has been done by checking which variables are picked, how many times, in the branch nodes of the model.
}
\value{
Vector of size (1 x NumVariable), showing the number of repetition of variables (serially) in the branch nodes of the model.
}

\examples{
library(MultivariateRandomForest)
trainX=matrix(runif(50*100),50,100)
trainY=matrix(runif(50*5),50,5) 
n_tree=2
m_feature=5
min_leaf=5
testX=matrix(runif(10*100),10,100) 

theta <- function(trainX){trainX}
results <- bootstrap::bootstrap(1:nrow(trainX),n_tree,theta) 
b=results$thetastar

Variable_number=ncol(trainY)
if (Variable_number>1){
  Command=2
}else if(Variable_number==1){
  Command=1
} 
NumVariable=ncol(trainX)
NumRepeatation=matrix(rep(0,n_tree*NumVariable),nrow=n_tree)

for (i in 1:n_tree){
  Single_Model=NULL
  X=trainX[ b[ ,i],  ]
  Y=matrix(trainY[ b[ ,i],  ],ncol=Variable_number)
  Inv_Cov_Y = solve(cov(Y)) # calculate the V inverse
  if (Command==1){
    Inv_Cov_Y=matrix(rep(0,4),ncol=2)
  }
  Single_Model=build_single_tree(X, Y, m_feature, min_leaf,Inv_Cov_Y,Command)
  NumRepeatation[i,]=variable_importance_measure(Single_Model,NumVariable)
}
}

