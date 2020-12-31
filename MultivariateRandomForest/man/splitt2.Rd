\name{splitt2}
\alias{splitt2}
\title{Split of the Parent node}
\usage{
splitt2(X, Y, m_feature, Index, Inv_Cov_Y, Command, ff)
}
\arguments{
\item{X}{Input Training matrix of size M x N, M is the number of training samples and N is the number of features}

\item{Y}{Output Training response of size M x T, M is the number of samples and T is the number of output responses}

\item{m_feature}{Number of randomly selected features considered for a split in each regression tree node.}

\item{Index}{Index of training samples}

\item{Inv_Cov_Y}{Inverse of Covariance matrix of Output Response matrix for MRF (Input [0 0; 0 0] for RF)}

\item{Command}{1 for univariate Regression Tree (corresponding to RF) and 2 for Multivariate Regression Tree (corresponding to MRF)}

\item{ff}{Vector of m_feature from all features of X. This varies with each split}
}
\value{
List with the following components:
\item{index_left}{Index of the samples that are in the left node after splitting}
\item{index_right}{Index of the samples that are in the right node after splitting}
\item{which_feature}{The number of the feature that produces the minimum splitting cost}
\item{threshold_feature}{The threshold value for the node split. 
A feature value less than or equal to the threshold will go to the left node and it will go to the right node otherwise.}
}
\description{
Split of the training samples of the parent node into the child nodes based on the feature and threshold that produces the minimum cost
}
\details{
At each node of a regression a tree, a fixed number of features (m_feature) are selected randomly to be 
considered for generating the split. Node cost for all selected features along with possible n-1 thresholds for 
n samples are considered to select the feature and threshold with minimum cost.
}
\examples{
library(MultivariateRandomForest)
X=matrix(runif(20*100),20,100)
Y=matrix(runif(20*3),20,3)
m_feature=5
Index=1:20
Inv_Cov_Y=solve(cov(Y))
ff2 = ncol(X) # number of features
ff =sort(sample(ff2, m_feature)) 
Command=2#MRF, as number of output feature is greater than 1
Split_criteria=splitt2(X,Y,m_feature,Index,Inv_Cov_Y,Command,ff) 
}

