\docType{package}
\name{SelvarMix-package}
\alias{SelvarMix-package}
\title{Regularization for variable selection in model-based clustering and discriminant analysis}
\description{
  SelvarMix is a package where a regularization approach 
  of variable selection is considered in model-based 
  clustering and discriminant analysis frameworks. 
  First, this procedure consists of ranking the variables 
  with a lasso-like procedure. Second, the method of
  Maugis et al (2009, 2011) is adapted to define the 
  role of variables in the two frameworks. 
  SelvarMix provides a faster variable selection 
  algorithm than the backward stepwise or forward 
  stepwise algorithms of Maugis et al (2009), allowing us 
  to study high-dimensional datasets.
}
\details{
  \tabular{ll}{ Package: \tab SelvarMix\cr Type: \tab
  Package\cr Version: \tab 1.0\cr Date: \tab
  2014-04-03\cr License: \tab GPL-3 + file LICENSE\cr
  LazyLoad: \tab yes\cr }

  The general purpose of the package is to perform variable 
  selection in model-based clustering and
  discriminant analysis. It focus on model-based 
  clustering, where the clusters are assumed to 
  arise from Gaussian distributions. 
  The most achieved model in model-based clustering
  has been proposed by Maugis et al (2009). This 
  so-called \eqn{SRUW} modeling considers three roles
  of variables: one variable may belong to the relevant 
  clustering set \eqn{S}, the redundant variable set
  \eqn{U} or the independent variable set \eqn{W}. 
  Moreover, the redundant variables may be explained 
  by a subset \eqn{R} of the relevant variables \eqn{S}.
  In order to avoid the greedy algorithms when high-dimensional data 
  are studied, the SelvarMix procedure is proposed. 
  It proceeds in two steps: First, the variables
  are ranked using a lasso-like procedure analogous to the one of 
  Zhou et al (2009); second, the \eqn{SRUW} procedure is run 
  on this  ranked set of  variables. 
  

}
\examples{
\dontrun{
## wine data set 
## n = 178 observations, p = 27 variables 
data(wine)
## variable selection in model-based clustering  
set.seed(123)
obj <- SelvarClustLasso(x=wine[,1:27], nbcluster=1:5, nbcores=4)
summary(obj)
print(obj)

## variables selection in discriminant analysis
set.seed(123)
a <- seq(1, 178, 10)
b <- setdiff(1:178, a)
obj <- SelvarLearnLasso(x=wine[b,1:27], z=wine[b,28], xt=wine[a,1:27], zt=wine[a,28], nbcores=4)
summary(obj)
print(obj)
}
}
\author{
  Author: Mohammed Sedki, Gilles Celeux and Cathy Maugis-Rabusseau  
  
}
\references{
  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2009. 
  "Variable selection in model-based clustering: 
  A general variable role modeling". Computational 
  Statistics and Data Analysis, vol. 53/11, pp. 3872-3882.

  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2011. 
  "Variable selection in model-based discriminant analysis". 
  Journal of Multivariate Analysis, vol. 102, pp. 1374-1387.

  Zhou, H., Pan, W., and Shen, X., 2009. "Penalized model-based 
  clustering with unconstrained covariance matrices". 
  Electronic Journal of Statistics, vol. 3, pp.1473-1496.
  
  Sedki, M., Celeux, G., Maugis-Rabusseau, C., 2014. 
  "SelvarMix: A R package for variable selection in 
  model-based clustering and discriminant analysis with 
  a regularization approach". Inria Research Report 
  available at \url{http://hal.inria.fr/hal-01053784}

}
\keyword{package}

