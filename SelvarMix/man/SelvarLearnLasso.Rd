\name{SelvarLearnLasso}
\alias{SelvarLearnLasso}
\title{
Regularization for variable selection in discriminant analysis 
}
\description{
This function implements the variable selection in discriminant 
analysis using a lasso ranking on the variables as described in 
Sedki et al (2014). The variable ranking step uses the penalized 
EM algorithm of Zhou et al (2009) (adapted in Sedki et al (2014) for the 
discriminant analysis settings).
A testing sample can be used to compute the averaged 
classification error rate.
}
\usage{
SelvarLearnLasso(x, z, lambda, rho, type, rank, hsize,  models, 
                 rmodel, imodel, xtest, ztest, nbcores)
}

\arguments{

  \item{x}{
  matrix  containing quantitative data. 
  Rows correspond to observations and
  columns correspond to variables
}

\item{z}{
  an integer vector or a factor corresponding to labels of data.
}

 \item{lambda}{
  numeric listing of tuning parameters for \eqn{\ell_1} mean penalty
}
  \item{rho}{
  numeric listing of tuning parameters for \eqn{\ell_1} precision matrix penalty
}

 \item{type}{character defining the type of ranking procedure, must be "lasso" or "likelihood". Default is "lasso"}


  \item{rank}{
  integer listing the rank of variables with (the length of this vector must be equal to the number of variables in the dataset)
}

\item{hsize}{
optional parameter make less strength  the forward and backward 
algorithms to select \eqn{S} and \eqn{W} sets
}

\item{models}{a  Rmixmod [\code{\linkS4class{Model}}] object
  defining the list of models to run. The models 
  Gaussian_pk_L_C, Gaussian_pk_Lk_C, Gaussian_pk_L_Ck, 
  and Gaussian_pk_Lk_Ck are called by default 
  (see mixmodGaussianModel() in Rmixmod package to
  specify other models)}

\item{rmodel}{list of character defining the covariance matrix form for
  the linear regression of \eqn{U} on the \eqn{R} set of variable. 
  Possible values: "LI" for spherical form, "LB" for diagonal form and  
  "LC" for general form. Possible values: "LI", "LB", "LC", c("LI", "LB")
  , c("LI", "LC"), c("LB", "LC") and c("LI", "LB", "LC"). 
  Default is c("LI", "LB", "LC")}

\item{imodel}{list of character defining the covariance matrix form for
  independent variables \eqn{W}. Possible values: 
  "LI" for spherical form and "LB" for diagonal form.
  Possible values: 
  "LI", "LB", c("LI", "LB"). Default is c("LI", LB")}


\item{xtest}{matrix  containing quantitative testing data. 
  Rows correspond to observations and
  columns correspond to variables 
}
\item{ztest}{
  an integer vector or a factor of size
  number of testing observations. Each cell corresponds to a
  cluster affectation
}

\item{nbcores}{
number of CPUs to be used when parallel computing is used (default is 2)
}
}

\value{ 
%\item{nbCluster}{number of clusters}
\item{S }{The selected set of relevant clustering variables}
\item{R }{The selected subset of regressors}
\item{U }{The selected set of redundant variables}
\item{W }{The selected set of independent variables}
\item{criterionValue}{The criterion value for the selected model}
%\item{nbcluster}{The selected number of clusters}
\item{model}{The selected covariance model} %% the selected gaussian mixture form
\item{rmodel}{The  selected covariance form for the regression}
\item{imodel}{The selected covariance form for the independent variables}
\item{parameters}{Rmixmod [\code{\linkS4class{Parameter}}] object containing all mixture parameters}
\item{regparameters}{Matrix containing all regression coefficients, each column is the regression coefficients of one 
redundant variable on the selected R set}
\item{proba}{Optional : matrix containing the conditional probabilities of belonging to each cluster for the testing observations}
\item{partition}{Optional: vector containing the cluster assignments of the testing observations according to the Maximum-a-Posteriori rule. When testing dataset is missed, we use the training dataset as testing one}
\item{error }{Optional : error rate done by the predicted partition (obtained using Maximum-A-Posteriori rule). When testing dataset is missed, we use the training dataset as testing one}
}


\author{
Mohammed Sedki <\url{mohammed.sedki@u-psud.fr}>
}
\references{
  Zhou, H., Pan, W., and Shen, X., 2009. "Penalized model-based 
  clustering with unconstrained covariance matrices". 
  Electronic Journal of Statistics, vol. 3, pp.1473-1496.
  
  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2009. 
  "Variable selection in model-based clustering: 
  A general variable role modeling". Computational 
  Statistics and Data Analysis, vol. 53/11, pp. 3872-3882.
 
 
  Sedki, M., Celeux, G., Maugis-Rabusseau, C., 2014. 
  "SelvarMix: A R package for variable selection in 
  model-based clustering and discriminant analysis with 
  a regularization approach". Inria Research Report 
  available at \url{http://hal.inria.fr/hal-01053784}
}

\keyword{discriminant analysis, variable selection, lasso ranking and graphical lasso}
\seealso{
\link{SelvarClustLasso}
\link{SortvarLearn}
\link{SortvarClust}
\link{wine}
}
\examples{
\dontrun{
## wine data set 
## n = 178 observations, p = 27 variables 
data(wine)
set.seed(123)
a <- seq(1, 178, 10)
b <- setdiff(1:178, a)
obj <- SelvarLearnLasso(x=wine[b,1:27], z=wine[b,28], xt=wine[a,1:27], zt=wine[a,28], nbcores=4)
summary(obj)
print(obj)
}
}
