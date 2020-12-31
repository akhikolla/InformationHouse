\name{MMEst}
\alias{MMEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MM inference method for variance component mixed models 
%%  ~~function to do ... ~~
}
\description{
This is the main function of the \cite{MM4LMM} package. It performs inference  in  a variance component mixed model using a Min-Max algorithm. Inference in multiple models (e.g. for GWAS analysis) can also be performed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MMEst(Y, Cofactor = NULL, X = NULL, formula=NULL, VarList, ZList = NULL, Method = "Reml",
	Henderson=NULL, Init = NULL, CritVar = 0.001, CritLogLik = 0.001,
	MaxIter = 100, NbCores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	A vector of response values.
%%     ~~Describe \code{Y} here~~
}
  \item{Cofactor}{
	An incidence matrix corresponding to fixed effects common to all models to be adjusted. If \code{NULL}, a single intercept is used as cofactor.
%%     ~~Describe \code{Cofactor} here~~
}
  \item{X}{
	An incidence matrix or a list of incidence matrices corresponding to fixed effects specific to each model. If \code{X} is a matrix, one model per column will be fitted. If \code{X} is a list, one model per element of the list will be fitted (default is \code{NULL}).
%%     ~~Describe \code{X} here~~
}
  \item{formula}{
	A formula object specifying the fixed effect part of all models separated by + operators. To specify an interaction between \code{Cofactor} and \code{X} use the colnames of \code{X} when it is a list or use "Xeffect" when \code{X} is a matrix.
%%     ~~Describe \code{formula} here~~
}

  \item{VarList}{
	A list of correlation matrices associated with random and residual effects.
%%     ~~Describe \code{VarList} here~~
}
  \item{ZList}{
	A list of incidence matrices associated with random and residual effects (default is \code{NULL}).
%%     ~~Describe \code{ZList} here~~
}
  \item{Method}{
	The method used for inference. Available methods are "Reml" (Restricted Maximum Likelihood) and "ML" (Maximum Likelihood).
%%     ~~Describe \code{Method} here~~
}
  \item{Henderson}{
	If \code{TRUE} the Henderson trick is applied. If \code{FALSE} the Henderson trick is not applied. If \code{NULL} the algorithm chooses wether to use the trick or not.
}
  \item{Init}{
	A vector of initial values for variance parameters (default is \code{NULL})
%%     ~~Describe \code{Init} here~~
}
  \item{CritVar}{
	Value of the criterion for the variance components to stop iteration. (see Details)
%%     ~~Describe \code{CritVar} here~~
}
  \item{CritLogLik}{
	Value of the criterion for the log-likelihood to stop iteration. (see Details)
%%     ~~Describe \code{CritLogLik} here~~
}
  \item{MaxIter}{
	Maximum number of iterations per model.
%%     ~~Describe \code{MaxIter} here~~
}
  \item{NbCores}{
	Number of cores to be used. 
%%     ~~Describe \code{NbCores} here~~
}
}
\details{
	If \code{X} is \code{NULL}, the following model is fitted:
		\deqn{Y = X_C \beta_C + \sum_{k=1}^K Z_k u_k}
	with \eqn{X_C} the matrix provided in \code{Cofactor}, \eqn{\beta_C} the unknown fixed effects, \eqn{Z_k} the incidence matrix provided for the kth component of \code{ZList} and \eqn{u_k} the kth vector of random effects. If \code{ZList} is unspecified, all incidence matrices are assumed to be the Identity matrix. Random effects are assumed to follow a Gaussian distribution with mean 0 and covariance matrix \eqn{R_k \sigma_k^2}, where \eqn{R_k} is the kth correlation matrix provided in \code{VarList}.


	If \code{X} is not \code{NULL}, the following model is fitted for each i:
		\deqn{Y = X_C \beta_C + X_{[i]} \beta_{[i]} + \sum_{k=1}^K Z_k u_k}
	where \eqn{X_{[i]}} is the incidence matrix corresponding to the ith component (i.e. column if \eqn{X} is a matrix, element otherwise) of \eqn{X}, and \eqn{\beta_{[i]}} is the (unknow) fixed effect associated to \eqn{X_{[i]}}.
	
	All models are fitted using the MM algorithm. If \code{Henderson}=\code{TRUE}, at each step the quantities required for updating the variance components are computed using the Mixed Model Equation (MME) trick. See Johnson et al. (1995) for details. 

%%  ~~ If necessary, more details than the description above ~~
}
\value{
The result is a list where each element corresponds to a fitted model. Each element displays the following:
\item{Beta }{Estimated values of \eqn{\beta_C} and \eqn{\beta_{i}}}
\item{Sigma2 }{Estimated values of \eqn{\sigma_1^2,...,\sigma_K^2}}
\item{VarBeta }{Variance matrix of \code{Beta}}
\item{LogLik (\code{Method}) }{The value of the (restricted, if \code{Method} is "Reml") log-likelihood}
\item{NbIt }{The number of iterations required to reach the optimum}
\item{Method }{The method used for the inference}
\item{attr }{An integer vector with an entry for each element of \code{Beta} giving the term in \code{Factors} which gave rise to this element (for internal use in the function \code{\link{AnovaTest}})}
\item{Factors }{Names of each term in the formula}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Johnson, D. L., & Thompson, R. (1995). Restricted maximum likelihood estimation of variance components for univariate animal models using sparse matrix techniques and average information. Journal of dairy science, 78(2), 449-456.

Hunter, D. R., & Lange, K. (2004). A tutorial on MM algorithms. The American Statistician, 58(1), 30-37.

Zhou, H., Hu, L., Zhou, J., & Lange, K. (2015). MM algorithms for variance components models. arXiv preprint arXiv:1509.07426.
%% ~put references to the literature/web site here ~
}
\author{F. Laporte and T. Mary-Huard
%%  ~~who you are~~
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
require('MM4LMM')

#### Example 1: variance component analysis, 1 model
data(VarianceComponentExample)
DataHybrid <- VarianceComponentExample$Data
KinF <- VarianceComponentExample$KinshipF
KinD <- VarianceComponentExample$KinshipD

##Build incidence matrix for each random effect
Zf <- t(sapply(as.character(DataHybrid$CodeFlint), function(x)
  as.numeric(rownames(KinF)==x)))
Zd <- t(sapply(as.character(DataHybrid$CodeDent), function(x)
  as.numeric(rownames(KinD)==x)))

##Build the VarList and ZList objects
VL = list(Flint=KinF , Dent=KinD , Error = diag(1,nrow(DataHybrid)))
ZL <- list(Flint=Zf , Dent=Zd , Error = diag(1,nrow(DataHybrid)))

##Perform inference
#A first way to call MMEst
ResultVA <- MMEst(Y=DataHybrid$Trait  , Cofactor = matrix(DataHybrid$Trial)
                  , ZList = ZL  ,  VarList = VL)
length(ResultVA)
print(ResultVA)

#A second way to call MMEst (same result)
Formula <- as.formula('~ Trial')
ResultVA2 <- MMEst(Y=DataHybrid$Trait  , Cofactor = DataHybrid,
                   formula = Formula
                  , ZList = ZL  ,  VarList = VL)
length(ResultVA2)
print(ResultVA2)



#### Example 2: Marker Selection with interaction between Cofactor and X matrix
Formula <- as.formula('~ Trial+Xeffect+Xeffect:Trial')
ResultVA3 <- MMEst(Y=DataHybrid$Trait  , Cofactor = DataHybrid,
                  X = VarianceComponentExample$Markers,
                   formula = Formula
                  , ZList = ZL  ,  VarList = VL)
length(ResultVA3)
print(ResultVA3[[1]])


#### Example 3: QTL detection with two variance components
data(QTLDetectionExample)
Pheno <- QTLDetectionExample$Phenotype
Geno <- QTLDetectionExample$Genotype
Kinship <- QTLDetectionExample$Kinship

##Build the VarList object
VLgd <- list(Additive=Kinship , Error=diag(1,length(Pheno)))

##Perform inference
ResultGD <- MMEst(Y=Pheno , X=Geno
                  , VarList=VLgd , CritVar = 10e-5)

length(ResultGD)
print(ResultGD[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ LMM }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
