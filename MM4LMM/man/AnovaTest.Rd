\name{AnovaTest}
\alias{AnovaTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Type I and Type III Tests for mixed models.
%%  ~~function to do ... ~~
}
\description{
This function computes Type I and Type III tests for each fixed effect of a model, as displayed by the \code{\link{MMEst}} function. Alternatively, a specific linear combination of the fixed parameters may be tested (at 0). 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AnovaTest(ResMMEst , TestedCombination=NULL , Type = "TypeIII" ,
    Cofactor = NULL , X = NULL , formula = NULL , VarList = NULL ,
    NbCores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ResMMEst}{
	A list as displayed by the \code{\link{MMEst}} function.
%%     ~~Describe \code{Y} here~~
}
  \item{TestedCombination}{
	A contrast matrix or a list of contrast matrices \eqn{C_m}. Each matrix corresponds to a (set of) linear combination to be (jointly) tested at 0.
}
  \item{Type}{
	"TypeI", "TypeIII" or "KR" (default is "TypeIII"). AnovaTest will compute tests of the given type for each fixed effect in the model. The option is ignored if a \code{TestedCombination} is provided. If \code{Type} is "KR" then AnovaTest will compute Type III test using Kenward Roger approximation, see Kenward and Roger (1997) for details.
}
  \item{Cofactor}{
	The incidence matrix corresponding to fixed effects common to all models used in \code{\link{MMEst}}. If \code{NULL}, a single intercept is used as cofactor. This entry is needed when \code{Type} is "KR".
%%     ~~Describe \code{Cofactor} here~~
}
  \item{X}{
	The incidence matrix or a list of incidence matrices corresponding to fixed effects specific to each model used in \code{\link{MMEst}} (default is \code{NULL}). This entry is needed when \code{Type} is "KR".
%%     ~~Describe \code{X} here~~
}
  \item{formula}{
	The formula object specifying the fixed effect part of all models separated by + operators used in \code{\link{MMEst}} (default is \code{NULL}). This entry is needed when \code{Type} is "KR".
%%     ~~Describe \code{formula} here~~
}

  \item{VarList}{
	The list of correlation matrices associated with random and residual effects used in \code{\link{MMEst}} (default is \code{NULL}). This entry is needed when \code{Type} is "KR".
%%     ~~Describe \code{VarList} here~~
}

  \item{NbCores}{
	The number of cores to be used. 
%%     ~~Describe \code{NbCores} here~~
}
}
\details{
	If no \code{TestedCombination} is provided, the function performs either Type I or Type III tests for each fixed effect in the model (default is Type III). If \code{TestedCombination} is provided, each linear combination in \code{TestedCombination} is tested at 0 using a Wald test. No check is performed regarding the estimability of the linear combination to be tested.  If the dimension of the combination does not match with the dimension of \code{ResMMEst}, a \code{NA} is returned. 

	
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The output of the function is a list with as many items as in the original input list \code{ResMMEst}. For each item of \code{ResMMEst}, a table is provided that contains the Wald test statistics, p-values and degrees of freedom for all tested hypotheses. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. Biometrics, 983-997.
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
  data(QTLDetectionExample)
  Pheno <- QTLDetectionExample$Phenotype
  Geno <- QTLDetectionExample$Genotype
  Kinship <- QTLDetectionExample$Kinship

  ##Build the VarList object
  VL <- list(Additive = Kinship , Error = diag(1,length(Pheno)))
  
  ##Perform inference
  Result <- MMEst(Y=Pheno , X = Geno , VarList = VL)

  ##Compute tests
  AOV <- AnovaTest(Result,Type="TypeI")
  
  ##Test specific contrast matrix
  TC = matrix(c(0,1),nrow=1)
  AOV2 <- AnovaTest(Result, TestedCombination = TC)
  
  AOV3 <- AnovaTest(Result, TestedCombination = TC , Type="KR" , X = Geno , VarList = VL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Test }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
