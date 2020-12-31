\name{SortvarLearn}
\alias{SortvarLearn}
\title{
Variable ranking with LASSO in discriminant analysis 
}
\description{
This function implements variable ranking procedure in discriminant analysis 
using the penalized EM algorithm of Zhou et al (2009) (adapted in Sedki et al (2014) 
for the discriminant analysis settings).
}
\usage{
SortvarLearn(x, z, type, lambda, rho, nbcores)
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

\item{type}{character defining the type of ranking procedure, must be "lasso" or "likelihood". Default is "lasso"}


  \item{lambda}{
  numeric listing of tuning parameters for \eqn{\ell_1} mean penalty
}
  \item{rho}{
  numeric listing of tuning parameters for \eqn{\ell_1} precision matrix penalty
}

\item{nbcores}{
  number of CPUs to be used when parallel computing is utilized (default is 2)
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 vector of integers corresponding to variable ranking.
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
\author{
Mohammed Sedki \email{mohammed.sedki@u-psud.fr}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{SortvarClust}
}
\examples{
\dontrun{
## wine data set 
## n = 178 observations, p = 27 variables 
require(Rmixmod)
require(glasso)
data(wine)
set.seed(123)
obj <- SortvarLearn(x=wine[,1:27], z=wine[,28], nbcores=4)
}
}
\keyword{Variable ranking}
\keyword{Penalized discriminant analysis}
