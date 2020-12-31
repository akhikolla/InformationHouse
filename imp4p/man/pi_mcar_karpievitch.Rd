\name{pi.mcar.karpievitch}
\alias{pi.mcar.karpievitch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the proportion of MCAR values in biological conditions using the method of Karpievitch (2009).
}
\description{
This function allows estimating the proportion of MCAR values in biological conditions using the method of Karpievitch (2009).
}
\usage{
pi.mcar.karpievitch(tab,conditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list composed of:
\item{pi.mcar}{The proportion of MCAR values in each biological condition.}
\item{prop.na}{The proportion of missing values for each peptide in each condition.}
\item{moy}{The average of observed values for each peptide in each condition.}
}
\references{
Karpievitch, Y., Stanley, J., Taverner, T., Huang, J., Adkins, J. N., Ansong, C., ... & Smith, R. D. (2009). A statistical framework for protein quantitation in bottom-up MS-based proteomics. Bioinformatics, 25(16), 2028-2034.}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estim.mix}}
}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600);

#Proportion of MCAR values in each condition
pi.mcar.karpievitch(tab=res.sim$dat.obs,conditions=res.sim$conditions)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
