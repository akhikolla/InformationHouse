\name{sim.data}
\alias{sim.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Data from a Variety of Functional Scenarios
}
\description{
This function generates data according to the simulation scenarios considered in Section 5 and plotted in Figure 2 of Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391. Each scenario has four covariates that have some non-linear association with the outcome. There is the option to also generate a user-specified number of covariates that have no association with the outcome.}
\usage{
sim.data(n, scenario, zerof, noise = 1, family = "gaussian")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations.
}
  \item{scenario}{
simulation scenario to use. Options are 1, 2, 3, or 4, which correspond to the simulation scenarios of Section 5 in Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391. Each scenario has four covariates.
}
  \item{zerof}{
number of noise covariates (those that have no relationship to the outcome) to include. This can be used to replicate the high-dimensional scenarios of Section 5 in Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391. The total number of covariates will be \code{4 + zerof}.
}
  \item{noise}{
the variance of the errors. If \code{family = "gaussian"}, the errors of observations are generated using MVN(0, \code{noise}I) with default of 1, otherwise \code{noise} is not used if \code{family="binomial"}. 
}
  \item{family}{
the error distribution of observations (must be \code{family="gaussian"} or \code{family="binomial"}). If \code{family = "gaussian"}, the errors of observations are generated using MVN(0, \code{noise}I), otherwise observations are Bernoulli if \code{family="binomial"}. 
}
}
\value{
\item{x}{n x p covariate matrix.}
\item{y}{n-vector containing the outcomes for the n observations in \code{x}.}
\item{theta}{n x p mean matrix used to generate \code{y}.}
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data to fit FLAM model with squared-error loss
set.seed(1)
data <- sim.data(n = 50, scenario = 1, zerof = 10, noise = 1)
flam.out <- flam(x = data$x, y = data$y, family = "gaussian")

#alternatively, generate data for logistic FLAM model 
#note: 'noise' argument no longer needed
data2 <- sim.data(n = 50, scenario = 1, zerof = 0, family = "binomial")
flam.logistic.out <- flam(x = data2$x, y = data2$y, family = "binomial")


#vary generating functions
#choose large n because we want to plot generating functions
data1 <- sim.data(n = 500, scenario = 1, zerof = 0)
data2 <- sim.data(n = 500, scenario = 2, zerof = 0)
data3 <- sim.data(n = 500, scenario = 3, zerof = 0)
data4 <- sim.data(n = 500, scenario = 4, zerof = 0)
#and plot to see functional forms
par(mfrow=c(2,2))
col.vec = c("dodgerblue1","orange","seagreen1","hotpink")
for (i in 1:4) {
	if (i==1) data = data1 else if (i==2) data = data2 
		else if (i==3) data = data3 else data = data4
	plot(1,type="n",xlim=c(-2.5,2.5),ylim=c(-3,3),xlab=expression(x[j]),
		ylab=expression(f[j](x[j])),main=paste("Scenario ",i,sep=""))
	sapply(1:4, function(j) points(sort(data$x[,j]), 
		data$theta[order(data$x[,j]),j],col=col.vec[j],type="l",lwd=3))
}

#include large number of predictors that have no relationship to outcome
data <- sim.data(n = 50, scenario = 1, zerof = 100, noise = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
