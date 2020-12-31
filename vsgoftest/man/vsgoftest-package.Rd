\name{vsgoftest-package}
\alias{vsgoftest-package}
\alias{vsgoftest}
\docType{package}
\title{
\packageTitle{vsgoftest}
}
\description{
\packageDescription{vsgoftest}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{vsgoftest}
\packageIndices{vsgoftest}
}
\author{
\packageAuthor{vsgoftest}

Maintainer: \packageMaintainer{vsgoftest}
}
\references{

Vasicek, O., A test for normality based on sample entropy, \emph{Journal of the Royal Statistical Society,} \bold{38(1)}, 54-59 (1976).

Song, K. S., Goodness-of-fit tests based on Kullback-Leibler discrimination information, \emph{Information Theory, IEEE Transactions on,} \bold{48(5)}, 1103-1117 (2002).

Girardin, V., Lequesne, J. Entropy-based goodness-of-fit tests - a unifying framework. Application to DNA replication. \emph{Communications in Statistics: Theory and Methods} (2017). https://doi.org/10.1080/03610926.2017.1401084

Lequesne, J., Regnault, P. vsgoftest: An R Package for Goodness-of-Fit Testing Based on Kullback-Leibler Divergence. \emph{Journal of Statistical Software}, \bold{96} (2020).  doi:10.18637/jss.v096.c01 

}
\keyword{ package }

\examples{
set.seed(1)
samp <- rnorm(50, mean = 2, s = 3)

##Estimating entropy
entropy.estimate(x = samp, window = 8)
log(2*pi*exp(1))/2 #true value of entropy of normal distribution

##Testing normality
vs.test(x = samp, densfun = 'dnorm', param = c(2,3), B = 500) #Simple null hypothesis
vs.test(x = samp, densfun='dnorm', B = 500) #Composite null hypothesis

}
