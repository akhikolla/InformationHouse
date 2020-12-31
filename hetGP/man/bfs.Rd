\name{bfs}
\alias{bfs}
\alias{bfs.exp}
\alias{bfs.gamma}

\docType{data}
\title{ Bayes Factor Data }
\description{
  Data from a Bayes factor MCMC-based simulation experiment comparing Student-t to Gaussian errors in an RJ-based Laplace prior Bayesian linear regession setting
}
\usage{data(ato)}
\format{
  Calling \code{data(bfs)} causes the following objects to be loaded into the namespace.

  \describe{
    \item{\code{bfs.exp}}{ 20x11 \code{data.frame} whose first column is \eqn{\theta}{theta}, indicating the mean parameter of an exponential distribution encoding the prior of the Student-t degrees of freedom parameter \eqn{\nu}{nu}.  The remaining ten
    columns comprise of Bayes factor evaluations under that setting  }
    \item{\code{bfs.gamma}}{ 80x7 \code{data.frame} whose first two columns are \eqn{\beta}{beta} and \eqn{\alpha}{alpha}, indicating the second and first parameters to a 
    Gamma distribution encoding the prior of the Student-t degrees of freedom parameters \eqn{\nu}{nu}.  The remaining five columns comprise of Bayes factor evaluations under those settings }
  }
}
\details{
  Gramacy & Pantaleo (2010), Sections 3-3-3.4, describe an experiment
  involving Bayes factor (BF)  calculations to determine if data are
  leptokurtic (Student-t errors) or not (simply Gaussian) as a function of the
  prior parameterization on the Student-t degrees of freedom parameter
  \eqn{\nu}{nu}. Franck & Gramacy (2018) created a grid of hyperparameter
  values in \eqn{\theta}{\theta} describing the mean of an Exponential
  distribution, evenly spaced in \eqn{\log_{10}}{log10} space from
  \code{10^(-3)} to \code{10^6} spanning \dQuote{solidly Student-t} (even
  Cauchy) to \dQuote{essentially Gaussian} in terms of the mean of the prior
  over \eqn{\nu}{nu}.  For each \eqn{\theta}{theta} setting on the grid they
  ran the Reversible Jump (RJ) MCMC to approximate the BF of Student-t over Gaussian 
  by feeding in sample likelihood evaluations provided by \pkg{monomvn}'s
  \code{\link[monomvn]{blasso}} to compute the BF. In order to understand the
  Monte Carlo variability in those calculations, ten replicates of the BFs
  under each hyperparameter setting were collected.  These data are provided
  in \code{bfs.exp}.

  A similar, larger experiment was provided with \eqn{\nu}{nu} under a Gamma
  prior with parameters \eqn{\alpha}{alpha} and \eqn{\beta \equiv \theta}{beta
  = theta}.  In this higher dimensional space, a Latin hypercube sample 
  of size eighty was created, and five replicates of BFs were recorded.  
  These data are provided in \code{bfs.gamma}.

  The examples below involve \code{\link{mleHetTP}} fits (Chung, et al., 2018)
  to these data and a visualization of the predictive surfaces thus obtained.
  The code here follows an example provided, with more detail, in
  \code{vignette("hetGP")} 
}

\note{ For code showing how these BFs were calculated, see supplementary material
from Franck & Gramacy (2018) }

\author{ 
Mickael Binois, \email{mbinois@mcs.anl.gov}, and
Robert B. Gramacy, \email{rbg@vt.edu}
}

\references{ 

  Franck CT, Gramacy RB (2018). Assessing Bayes factor surfaces using
  interactive visualization and computer surrogate modeling. Preprint
  available on arXiv:1809.05580.

  Gramacy RB (2017). \pkg{monomvn}: Estimation for Multivariate Normal
  and Student-t Data with Monotone Missingness. R package version 1.9-7,
  \url{https://CRAN.R-project.org/package=monomvn}. 

  R.B. Gramacy and E. Pantaleo (2010). Shrinkage regression for multivariate
  inference with missing data, and an application to portfolio balancing.
  Bayesian Analysis 5(2), 237-262. Preprint available on arXiv:0907.2135

  Chung M, Binois M, Gramacy RB, Moquin DJ, Smith AP, Smith AM (2018).
  Parameter and Uncertainty Estimation for Dynamical Systems Using Surrogate
  Stochastic Processes. Preprint available on arXiv:1802.00852. 

}

\examples{
data(bfs)

##
## Exponential version first
##

thetas <- matrix(bfs.exp$theta, ncol=1)
bfs <- as.matrix(t(bfs.exp[,-1]))

## the data are heavy tailed, so t-errors help
bfs1 <- mleHetTP(X=list(X0=log10(thetas), Z0=colMeans(log(bfs)),
  mult=rep(nrow(bfs), ncol(bfs))), Z=log(as.numeric(bfs)), lower=10^(-4), 
  upper=5, covtype="Matern5_2")

## predictions on a grid in 1d
dx <- seq(0,1,length=100)
dx <- 10^(dx*4 - 3)
p <- predict(bfs1, matrix(log10(dx), ncol=1))

## visualization
matplot(log10(thetas), t(log(bfs)), col=1, pch=21, ylab="log(bf)", 
  main="Bayes factor surface")
lines(log10(dx), p$mean, lwd=2, col=2)
lines(log10(dx), p$mean + 2*sqrt(p$sd2 + p$nugs), col=2, lty=2, lwd=2)
lines(log10(dx), p$mean - 2*sqrt(p$sd2 + p$nugs), col=2, lty=2, lwd=2)
legend("topleft", c("hetTP mean", "hetTP interval"), lwd=2, lty=1:2, col=2)

##
## Now Gamma version
##

D <- as.matrix(bfs.gamma[,1:2])
bfs <- as.matrix(t(bfs.gamma[,-(1:2)]))

## fitting in 2fd
bfs2 <- mleHetTP(X=list(X0=log10(D), Z0=colMeans(log(bfs)), 
  mult=rep(nrow(bfs), ncol(bfs))), Z = log(as.numeric(bfs)), 
  lower = rep(10^(-4), 2), upper = rep(5, 2), covtype = "Matern5_2", 
  maxit=100000)

## predictions on a grid in 2d
dx <- seq(0,1,length=100)
dx <- 10^(dx*4 - 3)
DD <- as.matrix(expand.grid(dx, dx))
p <- predict(bfs2, log10(DD))

## visualization via image-contour plots
par(mfrow=c(1,2))
mbfs <- colMeans(bfs)
image(log10(dx), log10(dx), t(matrix(p$mean, ncol=length(dx))),  
  col=heat.colors(128), xlab="log10 alpha", ylab="log10 beta", 
  main="mean log BF")
text(log10(D[,2]), log10(D[,1]), signif(log(mbfs), 2), cex=0.5)
contour(log10(dx), log10(dx),t(matrix(p$mean, ncol=length(dx))),
  levels=c(-5,-3,-1,0,1,3,5), add=TRUE, col=4)
image(log10(dx), log10(dx), t(matrix(sqrt(p$sd2 + p$nugs), 
  ncol=length(dx))),  col=heat.colors(128), xlab="log10 alpha", 
  ylab="log10 beta", main="sd log BF")
text(log10(D[,2]), log10(D[,1]), signif(apply(log(bfs), 2, sd), 2), 
  cex=0.5)
}

\seealso{ 
  \code{\link{ato}}, \code{\link{sirEval}}, \code{\link{mleHetTP}}, 
  \code{vignette("hetGP")}, \code{\link[monomvn]{blasso}}
}

\keyword{datasets}
\keyword{datagen}