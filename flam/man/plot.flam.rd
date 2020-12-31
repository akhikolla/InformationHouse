\name{plot.flam}
\alias{plot.flam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots Function Estimates for Fit of Class "flam"
}
\description{
This function plots the estimated functions from a model estimated using \code{\link{flam}}. The user specifies the model of interest (i.e., the tuning parameters) and a plot is made for the estimated association between all non-sparse (or a chosen subset of) features and the outcome.}
\usage{
\method{plot}{flam}(x, index, n.plot = 10, predictor.indicators = NULL,
predictor.labels = NULL, outcome.label = "outcome", ticks = F, 
col = "dodgerblue", n.panel.width = NULL, n.panel.height = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "flam".
}
  \item{index}{
the index for the model of interest to be plotted. Note that \code{index} of i corresponds to the model with tuning parameters \code{x$all.alpha[i]} and \code{x$all.lambda[i]}.
}
  \item{n.plot}{
the number of predictors to be plotted (default of 10). Note that only non-sparse predictors are plotted, however, this argument is ignored if \code{predictor.indicators} is specified.
}
  \item{predictor.indicators}{
a vector indicating which predictor function estimates to plot. The vector should contain the column numbers of \code{x$x} whose function estimates are to be plotted. By default, the \code{n.plot} predictors with the largest L2 norm are plotted.
}
  \item{predictor.labels}{
a vector containing the predictor labels of the same length as \code{predictor.indicators} if specified, otherwise of length \code{ncol(x$x)}. By default, "Predictor 1", "Predictor 2", ... are used.
}
  \item{outcome.label}{
the name of the outcome used in the y-axis label. By default, the label is "outcome".
}
  \item{ticks}{
a logical (\code{TRUE} or \code{FALSE}) for whether tick marks indicating the distribution of the predictor should be included at the bottom of each plot.
}
  \item{col}{
a vector of colors used to plot the function estimates. If \code{col} has length 1, then the same color is used for all predictors. Otherwise, \code{col} should indicate the color for each predictor and have the same length as \code{predictor.indicators} if specified, otherwise the number of columns of \code{x$x}.
}
  \item{n.panel.width}{
the plots will be plotted with \code{par(mfrow=c(n.panel.height,n.panel.width))}. If specified, \code{n.panel.height} must also be specified. By default, an appropriate grid of plots is used.
}
  \item{n.panel.height}{
the plots will be plotted with \code{par(mfrow=c(n.panel.height,n.panel.width))}. If specified, \code{n.panel.width} must also be specified. By default, an appropriate grid of plots is used.
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}
\details{
The estimated function fits are drawn by connecting the predictions for all of the observations in \code{x$x}. This may result in fits that appear not to be piecewise consant if only a single observation is observed for a range of x-values.
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\seealso{
\code{\link{flam}}
}
\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data
set.seed(1)
data <- sim.data(n = 50, scenario = 1, zerof = 10, noise = 1)
#fit model for a range of tuning parameters
flam.out <- flam(x = data$x, y = data$y, alpha.seq = c(0.8, 0.9, 1))

#we plot the predictor fits for a specific index, e.g. 25
#that is, lambda and alpha of
flam.out$all.lambda[25]; flam.out$all.alpha[25]
plot(flam.out, index = 25)
#the fit only has 5 non-sparse features

#by default, up to 10 non-sparse features with the largest L2 norms are 
#plotted, but we can plot a different number of features if desired
plot(flam.out, index = 40, n.plot = 12)
#or we can plot specific predictors of interest
plot(flam.out, index = 40, predictor.indicators = c(1:4, 6, 8, 11, 12))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
