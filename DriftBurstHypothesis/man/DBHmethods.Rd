\name{methods for DBH objects}
\alias{print.DBH}
\alias{plot.DBH}
\alias{getMu.DBH}
\alias{getSigma.DBH}
\alias{getDB.DBH}
\alias{getMean.DBH}
\alias{getVar.DBH}
\alias{getMu}
\alias{getSigma}
\alias{getDB}
\alias{getMean}
\alias{getVar}
\alias{getCriticalValues}
\alias{getCriticalValues.DBH}
\title{
  Methods for the DBH class
}
\description{
These are the methods for the class DBH which are currently implemented.
}

\usage{
\method{print}{DBH}(x, ...)
\method{plot}{DBH}(x, ...)
\method{getDB}{DBH}(x)
\method{getCriticalValues}{DBH}(x, alpha)
\method{getMu}{DBH}(x, annualize = FALSE, nDays = 252)
\method{getSigma}{DBH}(x, annualize = FALSE, nDays = 252)
\method{getMean}{DBH}(x, which = 'all')
\method{getVar}{DBH}(x, which = 'all', annualize = FALSE, nDays = 252)
}



\arguments{
\item{x}{DBH object}
\item{...}{Additional arguments for the plotting and printing routines. See details}
\item{alpha}{\code{double} confidence level for the critical values to be extracted. Default = 0.95 designating 95\%}
\item{annualize}{\code{logical} determining whether or not to annualize the series.}
\item{nDays}{\code{numeric} determining how many days to use for annualization. The standard of 252 reflects the average 252 trading days in a year.}
\item{which}{\code{character} determining which series to retrieve the mean and variance of in \code{getMean} and \code{getVar}. This argument can also be used in the \code{print} and \code{plot} methods, see details for use in these methods.}
}


\details{
For the \code{print} method, the \code{...} argument can be used to pass the following arguments:
\itemize{
\item{\code{criticalValue}: \code{numeric} determining the critical value to use when determining whether a drift burst is present in the data. If this critical value is omitted, a method to determine the critical value based on a interpolations of simulated data-set of critical values for different confidence levels and auto-correlations of the test statistic. See appendix B in the article.}
\item{\code{alpha}: \code{numeric} of length 1 which determines the confidence level of critical values extracted using the method above.}
}
For the \code{plot} method, the \code{...} argument can be used to pass the following arguments:
\itemize{
  \item{\code{which}: \code{character} Used for choosing which series to plot. Valid choices are: \code{"DriftBursts"}, \code{"DB"}, \code{"Sigma"}, \code{"Mu"}, and     \code{c("Sigma","Mu")}, the order of the latter is irrelevant. The case of the input does not matter. Default = \code{"driftbursts"}}
  \item{\code{price}:The price series which, if provided, will be overlayed in a red dotted line and the level will be shown at the right y-axis.       (Only   used if \code{which} is \code{"DriftBursts"}). Default = \code{NULL}}
  \item{\code{time}:Timestamps for the trades in seconds after midnight, which will be used for the x-axis of the plot if the price is overlayed.       Default = \code{NULL}}
  \item{\code{startTime}:Start of the trading day in seconds after midnight. Default = \code{34200}}
  \item{\code{endTime}:End of the trading day in seconds after midnight. Default = \code{57600}}
  \item{\code{leg.x}:X-position of the legend in the case \code{which} is \code{"DriftBursts"} AND the price is overlayed. Default = \code{"topleft"}.   Usage is as in the base \code{R} engine.}
  \item{\code{leg.y}:Y-position of the legend in the case \code{which} is \code{"DriftBursts"} AND the price is overlayed. Default = \code{NULL}.     Usage is as in the base \code{R} engine.}
  \item{\code{tz}: \code{character} denoting the time-zone. Default = \code{"GMT"}}
  \item{\code{annualize}: \code{logical} denoting whether to annualize in case which contains \code{"mu"}, \code{"sigma"}, or both. Default = \code{FALSE}}
  \item{\code{nDays}: \code{numeric} denoting how many to use for annualization if \code{annualize} is \code{TRUE}. Default = \code{252}}
  }
}
\value{
For \code{print} and \code{plot}, nothing is returned.

For \code{getCriticalValues}, the critical value, and the normalized critical value is returned in a \code{list}.

For  \code{getDB}, \code{getMu}, and \code{getSigma}, vectors with the same length as testTimes, containing the test statistic, drift, or volatility respecitvely, is returned

For the methods \code{getMean} and \code{getVar}:

When the \code{which} argument is \code{'all'}, getMean and getVar returns lists containing the mean or variance of the test statistic, drift, and volatility respectively.

When the \code{which} argument is \code{'db'} (or \code{'driftbursts'}), \code{'mu'}, or \code{'sigma'}, doubles containing the mean or variance is returned.

Note that the reason for having the getMean and getVar methods is the possibility of testing being mandated while data is not present. This primarily happens when loading in tick data for an entire year in a loop, including for example Christmas Eve in the U.S. where the markets close early. In this case if the user reuses the same \code{testTimes}, a warning will be thrown in \code{driftBursts}, and no testing will be done, but the output series will be padded with zeros to keep the same size as \code{testTimes}. 
}
\author{
Emil Sjoerup
}

\examples{
library(DriftBurstHypothesis)
set.seed(1)
# Set mean and variance bandwidth parameters
meanBandwidth = 300L
varianceBandwidth = 900L

# Simulate noise-less price series with 23400 observations, denoting 1 trader per second
# and generate corresponding timestamps.
iT = 23399
r = rnorm(iT, mean = 0, sd = 1)/sqrt(iT)
p = c(0,cumsum(r))
timestamps = seq(34200, 57600, length.out = iT+1)

# Test every minute after the instability period is over.
testTimes  = seq(34260 + varianceBandwidth, 57600, 60L)

# Calculate drift burst test statistic
DBH = driftBursts(timestamps, p, testTimes, preAverage = 1, ACLag = -1,
                  meanBandwidth = meanBandwidth, varianceBandwidth = varianceBandwidth)




print(DBH)
# Plots the test statistic with prices overlaid.
plot(DBH, timestamps = timestamps, price = p)

# Plots the annualized volatility
plot(DBH, which = 'sigma', annualize = TRUE)

# Plots the annualized drift and volatility
plot(DBH, which = c('sigma', 'mu'), annualize = TRUE)

# Retrieve the critical values of the drift burst test statistic
getCriticalValues(DBH)

# Calculate the mean of the test statistic, drift, and volatility
getMean(DBH, which = 'all')

# Calculate the variance of the test statistic
getVar(DBH, which = 'db')

# Extracts the annualized drift
annualizedDrift = getMu(DBH, annualize = TRUE)

# Extracts the annualized volatility
annualizedVolatility = getSigma(DBH, annualize = TRUE)
}


