\name{driftBursts}
\alias{driftBursts}

\title{
Drift Bursts
}
\description{
Calculates the Test-Statistic for the Drift Burst Hypothesis.
}
\usage{
driftBursts(timestamps = NULL, logPrices, testTimes = seq(34260, 57600, 60),
                        preAverage = 5, ACLag = -1L, 
                        meanBandwidth = 300L, varianceBandwidth = 900L, 
                        sessionStart = 34200, sessionEnd = 57600,
                        parallelize = FALSE, nCores = NA, warnings = TRUE)
}

\arguments{
  \item{timestamps}{
Either:
A \code{numeric} of \code{timestamps} for the trades in seconds after midnight.
Or: \code{NULL}, when the time argument is \code{NULL}, the \code{logprices} argument must be an \code{xts} object.
}
\item{logPrices}{
  A \code{numeric} or \code{xts} object containing the log-prices.
}
\item{testTimes}{
  A \code{numeric} containing the times at which to calculate the tests. The standard of \code{seq(34260, 57600, 60)} denotes calculating the test      -statistic once per minute, i.e. 390 times for a typical 6.5 hour trading day from 9:31:00 to 16:00:00.
  See details.
}
\item{preAverage}{
  An \code{integer} denoting the length of pre-averaging estimates of the log-prices.

}
\item{ACLag}{
  A positive \code{integer} greater than 1 denoting how many lags are to be used for the HAC estimator of the variance - the standard of \code{-1} denotes using an automatic lag selection algorithm for each iteration.
}
\item{meanBandwidth}{
  An \code{integer} denoting the bandwidth for the left-sided exponential kernel for the mean.
}
\item{varianceBandwidth}{
  An \code{integer} denoting the bandwidth for the left-sided exponential kernel for the variance.
}
\item{sessionStart}{
  A \code{double} to denote the start of the trading session in seconds after midnight.
}
\item{sessionEnd}{
  A \code{double} to denote the end of the trading session in seconds after midnight.
}
\item{parallelize}{
A \code{Boolean} to determine whether to parallelize the underlying C++ code. (Using OpenMP)
}
\item{nCores}{
An \code{integer} denoting the number of cores to use for calculating the code when paralellized.
If this argument is not provided, sequential evaluation will be used even though \code{bParallelize} is TRUE
}
\item{warnings}{
A \code{logical} denoting whether warnings should be shown.
}
}
\details{
If the \code{testTimes} vector contains instructions to test before the first trade, or more than 15 minutes after the last trade, these entries will be deleted, as not doing so may cause crashes.

The test statistic is unstable before \code{max(meanBandwidth , varianceBandwidth)} seconds has passed.

If \code{timestamps} is provided and \code{logPrices} is an \code{xts} object, the indices of logPrices will be used regardless.

Note that using an \code{xts} logPrices argument is slower than using a \code{numeric} due to the creation of the timestamps from the index of the input. When using \code{xts} objects, be careful to use the correct time zones. For example, if I as a dane use the \code{"America/New_York"} time zone for my \code{xts} objects, I have to add 14400 to my testing times. Same correction will have to be made to the \code{startTime} and \code{endTime} arguments in the plotting methods.


The lags from the Newey-West algorithm is increased by \code{2 * (preAveage-1)} due to the pre-averaging we know atleast this many lags should be corrected for. The maximum of 20 lags is also increased by this same factor for the same reason.


The methods \code{getMu(DBH)} and \code{getSigma(DBH)} have ad-hoc adjustments made to the output that are based on simulations. These adjustments give a slight down-ward bias in the value of the drift and sigma if the \code{preAverage} argument is VERY high, i.e. 10 or 20 for a noise-less price series of e.g. 23400 trades. If you find any better adjustments, please contact me.

To get the unadjusted series, use \code{DBH$mu} or \code{DBH$sigma} for the drift, and volatility respectively.
}

\value{
An object of class \code{DBH} and \code{list} containing the series of the drift burst hypothesis test-statistic as well as the estimated spot drift and variance series. The list also contains some information such as the variance and mean bandwidths along with the pre-averaging setting and the amount of observations. Additionally, the list will contain information on whether testing happened for all \code{testTimes} entries.
}

\section{warning}{
When using \code{getMu(DBH)}, beware that the series may be padded with zeros (a warning will be shown when running the \code{driftBursts} in this case). This will for example happen when using \code{testTimes} for an entire trading day but the data is from a half trading day. Therefore, methods for calculating the mean and variance of the series that take this into account are provided and documented \link[=getMu]{here}.
}
\references{
Christensen, Oomen and Reno (2018) <DOI:10.2139/ssrn.2842535>.
}
\author{
Emil Sjoerup
}

\seealso{
The DBH methods for the DBH class are documented \link[=getMu]{here}.
}

\examples{
#Simulate a price series containing both a flash crash and a flash-rally
flashCrashSim = function(iT, dSigma, dPhi, dMu){
  vSim = numeric(iT)
  vEps = rnorm(iT , sd =dSigma)
  vEpsy = rnorm(iT)
  vEps[30001:32000] = rnorm(2000 ,sd =seq(from = dSigma , 
                                          to = 2*dSigma , length.out = 2000)) 
  vEps[32001:34000] = rnorm(2000 ,sd =seq(from = 2*dSigma , 
                                          to = dSigma , length.out = 2000))
  vEpsy[30001:32000] = -rnorm(2000 ,mean =seq(from = 0,
                                              to = 0.3 , length.out = 2000)) 
  vEpsy[32001:34000] = -rnorm(2000 ,mean =seq(from = 0.3,
                                              to = 0 , length.out = 2000))
  
  
  vEps[60001:63000] = rnorm(3000,sd = seq(from = dSigma , 
                                          to = 2*dSigma , length.out = 3000))
  vEps[63001:66000] = rnorm(3000,  sd = seq(from = 2*dSigma , 
                                            to =  dSigma, length.out = 3000))
  
  vEpsy[60001:63000] = rnorm(3000 ,mean =seq(from = 0,
                                             to = 0.2 , length.out = 3000))
  vEpsy[63001:66000] = rnorm(3000 ,mean =seq(from = 0.2,
                                             to = 0 , length.out = 3000))
  vSim[1] = dMu + dPhi *rnorm(1 , mean = dMu , sd = dSigma /sqrt(1-dPhi^2))
  for (i in 2:iT) {
    vSim[i] = dMu + dPhi * (vSim[(i-1)] - dMu) + vEps[i]
  }
  vY = exp(vSim/2) * vEpsy
  return(vY)
}
#Set parameter values of the simulation
iT = 66500; dSigma = 0.3; dPhi = 0.98; dMu = -10;
#set seed for reproducibility
set.seed(123)
#Simulate the series
y = 500+cumsum(flashCrashSim(iT, dSigma, dPhi, dMu))

#insert an outlier to illustrate robustness.
y[50000] = 500

#Here, the observations are equidistant, but the code can handle unevenly spaced observations.
timestamps = seq(34200 , 57600 , length.out = iT)
testTimes = seq(34260, 57600, 60)
logPrices = log(y)

library("DriftBurstHypothesis")

#calculating drift burst test statistic

DBH = driftBursts(timestamps,  logPrices,
                        testTimes, preAverage = 5, ACLag = -1L,
                        meanBandwidth = 300L, varianceBandwidth = 900L,
                        parallelize = FALSE)


#plot test statistic
plot(DBH)
#plot both test statistic and price
plot(DBH, price = y, timestamps = timestamps)
#Plot the mu series
plot(DBH, which = "Mu")
#plot the sigma series
plot(DBH, which = "Sigma")
#plot both
plot(DBH, which = c("Mu", "Sigma"))


#Means of the tstat, sigma, and mu series.
mean(getDB(DBH))
mean(getSigma(DBH))
mean(getMu(DBH))



\dontrun{
################## same example with xts object:
library("xts")
#Set parameter values of the simulation
iT = 66500; dSigma = 0.3; dPhi = 0.98; dMu = -10;
#set seed for reproducibility
set.seed(123)
#Simulate the series
y = 500+cumsum(flashCrashSim(iT, dSigma, dPhi, dMu))

#insert an outlier to illustrate robustness.
y[50000] = 500

#Here, the observations are equidistant, but the code can handle unevenly spaced observations.
timestamps = seq(34200 , 57600 , length.out = iT)
startTime = strptime("1970-01-01 00:00:00.0000", "\%Y-\%m-\%d \%H:\%M:\%OS", tz = "GMT")
testTimes = seq(34260, 57600, 60)


price = xts(vY, startTime + timestamps)


DBHxts = driftBursts(timestamps = NULL,  log(price), 
                   testTimes, preAverage = 5, ACLag = -1L,
                   meanBandwidth = 300L, varianceBandwidth = 900L, 
                   parallelize = FALSE)

#check for equality
all.equal(as.numeric(getDB(DBH)), as.numeric(getDB(DBHxts)))
}
}