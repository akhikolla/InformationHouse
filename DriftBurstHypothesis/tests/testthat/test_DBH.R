library(testthat)
library(DriftBurstHypothesis)

context("DBH C++ error tests")
test_that("DBH C++ error tests",{
  #these functions has broken before, hopefully they won't break again.
  expect_equal(DriftBurstHypothesis:::AsymptoticVarianceC(c(1:3), 3), NaN) 
  expect_equal(DriftBurstHypothesis:::AsymptoticVarianceC(c(1:3), 4), NaN)
  expect_equal(DriftBurstHypothesis:::AutomaticLagSelectionC(1:10, 30) , 7)
  }
)

context("DBH test")
test_that("DBH sim test", {
  set.seed(123)
  iT = 23400; meanBandwidth = 300L
  timestamps = seq(0, 23400, length.out = iT+1)
  testtimes  = seq(60, 23400, 60L)
  
  r = rnorm(iT, mean = 0.02, sd = 1)/sqrt(iT)
  p = c(0,cumsum(r))
  
  
  DBH = driftBursts(timestamps, p, testtimes, preAverage = 1, ACLag = -1, 
                     meanBandwidth = meanBandwidth, varianceBandwidth = 5*meanBandwidth, parallelize = FALSE)
  #The new mu conversion is 1/2*OldMu
  expect_equal(mean(getMu(DBH)),0.01253808)
  expect_equal(mean(getSigma(DBH)), 0.9748568)
})


context("Examples check")
test_that("DBH Examples check",{
  #Both a flash crash and flash rally are coded into the function.
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
  vY = 500+cumsum(flashCrashSim(iT, dSigma, dPhi, dMu))
  
  #insert an outlier to illustrate robustness.
  vY[50000] = 500
  
  #Here, the observations are equidistant, but the code can handle unevenly spaced observations.
  timestamps = seq(34200 , 57600 , length.out = iT)
  testtimes = seq(34260, 57600, 60)
  logprices = log(vY)
  
  library("DriftBurstHypothesis")
  
  #calculating drift burst hypothesis
  
  DBH = driftBursts(timestamps,  logprices,
                     testtimes, preAverage = 5, ACLag = -1L,
                     meanBandwidth = 300L, varianceBandwidth = 900L,
                     parallelize = FALSE)
  
  
  #plot test statistic
  plot = plot(DBH)
  #plot both test statistic and price
  plot2 = plot(DBH, price = vY, timestamps = timestamps)
  #Plot the mu series
  plot3 = plot(DBH, which = "Mu")
  #plot the sigma series
  plot4 = plot(DBH, which = "Sigma")
  expect_equal(plot3$usr[1:2], plot4$usr[1:2])
  #plot both
  plot5 = plot(DBH, which = c("Mu", "Sigma"))
  
  #Means of the tstat, sigma, and mu series.
  expect_equal(mean(getDB(DBH)), 0.01863012)
  expect_equal(mean(getSigma(DBH)), 0.007197401)
  expect_equal(mean(getMu(DBH)), -0.000008112) 
  
  
  
  
  ################## same example with xts object:
  suppressPackageStartupMessages(library("xts"))
  #Set parameter values of the simulation
  iT = 66500; dSigma = 0.3; dPhi = 0.98; dMu = -10;
  #set seed for reproducibility
  set.seed(123)
  #Simulate the series
  vY = 500+cumsum(flashCrashSim(iT, dSigma, dPhi, dMu))
  
  #insert an outlier to illustrate robustness.
  vY[50000] = 500
  
  #Here, the observations are equidistant, but the code can handle unevenly spaced observations.
  timestamps = seq(34200 , 57600 , length.out = iT)
  StartTime = strptime("1970-01-01 00:00:00.0000", "%Y-%m-%d %H:%M:%OS", tz = "GMT")
  Tradetime = StartTime + timestamps
  testTimes = seq(34260, 57600, 60)
  
  
  price = xts(vY, Tradetime)
  
  
  DBHxts = driftBursts(timestamps = NULL,  log(price), 
                       testTimes, preAverage = 5, ACLag = -1L,
                       meanBandwidth = 300L, varianceBandwidth = 900L, 
                       parallelize = FALSE)
  plot6 = plot(DBHxts)
  expect_true(all.equal.list(plot, plot6))
  plot7 = plot(DBHxts, price = price)
  expect_true(all.equal.list(plot2, plot7))
  #check for equality
  expect_true(all.equal(as.numeric(getDB(DBH)), as.numeric(getDB(DBHxts))))
})