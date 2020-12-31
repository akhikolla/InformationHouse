driftBursts = function(timestamps = NULL, logPrices, testTimes = seq(34260, 57600, 60),
                       preAverage = 5, ACLag = -1L, meanBandwidth = 300L, 
                       varianceBandwidth = 900L, sessionStart = 34200, sessionEnd = 57600,
                       parallelize = FALSE, nCores = NA, warnings = TRUE){

  #########  Initialization  ############
  k                    = preAverage 
  vDiff                = diff(logPrices)
  iT                   = length(logPrices)
  xts                  = FALSE
  pad = removedFromEnd = 0
  tt                   = testTimes #tt is returned in the Info list. 
  #########  init end  ############
  
  ###Checks###
  if (meanBandwidth<0 | meanBandwidth %% 1 != 0) {
    stop("meanBandwidth must be a positive integer")
  }
  if(varianceBandwidth<0 | varianceBandwidth %% 1 != 0){
    stop("varianceBandwidth must be a positive integer")
  }
  if(ACLag !=-1 && ACLag%%1!=0 | -1>ACLag | ACLag == 1){
    stop("ACLag must be a positive integer greater than 1, or -1. \n
         The standard of -1 designates usage of an automated lag selection algorithm.")
    #Specifically Newey-West 1994
  }
  if(preAverage <=0 | preAverage%%1!=0 ){
    stop("preAverage must be a positive integer. No preaveraging is done when preAverage = 1.")
  }
  if(inherits(logPrices, "xts")){
    tz              = tzone(logPrices)
    timestamps      = index(logPrices)
    timestamps      = as.numeric(timestamps) - (.indexDate(logPrices)[1] * 86400)
    vIndex          = as.POSIXct((.indexDate(logPrices)[1] * 86400) + testTimes, origin = "1970-01-01", tz = tz)
    logPrices       = as.numeric(logPrices)
    #logPrices       = as.numeric(t(logPrices)) ##need to transpose, otherwise the program will crash.
    vDiff              = as.numeric(vDiff)[-1] ### need to remove first entry because diff() on an xts object produces NA in first entry.
    xts             = TRUE
  }
  if((anyNA(timestamps) & !is.null(timestamps)) | anyNA(logPrices) | anyNA(testTimes)){
    stop("NA's in timestamps, logPrices or testTimes - might cause crashes and are thus disallowed")
  }
  if((length(timestamps) != length(logPrices) & !is.null(timestamps))){
    stop("timestamps and logPrices input not of same length, to prevent crashing this is not allowed.")
  }
  if((is.na(nCores) | nCores %% 1 != 0) & parallelize){
    if(warnings){
    warning("No iCores argument was provided, or the provided nCores argument is not an integer.\n
             Sequential evaluation is used.")
  }
    parallelize = FALSE
  }
  if(min(timestamps) + 5 >= min(testTimes)){ ## Safeguard to prevent faulty inputs that causes crashes!
    testTimes = testTimes[-1]
    pad = 1
    while(min(timestamps) + 5 >= min(testTimes)){
      testTimes = testTimes[-1] 
      pad = pad + 1
    }
    if(warnings){
      warning(paste("\nThe first testing timestamps is before any observations. This causes fatal errors (and crashes if paralellized).",
                    "\nItereatively removing first testing timestamps until this is no longer the case.",
                    "\nRemoved the first", pad, "entries from testTimes and replacing with a 0\n"))
    }
    
    if(length(testTimes) == 0){
      stop("No testing done, the mandated testing times were all too soon after market open. 
           A five second minimum wait is put in as a safe-guard to prevent inputs that may cause crashes.")
    }
    
    
  }
  if(max(testTimes)>max(timestamps) + 900){
    testTimes = testTimes[-length(testTimes)]
    removedFromEnd = 1
    while(max(testTimes) >= max(timestamps) + 900){
      testTimes = testTimes[-length(testTimes)]  
      removedFromEnd = removedFromEnd + 1
    }
    if(warnings){
      warning(paste("\nThe last testing timestamps is more than 15 minutes after the last trade, this may cause fatal errors or crashes.",
                    "\nIteratively removing the last testing timestamps until this is no longer the case.",
                    "\nRemoved the last", removedFromEnd, "entries from testTimes\n"))
      
    }
    
  }
  if(15 <= max(diff(timestamps))/60){
    stop("There is a period of greater than 15 minutes with no trades.\n
          This causes fatal errors (and crashes if paralellized) and is thus disallowed")
  }
  ###Checks end###
  vpreAveraged = rep(0, iT - 1)
  vpreAveraged[(k*2 - 1):(iT - 1)] = cfilter(x = logPrices, c(rep(1,k),rep( -1,k)))[k:(iT - k)]
  
  if(parallelize & !is.na(nCores)){ #Parallel evaluation or not?
    lDriftBursts = DriftBurstLoopCPAR(vpreAveraged, vDiff, timestamps, testTimes, meanBandwidth, 
                                     varianceBandwidth, preAverage, ACLag, nCores)
  }else{
    lDriftBursts = DriftBurstLoopC(vpreAveraged, vDiff, timestamps, testTimes, meanBandwidth, 
                                  varianceBandwidth, preAverage, ACLag)  
  }



  if(pad != 0 | removedFromEnd != 0){
    lDriftBursts[["driftBursts"]] = c(rep(0,pad), lDriftBursts[["driftBursts"]], rep(0,removedFromEnd))
    lDriftBursts[["sigma"]]       = c(rep(0,pad), lDriftBursts[["sigma"]], rep(0,removedFromEnd))
    lDriftBursts[["mu"]]          = c(rep(0,pad), lDriftBursts[["mu"]], rep(0,removedFromEnd))
  }
  
  if(xts){
    lDriftBursts[["driftBursts"]] = xts(lDriftBursts[["driftBursts"]], order.by = vIndex, tz = tz)
    lDriftBursts[["sigma"]]       = xts(lDriftBursts[["sigma"]], order.by = vIndex, tz = tz)
    lDriftBursts[["mu"]]          = xts(lDriftBursts[["mu"]], order.by = vIndex, tz = tz)
  }
  
  lInfo = list("varianceBandwidth" = varianceBandwidth, "meanBandwidth" = meanBandwidth,"preAverage" = preAverage,
               "nObs" = iT, "testTimes" = tt, "padding" = c(pad, removedFromEnd), "sessionStart" = sessionStart, "sessionEnd" = sessionEnd)
  lDriftBursts[["info"]] = lInfo
  #replace NANs with 0's
  NANS = is.nan(lDriftBursts[["sigma"]])
  lDriftBursts[["driftBursts"]][NANS] = 0
  lDriftBursts[["sigma"]][NANS]       = 0
  
  class(lDriftBursts) = c("DBH", "list")
  return(lDriftBursts)
}
