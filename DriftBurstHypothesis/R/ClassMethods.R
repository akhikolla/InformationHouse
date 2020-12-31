plot.DBH = function(x, ...){
  
  #### Get extra passed options and data
  options = list(...)
  #### List of standard options
  opt = list(which = "driftbursts", price = NULL, timestamps = NULL, 
             startTime = ifelse(is.null(x$info[['sessionStart']]), min(x$info[['testTimes']]), x$info[['sessionStart']]), 
             endTime = ifelse(is.null(x$info[['sessionEnd']]), max(x$info[['testTimes']]), x$info[['sessionEnd']]),
             leg.x = "topleft", leg.y = NULL,  tz = "GMT", annualize = FALSE, nDays = 252, legend.txt = "")
  #### Override standard options where user passed new options
  opt[names(options)] = options

  #### Extract options (better way to do it?) 
  which      = tolower(opt$which)
  startTime  = opt$startTime
  endTime    = opt$endTime
  main       = opt$main
  tz         = opt$tz
  leg.x      = opt$leg.x
  leg.y      = opt$leg.y
  annualize  = opt$annualize
  nDays      = opt$nDays
  price      = opt$price
  timestamps = opt$timestamps
  tstat      = getDB(x)
  sigma      = getSigma(x, annualize, nDays) 
  mu         = getMu(x, annualize, nDays)
  startpar   = par(no.readonly = TRUE)
  testTimes  = x$info$testTimes
  horizLines = seq(round(min(tstat)), round(max(tstat)), 1)
  ###Setup done
  if(!all(which %in% c("driftbursts", "mu", "sigma", "db"))){
    stop("The which argument must be a character vector containing either:\n
         Sigma, Mu, both of these or DriftBursts. 
         CasE doesn't matter.")
  }
  if(inherits(tstat, "xts")){
    tstat     = as.numeric(tstat)
    sigma     = as.numeric(sigma)
    mu        = as.numeric(mu)
    if(!is.null(price)){
      tz          = tzone(price)
      timestamps  = index(price, tz = tz)
      timestamps  = as.numeric(timestamps) - (.indexDate(price)[1] * 86400)
      price       = as.numeric(price)
    }
  }
  if(testTimes[1] == startTime){
    testTimes = testTimes[-1]
    sigma = sigma[-1]
    mu = mu[-1]
    tstat=tstat[-1]
  }
  if(min(testTimes) < startTime | max(testTimes) > endTime){
    cat('\nTesting was tried before sessionStart or after sessionEnd, thus some of the tests may be cut off from the plot.
        \nIf the plot looks weird, consider changing sessionStart and sessionEnd. 
        \nThese should reflect the start of trading and the end of trading respectively')
  }
  xtext = as.POSIXct(testTimes, origin = "1970-01-01", tz = tz)
  xlim = c(startTime, endTime)
  xlab = "Time"
  
  if(all(which %in% c("driftbursts", "db"))){ #use all() because this function should accept which arguments with length longer than 1
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    if(!is.null(price)) par(mar = c(4,3.5,4,4), mgp = c(2,1,0)) #makes room for values on the right y-axis
    main = "Drift Bursts test statistic"
    ylab = "test-statistic"
    plot(tstat, x = xtext, type = "l", xaxt = 'n', ylab = ylab, main = main, xlab = xlab, xlim = xlim)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = horizLines, col = "grey" , lty = 3, cex = 0.1)
    legend.txt = "t-stat"
    if(!is.null(price)){
      if(is.null(timestamps)){
        stop("The timestamps of the price must be passed in the timestamps argument")
      }
      par(new = TRUE)
      plot(price, x = timestamps , type = "l", axes = FALSE, col = "red", xlab = "", ylab = "", lty = 2, xlim = xlim)  
      axis(4)
      mtext(side = 4, text = "price", line = 2.5)
      legend.txt = c(legend.txt, "price")
      legend(x = leg.x, leg.y, legend = legend.txt, lty = c(1,2), col = c(1,2), bg = rgb(0,0,0,0), box.lwd = 0, 
             box.col = rgb(0,0,0,0))
    }
  }
  if(all(which == "sigma")){ #use all() because this function should accept which arguments with length longer than 1
    main = "volatility"
    ylab = "local volatility"
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    plot(sigma, x = xtext, type = "l",  xaxt = 'n', ylab = ylab, main = main, xlab = xlab)  
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
  }
  if(all(which == "mu")){ #use all() because this function should accept which arguments with length longer than 1
    main = "drift"
    ylab = "drift"
    if(annualize){ 
      ylab = paste("annualized", ylab)
    }
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    plot(mu, x = xtext, type = "l",  xaxt = 'n', ylab = ylab, main = main, xlab = xlab)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = 0, col = "grey" , lty = 3)
  }
  if("mu" %in% which & "sigma" %in% which){
    par(mfrow = c(2,1), omi = c(0,0,0,0), mgp = c(2,1,0), mai = c(0.75,0.75,0.3,0.25))
    main = "drift"
    ylab = "drift"
    if(annualize){ 
      ylab = paste("annualized", ylab)
    }
    plot(mu, x = xtext, type = "l", xlab = "",  xaxt = 'n', ylab = ylab, main = main)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = 0, col = "grey" , lty = 3)
    main = "volatility"
    ylab = "volatility"
    if(annualize){ 
      ylab = paste("annualized", ylab)
    }
    plot(sigma, x = xtext, type = "l", xlab = "", xaxt = 'n', ylab = ylab, main = main) 
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
  }
  par(startpar)
}

print.DBH = function(x, ...){
  usePolynomialInterpolation = TRUE
  options = list(...)
  if('criticalValue' %in% names(options)){
    usePolynomialInterpolation = FALSE
  }
  #### List of standard options
  opt = list(alpha = 0.95)
  #### Override standard options where user passed new options
  opt[names(options)] = options
  if(usePolynomialInterpolation){
    alpha = opt$alpha
    criticalValue = getCriticalValues(x, alpha)$quantile  
  }else{
    criticalValue = opt$criticalValue
  }
  
  varDB = getVar(x, which = 'db')
  padding = x$info$padding
  #We always remove the first entry, as this is used to denote the start of trading.
  whichToInclude = seq(padding[1] + 1, length(x$info$testTimes)- padding[2]) 
  allMeans = getMean(x, which = 'all')
  cat("\n-------------Drift Burst Hypothesis------------\n")
  cat("Tests performed:                     ", length(whichToInclude))
  #browser()
  if(usePolynomialInterpolation){
    cat("\nAny drift bursts (|T| > ",paste0(round(criticalValue,3)),"):    ", ifelse(any(abs(getDB(x))>criticalValue) , 'yes', 'no'))
  }else{
    cat("\nAny drift bursts (|T| > ",paste0(criticalValue[1]),"):        ", ifelse(any(abs(getDB(x))>criticalValue) , 'yes', 'no'))
  }
  cat("\nMax absolute value of test statistic:", round(max(abs(getDB(x))), digits=5))
  cat("\nMean test statistic:                 ", round(allMeans$meanDB, digits = 5))
  cat("\nVariance of test statistic:          ", round(varDB, digits = 5))
  cat("\nMean mu:                             ", round(allMeans$meanMu, digits = 5))
  cat("\nMean sigma:                          ", round(allMeans$meanSigma, digits = 5))
  cat("\n-----------------------------------------------\n")
}

getDB = function(x){
  UseMethod("getDB", x)
}

getDB.DBH = function(x){
  DB = x$driftBursts
  return(DB)
}


getSigma = function(x, annualize = FALSE, nDays = 252){
  UseMethod("getSigma", x)
}

getSigma.DBH = function(x, annualize = FALSE, nDays = 252){
  sigma = sqrt((x$sigma * 2 * x$info$nObs)  / (x$info$nObs / (x$info$sessionEnd - x$info$sessionStart)))/
    (x$info$preAverage^2)
  if(annualize){sigma = sigma * sqrt(nDays)}
  return(sigma)
}


getMu = function(x, annualize = FALSE, nDays = 252){
  UseMethod("getMu", x)
}

getMu.DBH = function(x, annualize = FALSE, nDays = 252){
  mu = (x$mu * x$info$meanBandwidth / (x$info$nObs / (x$info$sessionEnd - x$info$sessionStart)))/
    (x$info$preAverage^2 * 2)
  if(annualize){mu =  mu * nDays}
  return(mu)
}

getMean = function(x, which = 'all'){
  UseMethod('getMean', x)
}

getMean.DBH = function(x, which = 'all'){
  which     = tolower(which)
  padding = x$info$padding
  if(!(which %in% c('all', 'db', 'driftbursts', 'mu', 'sigma'))){
    stop("The which argument must be a character vector containing either:\n
         Sigma, Mu, both of these or DriftBursts. 
         CasE doesn't matter.")
  }
  #We always remove the first entry, as this is used to denote the start of trading.
  whichToInclude = seq(padding[1] + 1, length(x$info$testTimes)- padding[2]) 
  
  if(which == 'all'){
    meanDB = mean(getDB(x)[whichToInclude], na.rm = TRUE)
    meanMu = mean(getMu(x)[whichToInclude], na.rm = TRUE)
    meanSigma = mean(getSigma(x)[whichToInclude], na.rm = TRUE)  
    out = list('meanDB' = meanDB, 'meanMu' = meanMu, 'meanSigma' = meanSigma)
  }
  if(which %in% c('driftbursts', 'db')){
    out = mean(getDB(x)[whichToInclude], na.rm = TRUE)
  }
  if(which == 'mu'){
    out = mean(getMu(x)[whichToInclude], na.rm = TRUE)
  }
  if(which == 'sigma'){
    out = mean(getSigma(x)[whichToInclude], na.rm = TRUE)
  }
  
  return(out)
  
}


getVar = function(x, which = 'all', annualize = FALSE, nDays = 252){
  UseMethod('getVar', x)
}


getVar.DBH = function(x, which = 'all', annualize = FALSE, nDays = 252){
  which     = tolower(which)
  padding = x$info$padding
  if(!(which %in% c('all', 'db', 'driftbursts', 'mu', 'sigma'))){
    stop("The which argument must be a character vector containing either:\n
         Sigma, Mu, both of these or DriftBursts. 
         CasE doesn't matter.")
  }
  #We always remove the first entry, as this is used to denote the start of trading.
  whichToInclude = seq(padding[1] + 1, length(x$info$testTimes)- padding[2]) 
  if(which == 'all'){
    varDB = var(getDB(x)[whichToInclude])
    varMu = var(getMu(x, annualize, nDays)[whichToInclude])
    varSigma = var(getSigma(x, annualize, nDays)[whichToInclude])  
    out = list('varDB' = varDB, 'varMu' = varMu, 'varSigma' = varSigma)
  }
  if(which %in% c('driftbursts', 'db')){
    out = var(getDB(x)[whichToInclude])
  }
  if(which == 'mu'){
    out = var(getMu(x, annualize, nDays)[whichToInclude])
  }
  if(which == 'sigma'){
    out = var(getSigma(x, annualize, nDays)[whichToInclude])
  }
  
  return(out)

}



getCriticalValues = function(x, alpha = 0.95){
  UseMethod('getCriticalValues', x)
}

getCriticalValues.DBH = function(x, alpha = 0.95){
  return(DBHCriticalValues(x, alpha))
}