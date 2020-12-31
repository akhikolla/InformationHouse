## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  python.reticulate = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  calculateCounts = function(data=c(), parameters=c()) {
#    # data: 3 column matrix with acc data
#    # parameters: the sample rate of data
#    library("activityCounts")
#    if (ncol(data) == 4) data= data[,2:4]
#    mycounts = counts(data=data, hertz=parameters,
#                      x_axis=1, y_axis=2, z_axis=3,
#                      start_time = Sys.time())
#    mycounts = mycounts[,2:4] #Note: do not provide timestamps to GGIR
#    return(mycounts)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  source("~/calculateCounts.R")
#  myfun =  list(FUN=calculateCounts,
#                parameters= 30,
#                expected_sample_rate= 30,
#                expected_unit="g",
#                colnames = c("countsX","countsY","countsZ"),
#                outputres = 1,
#                minlength = 1,
#                outputtype="numeric",
#                aggfunction = sum,
#                timestamp=F,
#                reporttype="scalar")

## ----eval=FALSE---------------------------------------------------------------
#  library(GGIR)
#  g.shell.GGIR(datadir="~/myaccelerometerdata",
#               outputdir="~/myresults",
#               mode=1:2,
#               epochvalues2csv = TRUE,
#               do.report=2,
#               myfun=myfun) #<= this is where object myfun is provided to g.shell.GGIR

## ----eval=FALSE---------------------------------------------------------------
#  dominant_frequency = function(data=c(), parameters=c()) {
#    # data: 3 column matrix with acc data
#    # parameters: the sample rate of data
#    source_python("dominant_frequency.py")
#    sf=parameters
#    N = nrow(data)
#    ws = 5 # windowsize
#    if (ncol(data) == 4) data= data[,2:4]
#    data = data.frame(t= floor(seq(0,(N-1)/sf,by=1/sf)/ws),
#                      x=data[,1], y=data[,2], z=data[,3])
#    df = aggregate(data, by = list(data$t),
#                   FUN=function(x) {return(dominant_frequency(x,sf))})
#    df = df[,-c(1:2)]
#    return(df)
#  }
#  }

## ----eval=FALSE---------------------------------------------------------------
#    library("reticulate")
#    use_virtualenv("~/myvenv", required = TRUE) # Local Python environment
#    py_install("numpy", pip = TRUE)
#  

## ----eval=FALSE---------------------------------------------------------------
#  source("~/dominant_frequency.R")
#  myfun =  list(FUN=dominant_frequency,
#                parameters= 30,
#                expected_sample_rate= 30,
#                expected_unit="g",
#                colnames = c("domfreqX", "domfreqY", "domfreqZ"),
#                minlength = 5,
#                outputres = 5,
#                outputtype="numeric",
#                aggfunction = median
#                timestamp=F,
#                reporttype="scalar")

## ----eval=FALSE---------------------------------------------------------------
#  library(GGIR)
#  g.shell.GGIR(datadir="~/myaccelerometerdata",
#               outputdir="~/myresults",
#               mode=1:2,
#               epochvalues2csv = TRUE,
#               do.report=2,
#               myfun=myfun,
#               do.parallel = FALSE)

