# Other comparisons: multiple metrics for the same options - better to do this outside the function because each metric could have its own colLim & thresholds
# So it is better to retain 'plotOptions' as a base function that operates on a single metric
# if functionality to compare multiple performance metrics are required it is better to have a wrapper function of plotOptions (plotTwoOptions?)
# function plots performanceOpt2 minus performanceOpt1 and contours showing the movement of performance threshold contours

#' Plots the differences in performance metrics from two system options
#'
#' \code{plotOptions} uses the system model performances calculated using the function \code{runSystemModel} for two alternate system model options,
#' and the summary of the simulation generated using the functions \code{generateScenarios} & \code{getSimSummary} as input. The function plots the differences in 
#' the performance metrics between the two options, and the changes in performance thresholds in the space. 
#' The user may specify the attributes to be used as the axes of the plot. The function contains arguments to control the finer details of the plot.
#' @inheritParams plotPerformanceSpace
#' @param performanceOpt1 a named list; contains the system model performance calculated using \code{runSystemModel} for system model option 1. 
#' If the list contains more than one performance metric, the argument \code{metric} can be used to specify the metric to be used.
#' @param performanceOpt2 a named list; contains the system model performance calculated using \code{runSystemModel} for system model option 2. 
#' If the list contains more than one performance metric, the argument \code{metric} can be used to specify the metric to be used.
#' @param metric a string; the name of the performance metric to be plotted. The argument can be used to select the metric from 
#' \code{performanceOpt1} and \code{performanceOpt2} lists for plotting. If \code{NULL} (the default), the first metric in the lists will be used. 
#' @param opt1Label a string; the text to label \code{performanceOpt1}.
#' @param opt2Label a string; the text to label \code{performanceOpt2}.
#' @param titleText a string; text for the title of the plot. The default is \code{paste0(opt2Label, " - ", opt1Label)}.
#' @return The plot of the differences in the performance metrics (option 2 - option 1) in a ggplot object.
#' @seealso \code{runSystemModel}, \code{plotPerformanceSpace}, \code{generateScenarios}, \code{getSimSummary} 
#' @examples
#' # load example datasets
#' data("egSimSummary")
#' data("egSimPerformance")       # performance of option1
#' data("egSimPerformance_systemB")  # performance of option2
#' data("egClimData")
#' plotOptions(egSimPerformance[1], egSimPerformance_systemB [1], egSimSummary, 
#' attX = "P_ann_seasRatio_m", attY = "P_ann_tot_m", topReps = 7, perfThreshLabel = "Threshold (28L)",
#' perfThresh = 28,  opt1Label = "System A", opt2Label = "System B", climData = egClimData)
#' @export

plotOptions <- function(performanceOpt1,               # system model performance for option 1, matrix of size targets x replicates
                        performanceOpt2,               # system model performance for option 2, matrix of size targets x replicates
                        sim,                           # simulation containing all targets and replicates
                        metric = NULL,                 # name of the performance metric
                        attX = NULL,                   # attribute to be plotted on x-axis
                        attY = NULL,                   # attribute to be plotted on y-axis
                        topReps = NULL,                # number of topReps based on fitness (i.e., the objective function score which is -ve)
                        opt1Label = "Option 1",        # label of system option 1
                        opt2Label = "Option 2",        # label of system option 2
                        titleText = paste0(opt2Label, " - ", opt1Label), # title of the plot
                        perfThresh = NULL,             # desired performance threshold; plot would contain a contour to mark this threshold
                        perfThreshLabel = "Threshold", # label text for the threshold
                        attSlices = NULL,              # list containing the slices of attributes to use for plotting
                        climData = NULL,               # changes in climate attributes from other sources - can include label. If the performance measure being plotted is a column in the data.frame, the points will be coloured accordingly
                        colMap = NULL,                 # alternate colormap
                        colLim = NULL                  # if null, the full limit is used
) {
  
  # assuming that performance is a list with a name
  # it may also be a matrix without a name; will be named "performance"
  if (is.list(performanceOpt1) & is.list(performanceOpt2)) {
    if (!is.null(metric)) {
      if (!is.character(metric)) stop("`metric` should be specified as a string.")
      if (length(metric) > 1) stop("A single `metric` should be specified.")
      performanceOpt1 <- performanceOpt1[metric]
      performanceOpt2 <- performanceOpt2[metric]
      if (is.null(performanceOpt1[[1]])) stop(paste0("Cannot find metric `", metric, "` in performanceOpt1."))
      if (is.null(performanceOpt2[[1]])) stop(paste0("Cannot find metric `", metric, "` in performanceOpt2."))
      perfName <- metric
    } else {
      if (is.null(names(performanceOpt1)) & is.null(names(performanceOpt2))) {
        perfName <- rep("Performance", length(performanceOpt1))
      } else {
        if (!identical(names(performanceOpt1), names(performanceOpt2))) stop("The metrics in performanceOpt1 and performanceOpt2 should be the same.")
        perfName <- names(performanceOpt1)
      }
    }
  } else if (is.matrix(performanceOpt1) & is.matrix(performanceOpt2)) {
    # already a matrix
    perfName <- "Performance"
    # change to list
    perfMat <- performanceOpt1
    performanceOpt1 <- list()
    performanceOpt1[[1]] <- perfMat; rm(perfMat)
    perfMat <- performanceOpt2
    performanceOpt2 <- list()
    performanceOpt2[[1]] <- perfMat; rm(perfMat)
  } else {
    stop("performanceOpt1 and performaceOpt2 should be list or matrix, and should be of the same class.")
  }
  
  # unpacking sim metadata
  repNames <- names(sim[grep("Rep", names(sim))])
  if (is.null(repNames)) stop("There are no replicates in sim.")
  tarNames <- names(sim[[repNames[1]]])
  nRep <- length(repNames)
  nTar <- length(tarNames)
  targetMat <- sim[["expSpace"]][["targetMat"]]
  
  if (is.list(sim[["controlFile"]])) {
    # get the simulation fitness if stochastic
    simFitness <- getSimFitness(sim)
    
    # check that the performance is calculated from this simulation
    if (!identical(dim(simFitness), dim(performanceOpt1[[1]]))) {
      stop("The number of targets and replicates in sim and performanceOpt1 should match. Is the performanceOpt1 calculated using sim?")
    }
    if (!identical(dim(simFitness), dim(performanceOpt2[[1]]))) {
      stop("The number of targets and replicates in sim and performanceOpt2 should match. Is the performanceOpt2 calculated using sim?")
    }
    if (!is.null(topReps)) {
      if (topReps > nRep) {
        message(paste0("sim does not contain ", topReps, " replicates. Using all replicates available in sim.\n"))
      }
    }
  } else {
    if (!(sim[["controlFile"]] == "scaling")) stop("sim$controlFile is unrecognized.")
    # scaling
    simFitness <- NULL
    topReps <- NULL
  }
  
  # subset targets if attSlices are specified - TO DO - change attSlices to specify the number of the slice rather than the value
  if (!is.null(attSlices)) {
    # indices to subset the rows of the targetMat 
    tarInd <- getSliceIndices(sim[["expSpace"]], attSlices)
    # subset targetMat
    targetMat <- targetMat[tarInd, ]
  } else {
    tarInd <- NULL
  }
  
  # if null get appropriate tags
  if (is.null(attX) | is.null(attY)) {
    attPerturb <- sim[["expSpace"]][["attPerturb"]]
    attList <- getAttXY(attPerturb, attX, attY)
    attX <- attList[["attX"]]
    attY <- attList[["attY"]]
  }
  
  # check the perturbations in all attributes
  checkAttPert(targetMat, attX, attY)
  
  # identify x and y columns
  attNames <- colnames(targetMat)
  colAttX <- which(attNames == attX)
  colAttY <- which(attNames == attY)
  
  nPerf1 <- length(performanceOpt1)
  nPerf2 <- length(performanceOpt2)
  perfPlots <- list()
  
  if (!(nPerf1 == nPerf2)) message("The number of metrics in performanceOpt1 and performanceOpt2 are different. Using the first metric from both.")
  if (nPerf1 > 1 | nPerf2 > 1) message("performance options contain more than one performance metric, the first metric is plotted.")
  
  # loop over multiple metrics - can do this, but becomes complicated due for multiple colMap, colLim, and threshold labels
  #for (p in 1) {
  # subset performance based on slices
  if (!is.null(tarInd)) {
    perfMatrixOpt1 <- performanceOpt1[[1]][tarInd, ]
    perfMatrixOpt2 <- performanceOpt2[[1]][tarInd, ]
  } else {
    perfMatrixOpt1 <- performanceOpt1[[1]]
    perfMatrixOpt2 <- performanceOpt2[[1]]
  }
  
  # get the average performance & difference to be plotted, name appropriately
  performanceAvOpt1 <- getPerfStat(perfMatrixOpt1, simFitness, topReps, nRep)
  performanceAvOpt2 <- getPerfStat(perfMatrixOpt2, simFitness, topReps, nRep)
  performanceAvDiff <- performanceAvOpt2 - performanceAvOpt1
  names(performanceAvOpt1) <- perfName[1]
  names(performanceAvOpt2) <- paste0(opt2Label, perfName[1])
  names(performanceAvDiff) <- perfName[1]
  
  # create data.frame for plotting
  perfPlotData <- cbind(targetMat[colAttX], targetMat[colAttY], performanceAvDiff)
  
  if (is.null(colMap)) {
    # set to the default divergent colourmap if there are positive and negative values in the difference matrix
    if ((sum(performanceAvDiff > 0) > 0) & (sum(performanceAvDiff < 0) > 0)) {
      colMap = foreSIGHT.divColmap
      if (is.null(colLim)) colLim <- c((-1*max(abs(performanceAvDiff))), max(abs(performanceAvDiff)))
    }
  }
  # if the user specified divergent limits with colMap set to NULL
  if (!(is.null(colLim))) {
    if (colLim[1] < 0 & colLim[2] > 0) {
      if (is.null(colMap)) colMap <- foreSIGHT.divColmap
    }
  }
  
  # Base R: not used
  # plotPerfQuiltPlot(perfPlotData, nx, ny, colLim = colLim, colBar = colBar, perfThresh = perfThresh, climData = climData)
  # plot performance threshold difference heatmap without thershold contours - because thresholds are based on original metric, not the difference
  p1 <- heatPlot(perfPlotData, colLim = colLim, colMap = colMap, 
                        perfThresh = NULL, perfThreshLabel = perfThreshLabel, climData = climData)
  
  # Add title
  p1 <- p1 + ggplot2::ggtitle(titleText)
  
  if (!is.null(perfThresh)) {
    
    perfPlotData <- cbind(targetMat[colAttX], targetMat[colAttY], performanceAvOpt1)
    opt1Mean <- getAggMean(perfPlotData)
    p1 <- addThresholdLines(p1, opt1Mean, perfThresh, perfThreshLabel = paste0(opt1Label, " ", perfThreshLabel), names(opt1Mean)[1:2], names(opt1Mean)[3], 
                            lineCol = option1_threshCol, lineSize = option1_lineSize, lineAlpha = option1_lineAlpha, label = "end")
    
    # Add the thresholds for opt1 and opt2
    perfPlotData <- cbind(targetMat[colAttX], targetMat[colAttY], performanceAvOpt2)
    opt2Mean <- getAggMean(perfPlotData)
    p1 <- addThresholdLines(p1, opt2Mean, perfThresh, perfThreshLabel = paste0(opt2Label, " ", perfThreshLabel), names(opt2Mean)[1:2], names(opt2Mean)[3], 
                            lineCol = option2_threshCol, lineSize = option2_lineSize, lineAlpha = option2_lineAlpha)
  }
  #print(perfPlots[[1]])
  #}
  print(p1)
  return(invisible(p1))
  
}

getAggMean <- function(data) {
  # aggregate data
  tempMat <- data
  names(tempMat) <- c("x", "y", "z")
  dataMean <- aggregate(.~x+y, tempMat, mean)
  names(dataMean) <- names(data)
  return(dataMean)
}

# Function to add thresholds based on the original performance over the difference heat maps created by setting performance thresholds to NULL


