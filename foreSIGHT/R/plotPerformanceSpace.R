# POSSIBLE NAMES:
# plotPerformanceSpace
# plotSystemPerformace

# The other function for OAT performance plots can be:
# plotPerformanceOAT

#------------------------------------------------------------------------------------------
#Plotting grid

# visualizeSpaces()
#------------------------------------------------------------------------------------------

# User input that will become function arguments

# should the plotScenarios function also include topReps? so that the user knows what is the mean bias and SD associated with those replicates
# should the plot include a colorbar? For a wrapper function to plot multiple panels, this could be false


getAttXY <- function(attPerturb, attX, attY) {
  
  if (length(attPerturb) < 2) stop(paste0("sim should contain two or more perturbed attributes to plot a performance space."))
  
  if (is.null(attX) & is.null(attY)) {
    attX <- attPerturb[1]
    attY <- attPerturb[2]
    message(paste0("Using attX = ", attX, ", and attY = ", attY, ". Specify attX and attY in the function call to choose alternate attributes."))
    
  } else {
    
    attInd <- which(attPerturb %in% c(attX, attY))
    if (!identical(attInd, integer(0))) {
      attPerturb <- attPerturb[-attInd]
    }
    if (is.null(attX)) {
      attX <- attPerturb[1]
      message(paste0("Using attX = ", attX, ". Specify attX in the function call to choose an alternate attribute."))
    } else {
      attY <- attPerturb[1]
      message(paste0("Using attY = ", attY, ". Specify attY in the function call to choose an alternate attribute."))
    }
  }
  
  return(list(attX = attX,
              attY = attY))
  
}




# returns the sum of fitness as a matrix
getSimFitness <- function(sim) {
  
  # unpacking sim metadata
  repNames <- names(sim[grep("Rep", names(sim))])
  if (is.null(repNames)) stop("There are no replicates in sim.")
  tarNames <- names(sim[[repNames[1]]])
  
  nRep <- length(repNames)
  nTar <- length(tarNames)
  
  # get the fitness value
  fitName <- "score"
  simFitness <- matrix(NA, nrow = length(tarNames), ncol = length(repNames)) 
  for (r in 1:nRep) {
    simFitness[ ,r] <- sapply(sim[[repNames[r]]], function(x){sum(x[["score"]])})
  }
  return(simFitness)
}


getPerfStat <- function(performance, 
                      simFitness, 
                      topReps,
                      nRep,
                      statFUN = mean
                      ) {
  nTar <- nrow(performance)
  # rank in order of fitness (take top X)
  # this ordering is done for each target
  if (!is.null(topReps)) {
    if (topReps < nRep) {
      sortedPerf <- matrix(0, nrow = dim(performance)[1], ncol = dim(performance)[2])
      for(i in 1:nTar){
        ind <- order(simFitness[i, ], decreasing = TRUE)  # make the best performance 1st
        sortedPerf[i, ] <- performance[i, ind]
      }
      # performance as the average of nReps (in terms of scenario fit) of the set (USER INPUT)
      performanceStat <- as.data.frame(apply(sortedPerf[ , 1:topReps], 1, FUN = statFUN))
    } else {
      performanceStat <- as.data.frame(apply(performance, 1, FUN = statFUN))
    }
  } else {
    performanceStat <- as.data.frame(apply(performance, 1, FUN = statFUN))
  }
  return(performanceStat)
}

checkAttPert <- function(targetMat, attX, attY) {
  
  if (attX == attY) {
    stop("attY should be different from attX.")
  }
  
  # Check attX and attY
  attNames <- colnames(targetMat)
  colAttX <- which(attNames == attX)
  colAttY <- which(attNames == attY)
  
  if (length(colAttX) == 0 | colAttY == 0) {
    stop("attX or attY is not available in sim.")
  }
  
  # tarAttX <- targetMat[ ,colAttX]
  # tarAttY <- targetMat[ ,colAttY]
  
  # May not need two attributes - can plot a single tile using just one sample?
  # if (length(unique(tarAttX)) < 2) stop("Perturbation of attX should contain atleast two samples to plot the performance space.")
  # if (length(unique(tarAttY)) < 2) stop("Perturbation of attY should contain atleast two samples to plot the performance space.")
  
  # Check the rest of the attributes - commented since it does not work for some reason
  # tarAttRest <- targetMat[ , -c(colAttX, colAttY)]
  # for (a in 1:ncol(tarAttRest)) {
  #   allPointsAtt <- sort(unique(tarAttRest[ ,a]))
  #   if (!(all( abs(allPointsAtt - mean(allPointsAtt)) < 0.000005 ))) {
  #     aName <- colnames(tarAttRest)[a]
  #     warning(paste0(aName, " is not perturbed on a regular grid. The performance space will be averaged across the irregular perturbations in this attribute."))
  #   }
  # }
  return(invisible())
}

# function to return the indices to subset the exposure space
# will be used to subset targetMat and performance
getSliceIndices <- function(expSpace, attSlices) {
  
  if (!(is.list(attSlices))) stop("attSlices should be a list.")
  attIn <- names(attSlices)
  if (is.null(attIn)) stop("attSlices should be named using the attribute tags to be sliced.")
  
  if (sum(attIn %in% c(expSpace[["attPerturb"]], expSpace[["attHold"]])) != length(attIn)) stop("attSlices should contain only attributes present in the expSpace of sim.")
  
  # identify whether perturbed or held
  attPInd <- which(attIn %in% expSpace[["attPerturb"]])
  attHInd <- which(attIn %in% expSpace[["attHold"]])
  
  # slices on attHold can only be 0(for Temp) or 1(other variables)
  # user is not expected to input these, in case they do input it the following code ensures that it does not result in an error.
  if (!identical(attHInd, integer(0))) {
    for (i in 1:length(attHInd)) {
      iAtt <- attIn[attHInd[i]]
      iAttVar <- strsplit(iAtt, "_")[[1]][1]
      if (iAttVar == "Temp") {
        if (attSlices[[iAtt]] != 0) stop(paste0(iAtt, " is a held attribute. attSlices[[\"", iAtt, "\"]] should be 0 or NULL."))
      } else {
        if (attSlices[[iAtt]] != 1) stop(paste0(iAtt, " is a held attribute. attSlices[[\"", iAtt, "\"]] should be 1 or NULL."))
      }
    }
  }
  
  targetMat <- expSpace[["targetMat"]]
  
  if (!identical(attPInd, integer(0))) {
    # perturbed attributes, this is where the actual slicing happens
    sliceIndicesList <- list()
    for (i in 1:length(attPInd)) {
      iAtt <- attIn[attPInd[i]]
      iSlice <- attSlices[[iAtt]]
      # single value to slice
      if (length(iSlice) == 1) {
        iSliceInd <- which(targetMat[[iAtt]] == iSlice)
      # a range to slice
      } else if (length(iSlice == 2)) {
        sMin <- min(iSlice)
        sMax <- max(iSlice)
        iSliceInd <- which(targetMat[[iAtt]]>=sMin & targetMat[[iAtt]]<=sMax)
      }
      if (identical(iSliceInd, integer(0))) stop(paste0("attSlices[[\"", iAtt, "\"]] is not valid for this expSpace."))
      sliceIndicesList[[i]] <- iSliceInd
    }
    sliceIndices <-  Reduce(intersect, sliceIndicesList)
  } else {
    sliceIndices <- c(1:nrow(targetMat))
  }
  return(sliceIndices)
}

#' Plots a performance space using the system performance and scenarios as input
#'
#' \code{plotPerformanceSpace} uses the system model performance calculated using the function \code{runSystemModel} and 
#' the summary of the simulation generated using the functions \code{generateScenarios} & \code{getSimSummary} as input to plot the performance space of the system.
#' The user may specify the attributes to be used as the axes of the performance space.
#' @param performance a named list; contains the system model performance calculated using \code{runSystemModel}. 
#' If the list contains more than one performance metric, the argument \code{metric} can be used to specify the metric to be used. 
#' @param sim a list; summary of the simulation containing the scenarios generated using the function \code{generateScenarios} that is used to run the system model using \code{runSystemModel}.
#' The summary of the simulation may be obtained by using the function \code{getSimSummary} on the full simulation. The summary object is much smaller in size for ease of storage and use with the performance 
#' plotting functions like \code{plotPerformanceSpace}.
#' @param metric a string; the name of the performance metric to be plotted. The argument can be used to select a metric from \code{performance} for plotting. 
#' If \code{NULL} (the default), the first metric in the list will be used.
#' @param attX a string; the tag of the perturbed attribute to plot on the xaxis. The attribute must be one of the perturbed attributes of \code{sim}. 
#' Type \code{sim$expSpace$attPerturb} to view all perturbed attributes of \code{sim}. If \code{NULL} (default), the first perturbed attribute of \code{sim} will be used.
#' @param attY a string; the tag of the perturbed attribute to plot on the yaxis. The attribute must be another perturbed attribute of \code{sim}.
#' If \code{NULL}, the second perturbed attribute of \code{sim} will be used.
#' @param topReps an integer (default is \code{NULL}); the number of "top" replicates in terms of simulation fitness to be used. If \code{topReps} is specified, \code{topReps}
#' number of replicates will be identified for each target and the average performance across these replicates will be plotted. If \code{NULL}, the average performance across all the replicates will be plotted.
#' @param perfThresh a number; the minimum or maximum threshold value of the performance metric. A line will be drawn to mark this threshold value in the performance space.
#' @param perfThreshLabel a string; the text to label \code{perfThresh}.
#' @param attSlices a list; used to subset perturbed attributes in \code{sim} for the plot. This argument would typically be used in cases where there are more than two perturbed attributes.
#' The elements of the list correspond to the perturbed attributes to be subsetted and must be named using the attribute tag. Each element may contain a single value or a two-element vector specifying the minimum-maximum values. 
#' If the element is a single value, the exposure space is sliced on this single value of the attribute. If minimum-maximum values are specified, the exposure space will be sliced to subset this range.
#' If \code{attSlices} includes \code{attX} or \code{attY}, these attributes will be sliced and the resulting plot will be a "zoomed-in" space.
#' @param climData data.frame; the values of attX and attY from other sources like climate models. This data will be plotted as points in the performance space.
#' The data frame may contain columns with values of the performance metric to be plotted and the "Name" of the dataset. 
#' If the performance metric is available in the data.frame, the points will be coloured based on the performance \code{colMap} scale.
#' If the \code{Name} of the data is available in the data.frame, the points will be identified using the \code{Name}. 
#' Please refer data provided with the package that may be loaded using \code{data("egClimData")} for an example of the expected format of \code{climData}.
#' @param colMap a vector of colours; to specify the colourmap to be used. If \code{NULL}, the default foreSIGHT colourmap is used.
#' @param colLim a vector of 2 values; the minimum and maximum limits of the colour scale.
#' @details If the space contains more than two perturbed attributes, the performance values are averaged across the perturbations in the attributes other than \code{attX} and \code{attY}.
#' The user may specify argument \code{attSlices} to slice the performance space at specific values of the other perturbed attributes. If \code{attSlices} are used to 
#' specify minimum-maximum values to subset other perturbed attributes, the performance values are averaged across the subsetted perturbations in these attributes.
#' If the input performance list contains multiple performance metrics, the function plots the first metric.
#' The function may be called with \code{performance} argument specifying the metric to be plotted \code{plotPerformanceSpace(performance[2], sim)} to plot other metrics.
#' @return The plot of the performance space and the ggplot object.
#' @seealso \code{runSystemModel}, \code{generateScenarios}, \code{getSimSummary}, \code{plotPerformanceOAT}
#' @examples
#' # load example datasets
#' data("egSimSummary")       # summary of stochastic simulation
#' data("egSimPerformance")   # system performance calculated using the stochastic simulation
#' data("egClimData")         # alternate climate data and system performance
#'
#' plotPerformanceSpace(egSimPerformance[2], egSimSummary)
#' 
#' # adding climate data, using top 10 replicates
#' plotPerformanceSpace(egSimPerformance[1], egSimSummary, topReps = 10, climData = egClimData)
#' 
#' # adding a threshold
#' plotPerformanceSpace(egSimPerformance, egSimSummary, metric = "Avg. Deficit (L)", 
#' climData = egClimData, perfThresh = 27.5, perfThreshLabel = "Max Avg. Deficit")
#' 
#' # user specified colMap
#' plotPerformanceSpace(egSimPerformance[1], egSimSummary, climData = egClimData, perfThresh = 27.5, 
#' perfThreshLabel = "Max Avg. Deficit", colMap = viridisLite::inferno(100))
#' 
#' # example of performance generated using simple scaled simulation
#' data("egScalPerformance")
#' data("egScalSummary")
#' data("egClimData")
#' plotPerformanceSpace(egScalPerformance[1], egScalSummary, climData = egClimData, 
#' perfThresh = 28.25, perfThreshLabel = "Max Avg. Deficit")
#' @export

plotPerformanceSpace <- function(performance,                   # system model performance, matrix of size targets x replicates 
                                 sim,                           # simulation containing all targets and replicates
                                 metric = NULL,                 # name of the performance metric
                                 attX = NULL,                   # attribute to be plotted on x-axis
                                 attY = NULL,                   # attribute to be plotted on y-axis
                                 topReps = NULL,                # number of topReps based on fitness (i.e., the objective function score which is -ve)
                                 perfThresh = NULL,             # desired performance threshold; plot would contain a contour to mark this threshold
                                 perfThreshLabel = "Threshold", # label text for the threshold
                                 attSlices = NULL,              # list containing the slices of attributes to use for plotting
                                 climData = NULL,               # changes in climate attributes from other sources - can include label. If the performance measure being plotted is a column in the data.frame, the points will be coloured accordingly
                                 colMap = NULL,                 # alternate colormap
                                 colLim = NULL                  # if null, the full limit is used
                                 ) {
  
  # assuming that performance is a list with a name
  # it may also be a matrix without a name; will be named "performance"
  if (is.list(performance)) {
    if (!is.null(metric)) {
      if (!is.character(metric)) stop("`metric` should be specified as a string.")
      if (length(metric) > 1) stop("A single `metric` should be specified.")
      performance <- performance[metric]
      if (is.null(performance[[1]])) stop(paste0("Cannot find metric `", metric, "` in performance."))
      perfName <- metric
    } else {
      if (is.null(names(performance))) {
        perfName <- rep("Performance", length(performance))
      } else {
        perfName <- names(performance)
      }
    }
  } else if (is.matrix(performance)) {
    # already a matrix
    perfName <- "Performance"
    # change to list
    perfMat <- performance
    performance <- list()
    performance[[1]] <- perfMat
    rm(perfMat)
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
    if (!identical(dim(simFitness), dim(performance[[1]]))) {
      stop("The number of targets and replicates in sim and performance should match. Is the performance calculated using sim?")
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
  
  # subset targets if attSlices are specified
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
  
  nPerf <- length(performance)
  perfPlots <- list()
  
  if (nPerf > 1) message("performance contains more than one performance metric, the first metric is plotted.")
  
  # loop over multiple metrics - can do this, but becomes complicated due for multiple colMap, colLim, and threshold labels
  #for (p in 1) {
    # subset performance based on slices
    if (!is.null(tarInd)) {
      perfMatrix <- performance[[1]][tarInd, ]
    } else {
      perfMatrix <- performance[[1]]
    }
    
    # get the average performance to be plotted
    performanceAv <- getPerfStat(perfMatrix, simFitness, topReps, nRep)
    
    # name appropriately
    names(performanceAv) <- perfName[1]
    
    # create data.frame for plotting
    perfPlotData <- cbind(targetMat[colAttX], targetMat[colAttY], performanceAv)
    
    # Base R: not used
    # plotPerfQuiltPlot(perfPlotData, nx, ny, colLim = colLim, colBar = colBar, perfThresh = perfThresh, climData = climData)
    perfPlots <- heatPlot(perfPlotData, colLim = colLim, colMap = colMap, 
                   perfThresh = perfThresh, perfThreshLabel = perfThreshLabel, climData = climData)
    #print(perfPlots[[1]])
    
  #}
  print(perfPlots)
  return(invisible(perfPlots))

}


plotPerfQuiltPlot <- function(plotData, nx, ny, colLim = NULL, colBar = TRUE, perfThresh = NULL, climData = NULL){
  
  xyAtts <- colnames(plotData)[1:2]
  xyAttDefs <- mapply(tagBlender_noUnits, xyAtts, USE.NAMES = FALSE)
  
  if (is.null(colLim)) {
    tempMat <- plotData
    names(tempMat) <- c("x", "y", "z")
    tempMatAvg <- aggregate(.~x+y, tempMat, mean)
    colLim = c(min(tempMatAvg[ ,3]), max(tempMatAvg[ ,3]))
    rm(tempMat, tempMatAvg)
  }
  
  varNames <- sapply(strsplit(xyAtts, "_"), `[[`, 1)
  varUnits <- getVarUnits(varNames)
  xyLabels <- paste0(xyAttDefs, " (", varUnits, ")")
  
  
  par(mar=c(3.8,3.8,1,1),oma=c(5,0.3,0.3,0.3), mgp = c(0,0.5,0))
  fields::quilt.plot(x = plotData[ ,1], y = plotData[ ,2], z = plotData[ ,3], nx = nx ,ny = ny,
             add.legend = FALSE, nlevel = perfSpace_nlevel, col = foreSIGHT.colmap(perfSpace_nlevel), 
             xlim = c(min(plotData[ ,1]), max(plotData[ ,1])), ylim = c(min(plotData[ ,2]), max(plotData[ ,2])),
             zlim = colLim, xlab = "", ylab =)
  box(col = "black")

  mtext(side=1,text=xyLabels[1],line=1.8)
  mtext(side=2,text=xyLabels[2],line=1.8)

  if (perfSpace_contours) {
    # get image for contours
    look <- fields::as.image(plotData[ ,3], ind = cbind(plotData[ ,1], plotData[ ,2]), nx = nx, ny = ny)
    contour(add = TRUE, x = look$x, y = look$y, z = look$z, method="edge", labcex = 1, nlevels = perfSpace_nContour)
  }
  
  if (!is.null(perfThresh)) {
    contour(add = TRUE, x = look$x, y = look$y, z = look$z, lty=threshLty, lwd = threshLwd, col = perfSpace_threshCol, drawlabels = FALSE, levels = perfThresh)
    legend("topright", inset = c(0.005, 0.005), legend=paste0("Performance\nThreshold (", perfThresh, ")"),
           col=perfSpace_threshCol, lwd=threshLwd, lty=threshLty, cex=1, bty="n")
  }
  
  if (!is.null(climData)) {
    if (sum(colnames(climData) %in% xyAtts) == 2) {
      points(x = climData[[xyAtts[1]]], y = climData[[xyAtts[2]]], pch = perfSpace_climDataPch, col = perfSpace_climDataCol, bg = perfSpace_climDataBg)
    } else {
      warning(paste0("climData is not plotted since it does not contain ", paste(xyAtts[(colnames(climData) %in% xyAtts)], sep = ","), "."))
    }
  }
  
  if(colBar){
    fields::image.plot(legend.only = TRUE, zlim = colLim, col = foreSIGHT.colmap(perfSpace_nlevel), 
                horizontal = TRUE, smallplot=c(0.2,0.9,0.0001,0.02))
  }
  mtext(tag_text, side=1, line=1.5, adj=1.0, cex=0.8, col=tag_textCol, outer=TRUE)
  
}

heatPlot <- function(plotData, colLim = NULL, colMap = NULL, perfThresh = NULL, perfThreshLabel = "Threshold", climData = NULL) {
  
  level <- NULL
  
  xyAtts <- colnames(plotData)[1:2]
  perfName <- colnames(plotData)[3]
  
  xyAttDefs <- mapply(tagBlender_noUnits, xyAtts, USE.NAMES = FALSE)
  
  # aggregate data
  tempMat <- plotData
  names(tempMat) <- c("x", "y", "z")
  plotDataMean <- aggregate(.~x+y, tempMat, mean)
  names(plotDataMean) <- names(plotData)

  if (is.null(colLim)) {
    colLimIn = c(min(plotDataMean[ ,3]), max(plotDataMean[ ,3]))
  } else {
    colLimIn <- colLim
  }
  
  if (!is.null(perfThresh)) {
    tempData <- plotDataMean[ ,3]
    threshName <- paste0("Thresh", names(plotDataMean)[3])
    names(tempData) <- threshName
    plotDataMean <- cbind(plotDataMean, tempData)
    names(plotDataMean) <- c(names(plotData), threshName)
  }
  
  varNames <- sapply(strsplit(xyAtts, "_"), `[[`, 1)
  varUnits <- getVarUnits(varNames)
  xyLabels <- paste0(xyAttDefs, " (", varUnits, ")")
  xlimits <- c(min(plotDataMean[ ,1]), max(plotDataMean[ ,1]))
  ylimits <- c(min(plotDataMean[ ,2]), max(plotDataMean[ ,2]))
  
  p1 <- ggplot() +
    geom_tile(data = plotDataMean, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], fill = .data[[perfName]])) +
    labs(x = xyLabels[1], y = xyLabels[2]) +
    scale_x_continuous(expand=c(0, 0)) +    # no extra space on x and y axes
    scale_y_continuous(expand=c(0, 0)) +
    coord_cartesian(xlim=xlimits, ylim=ylimits) + 
    theme_heatPlot()
  
  # contours
  if(perfSpace_contours) {  # TRUE or FALSE default setting
    
    contourBreaks <- pretty(colLimIn, perfSpace_nContour)
    # not able to add bins here - fix later
    p1 <- p1 + geom_contour(data = plotDataMean, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], z = .data[[perfName]]), colour = "black", breaks = contourBreaks)  + #,breaks=seq(0.6,0.8,0.01)
    directlabels::geom_dl(data = plotDataMean, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], z = .data[[perfName]], label = stat(level)), #edited 20/06/2018
                           method = list("first.points", "calc.boxes", "enlarge.box", box.color = NA, fill = "transparent", vjust=-0.5,hjust=-0.5, "draw.rects"),stat="contour", breaks = contourBreaks)  #,breaks=seq(0.6,0.8,0.01)
  }
  
  if (!is.null(perfThresh)) {
    p1 <- addContourThreshold(p1, plotDataMean, perfThresh, perfThreshLabel, xyAtts, perfName)

  }
                                                                  # original lim  # modified based on climData 
  p2List <- addClimData(p1, climData, perfName, xyAtts, xlimits, ylimits, colLim, colLimIn)
  p2 <- p2List[[1]]
  colLimIn <- p2List[[2]]
  
  # Does not work if there are more than one contour line
  # if (!is.null(perfThresh)) {
  #   p2 <- p2 +
  #   directlabels::geom_dl(data = plotDataMean, aes(x = !!as.symbol(xyAtts[1]), y = !!as.symbol(xyAtts[2]), z = !!as.symbol(perfName), label = perfThreshLabel), #edited 20/06/2018
  #                         #method = list("first.points", "calc.boxes", "enlarge.box", box.color = "black", fill = "white", "draw.rects"),stat="contour", breaks = perfThresh)  #,breaks=seq(0.6,0.8,0.01)
  #                         # method = list(box.color = NA, "angled.boxes"),stat="contour", breaks = perfThresh)  #,breaks=seq(0.6,0.8,0.01)
  #                         method = list("smart.grid", box.color = "black", fill = "white", "draw.rects"), stat = "contour", breaks = perfThresh)
  #   #p2 <- directlabels::direct.label(p2, method = "angled.boxes", stat = "contour", breaks = perfThresh)
  # }
  
  if (is.null(colMap)) {
    coloursIn <- foreSIGHT.colmap(perfSpace_nlevel)
  } else {
    coloursIn <- colMap
  }
  
  # p2 <- p2 + scale_fill_gradientn(colours = grDevices::adjustcolor(coloursIn, alpha.f = perfSpace_alpha), limits = colLimIn,
  #                                 guide = guide_colorbar(title = perfName, title.position = "right", order = 1, barwidth = 12, barheight = 0.6)) + labs(tag = tag_text)
  
  p2 <- p2 + scale_fill_gradientn(colours = coloursIn, limits = colLimIn,
                                  guide = guide_colorbar(title = perfName, title.position = "right", order = 1, barwidth = 12, barheight = 0.6)) + labs(tag = tag_text)
  
  #print(p2) 
  return(p2)
  
}


addContourThreshold <- function(p1, plotDataMean, perfThresh, perfThreshLabel, xyAtts, perfName){
  
  
  p1 <- p1 + geom_contour(data = plotDataMean, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], z = .data[[perfName]]), 
                          breaks = perfThresh, colour = perfSpace_threshCol, size = threshLineSize)  #+ #,breaks=seq(0.6,0.8,0.01)
  # directlabels::geom_dl(data = plotDataMean, aes(x = !!as.symbol(xyAtts[1]), y = !!as.symbol(xyAtts[2]), z = !!as.symbol(perfName), label = perfThreshLabel), #edited 20/06/2018
  #                       #method = list("first.points", "calc.boxes", "enlarge.box", box.color = "black", fill = "white", "draw.rects"),stat="contour", breaks = perfThresh)  #,breaks=seq(0.6,0.8,0.01)
  #                       # method = list(box.color = NA, "angled.boxes"),stat="contour", breaks = perfThresh)  #,breaks=seq(0.6,0.8,0.01)
  #                       method = list("smart.grid", box.color = "black", fill = "white", "draw.rects"), stat = "contour", breaks = perfThresh)

  # Had to take this roundabout way to label the threshold - rethink
  # performance threshold
  nx <- length(unique(plotDataMean[, 1]))
  ny <- length(unique(plotDataMean[, 2]))
  look <- fields::as.image(plotDataMean[ ,3], ind = cbind(plotDataMean[ ,1], plotDataMean[ ,2]), nx = nx, ny = ny)
  lineThresh <- contourLines(x = look$x, y = look$y, z = look$z, levels = perfThresh)
  if (!(length(lineThresh) == 0)) {
    # lineNo <- NULL
    # lineLen <- NULL
    # if (length(lineThresh) > 1) {
    #   for (i in 1:length(lineThresh)) {
    #     nPts <- length(lineThresh[[i]])
    #     if (!(lineThresh[[i]]$x[1] == lineThresh[[i]]$x[nPts] | lineThresh[[i]]$y[1] == lineThresh[[i]]$y[nPts])) {
    #       lineLenNew <- length(lineThresh[[i]]$x)
    #       if (!is.null(lineNo)) {
    #         if (lineLenNew > lineLen) {
    #           lineNo <- i
    #           lineLen <- lineLenNew
    #         }
    #       } else {
    #         lineNo <- i
    #         lineLen <- lineLenNew
    #       }
    #     }
    #   }
    # } else {
    #   lineNo <- 1
    # }
    lineNo <- 1
    
    perfThreshData <- data.frame(x = lineThresh[[lineNo]]$x, y = lineThresh[[lineNo]]$y)
    perfThreshData$Tname <- perfThreshLabel
    nLines <- nrow(perfThreshData)
    
    p1 <- p1 + #geom_line(perfThreshData, mapping = aes(x = x, y = y, colour = Tname), colour = perfSpace_threshCol, size = threshLineSize) +
      geom_label(data = perfThreshData[nLines, ], mapping = aes(x = .data$x, y = .data$y, label = .data$Tname), size = threshLabelSize, vjust = "inward", hjust = "inward") #, colour = threshCol, fill = "black")
  } else {
    message("perfThresh is not plotted becasue it is outside the performance space.")
  }
  return(p1)
  
}

# Add contour function for plotOptions. The contour based on the original performance metric values is added on top of heatPlots created using the difference data 
# So the threshold is a line based on the original performance metric - added on heatMaps plotted using the differences (option 2 - option 1)
# Similar to this function, ggplot_build can be used instead of fields::look to get the contour lines from a single contour line
addThresholdLines <- function(p1, plotDataMean, perfThresh, perfThreshLabel, xyAtts, perfName, lineCol, lineSize, lineAlpha, label = "beginning"){
  
  p1 <- p1 + geom_contour(data = plotDataMean, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], z = .data[[perfName]]), 
                          breaks = perfThresh, colour = lineCol, size = lineSize, alpha = lineAlpha)  #+ #,breaks=seq(0.6,0.8,0.01)
  p2 <- ggplot()+geom_contour(data = plotDataMean, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], z = .data[[perfName]]), 
                              breaks = perfThresh, colour = lineCol, size = lineSize, alpha = lineAlpha)  #+ #,breaks=seq(0.6,0.8,0.01)
  # Code to get the x-y co-ordinates of the contour line from p2 and plot it on p1
  p2Build <- ggplot_build(p2)
  threshContour <- data.frame(x = p2Build[["data"]][[1]]$x, y = p2Build[["data"]][[1]]$y)
  #threshContour <- as.data.frame(spline(p2Build[["data"]][[1]]$x, p2Build[["data"]][[1]]$y))
  
  if (nrow(threshContour) > 0) {
    if (label == "beginning") {
      nLine <- 1
    } else {
      nLine <- nrow(threshContour)-1
      if (nLine == 0) nLine = 1
    }
    
    threshContour$Tname <- perfThreshLabel
    p1 <- p1 + #geom_line(threshContour, mapping = aes(x = x, y = y), colour = lineCol, size = threshLineSize) +
      geom_label(data = threshContour[nLine, ], mapping = aes(x = .data$x, y = .data$y, label = .data$Tname), size = threshLabelSize, vjust = "inward", hjust = "inward") #, colour = threshCol, fill = "black")
  
  # # threshold label
  # nx <- length(unique(plotDataMean[, 1]))
  # ny <- length(unique(plotDataMean[, 2]))
  # look <- fields::as.image(plotDataMean[ ,3], ind = cbind(plotDataMean[ ,1], plotDataMean[ ,2]), nx = nx, ny = ny)
  # lineThresh <- contourLines(x = look$x, y = look$y, z = look$z, levels = perfThresh)
  # if (!(length(lineThresh) == 0)) {
  #   lineNo <- 1
  #   
  #   perfThreshData <- data.frame(x = lineThresh[[lineNo]]$x, y = lineThresh[[lineNo]]$y)
  #   perfThreshData$Tname <- perfThreshLabel
  #   nLines <- nrow(perfThreshData)
  #   
  #   p1 <- p1 + #geom_line(perfThreshData, mapping = aes(x = x, y = y, colour = Tname), colour = perfSpace_threshCol, size = threshLineSize) +
  #     geom_label(aes(x = x, y = y, label = Tname), data = perfThreshData[nLines, ], size = threshLabelSize, vjust = "inward", hjust = "inward") #, colour = threshCol, fill = "black")
  } else {
    message("perfThresh is not plotted becasue it is outside the performance space.")
  }
  return(p1)
  
}



