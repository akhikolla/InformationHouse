
#' Plots performance for one-at-a-time (OAT) perturbations in attributes
#'
#' \code{plotPerformanceOAT} uses the system model performance calculated using the function \code{runSystemModel} and 
#' the summary of the simulation generated using the function \code{generateScenarios} & \code{getSimSummary} as input.
#' The function creates line plots, each panel shows the variations in performance with perturbations in a single attribute. The function is intended for 
#' use with simulations with attributes perturbed on a one-at-a-time (OAT) grid.
#' @param performance a named list; contains the system model performance calculated using \code{runSystemModel}. 
#' If the list contains more than one performance metric, the first metric will be plotted.
#' @param metric a string; the name of the performance metric to be plotted. The argument can be used to select a metric from performance for plotting. 
#' @param sim a list; a summary of a simulation containing the scenarios generated using the function \code{generateScenarios} that is used to run the system model using \code{runSystemModel}.
#' The summary may be obtained using the function \code{getSimSummary}
#' @param topReps an integer (default = NULL); the number of "top" replicates to be used. The "top" replicates will be identified for each target based on the simulation fitness. 
#' The average performance across \code{topReps} replicates will be plotted.
#' @param col a colour; the colour of the lines. If \code{NULL}, the a default colour is used.
#' @param ylim a vector of 2 values; the minimum and maximum limits of the y-axis (performance) scale.
#' @details The plots show the mean value of performance across replicates. The ranges between the minimum and maximum values of performance across replicates are shaded.
#' The function is intended for use with simulations containing attributes perturbed on an "OAT" grid. If the perturbations are on a "regGrid", this function will subset 
#' OAT perturbations, if available, to create the plots. The function creates separate plots for perturbations in attributes of temperature and other variables.
#' The function may be called with \code{performance} argument specifying the metric to be plotted to plot other metrics.
#' @return The plot of the performance space and the ggplot object.
#' @seealso \code{runSystemModel}, \code{generateScenarios}, \code{plotPerformanceSpace}, \code{getSimSummary}
#' @examples
#' # load example datasets
#' data("egSimSummary")
#' data("egSimPerformance")
#' plotPerformanceOAT(egSimPerformance[2], egSimSummary)
#' plotPerformanceOAT(egSimPerformance[1], egSimSummary)
#' # using the metric argument
#' plotPerformanceOAT(egSimPerformance, egSimSummary, metric = "Reliability (-)")
#' @export

# Need to think about whether the function would work for "regGrid" perturbations - in this case perturbations would be averaged across perturbations in other attributes
# Impose restriction for "OAT" later if required

plotPerformanceOAT <- function(performance,                   # system model performance, matrix of size targets x replicates
                               sim,                           # simulation containing all targets and replicates
                               metric = NULL,                 # the name of the performance metric to be plotted
                               topReps = NULL,                   # number of topReps based on fitness (i.e., the objective function score which is -ve)
                               # perfThresh = NULL,             # desired performance threshold; plot would contain a contour to mark this threshold
                               # perfThreshLabel = "Threshold", # label text for the threshold
                               # do we need this - attSlices = NULL,              # list containing the slices of attributes to use for plotting
                               # climData = NULL,               # changes in climate attributes from other sources - can include label. If the performance measure being plotted is a column in the data.frame, the points will be coloured accordingly
                               col = NULL,                    # colour of the ribbon
                               ylim = NULL                    # ylim of the data, xlim is determined by the perturbation range
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
  
  if (is.list(sim[["controlFile"]])) {
    # get the simulation fitness if stochastic
    simFitness <- getSimFitness(sim)
    
    # check that the performance is calculated from this simulation
    if (!identical(dim(simFitness), dim(performance[[1]]))) {
      stop("The number of targets and replicates in sim and performance should match. Is the performance calculated using sim?")
    }
  } else {
    # scaling
    simFitness <- NULL
  }
  
  # unpacking sim metadata
  repNames <- names(sim[grep("Rep", names(sim))])
  #if (is.null(repNames)) stop("There are no replicates in sim.")
  #tarNames <- names(sim[[repNames[1]]])
  nRep <- length(repNames)
  #nTar <- length(tarNames)
  targetMat <- sim[["expSpace"]][["targetMat"]]
  
  if (!is.null(topReps)) {
    if (topReps > nRep) {
      message(paste0("sim does not contain ", topReps, " replicates. Using all the replicates available in sim.\n"))
    }
  }
  
  # Use code later if adding an attSlices argument
  #======================================================
  # # subset targets if attSlices are specified
  # if (!is.null(attSlices)) {
  #   # indices to subset the rows of the targetMat 
  #   tarInd <- getSliceIndices(sim[["expSpace"]], attSlices)
  #   
  #   # subset targetMat
  #   targetMat <- targetMat[tarInd, ]
  # } else {
  #   tarInd <- NULL
  # }
  
  
  # get attPerturb with atleast 2 samples
  # doing this here instead of using attPerturbSamp directly since the targetMat may be subsetted
  attPerturb <- getAttPerturb(targetMat)
  if (is.null(attPerturb)) stop("The simulation does not contain OAT perturbed attributes to plot.")
  
  # identify x and y columns
  # attNames <- colnames(targetMat)
  # colAttX <- which(attNames == attX)
  # colAttY <- which(attNames == attY)
  
  nPerf <- length(performance)
  
  if (nPerf > 1) message("performance contains more than one performance metric, the first metric is plotted.")
  
  # loop over multiple metrics - can do this, but becomes complicated due for multiple colMap, colLim, and threshold labels
  #for (p in 1) {
    # subset performance based on slices
    # if (!is.null(tarInd)) {
    #   perfMatrix <- performance[[p]][tarInd, ]
    # } else {
    #   perfMatrix <- performance[[p]]
    # }
    
    perfMatrix <- performance[[1]]
    
    # get the average performance to be plotted
    performanceAv <- getPerfStat(perfMatrix, simFitness, topReps, nRep, statFUN = mean)
    pMin <- getPerfStat(perfMatrix, simFitness, topReps, nRep, statFUN = min)
    pMax <- getPerfStat(perfMatrix, simFitness, topReps, nRep, statFUN = max)
    
    # name appropriately
    names(performanceAv) <- perfName[1]
    # "pMin" and "pMax" are used in the OATPlot code
    # If changing, needs to be changed there
    names(pMin) <- "pMin"
    names(pMax) <- "pMax"
    
    # too few replicates to calculate min-max
    if (nRep < 3) {
      pMin <- NULL
      pMax <- NULL
    }
    
    # create data.frame for plotting
    plotData <- getOATData(attPerturb, targetMat, performanceAv, pMin, pMax)
    
    perfPlots <- lapply(plotData, OATPlot, col = col, ylimits = ylim)
    # for (i in 1:length(plotData)) {
    #   perfPlots[[i]] <- OATPlot(plotData[[i]], col = col, ylimits = ylim)
    #                             # Add later if required
    #                             # perfThresh = perfThresh, perfThreshLabel = perfThreshLabel, climData = climData)
    #   print(perfPlots[[i]])
    # }
  #}
  return(perfPlots)
}


getAttPerturb <- function(targetMat) {   # data.frame as input
  attNames <- colnames(targetMat)
  nAtt <- ncol(targetMat)
  attPerturb <- NULL
  
  for (i in 1:nAtt) {
    attVal <- targetMat[ ,i]
    nPert <- length(unique(attVal))
    if (nPert > 1) attPerturb <- c(attPerturb, attNames[i])
  }
  return(attPerturb)
}


getOATData <- function(attPerturb,   # vector; perturbed attNames
                      targetMat,    # data.frame; targets x attributes
                      pAv,          # vector; average performance at targets
                      pMin = NULL,  # vector; minimum performance at t      
                      pMax = NULL   # vector; maximum performance at t
) {
  attNames <- colnames(targetMat)
  nAtt <- ncol(targetMat)
  
  # no. of dfs = no.of atts that have OAT perturbations
  count <- 0
  OATDataList <- list()
  OATAttVar <- NULL
  for (i in 1:length(attPerturb)) {
    icol <- which(attNames == attPerturb[i])
    attVar <- strsplit(attPerturb[i], "_")[[1]][1]
    attOther <- attNames[-icol]
    
    iIndList <- list()
    # Index where other attributes are not perturbed
    for (j in 1:length(attOther)) {
      othVar <- strsplit(attOther[j], "_")[[1]][1]
      if (othVar == "Temp") {
        noPert = 0
      } else {
        noPert = 1
      }
      jcol <- which(attNames == attOther[j])
      targetJ <- targetMat[ ,jcol]
      jInd <- which(targetJ == noPert)
      iIndList[[j]] <- jInd
    }
    
    # Intersect all index where the other attributes are not perturbed
    iInd <-  Reduce(intersect, iIndList)
    
    if(!(identical(iInd, numeric(0)))) {
      if(!identical(iInd, integer(0))) {
        count <- count + 1
        if (!is.null(pMin) & !is.null(pMax)) {
          OATDataList[[count]] <- data.frame(targetMat[ ,icol][iInd], pAv[iInd, 1], pMin[iInd, 1], pMax[iInd, 1])
          names(OATDataList[[count]]) <- c("perturbation", paste0("Performance metric: ", names(pAv)), names(pMin), names(pMax))
        } else {
          OATDataList[[count]] <- data.frame(targetMat[ ,icol][iInd], pAv[iInd, 1])
          names(OATDataList[[count]]) <- c("perturbation", paste0("Performance metric: ", names(pAv)))
        }
        OATDataList[[count]][["attribute"]] <- attPerturb[i]
        OATAttVar <- c(OATAttVar, attVar)
      }
    }
  }
  
  if (count == 0) stop("The simulation does not contain OAT perturbed attributes to plot.")
  
  #***********
  # CODE TO DIVIDE DATA INTO APPROPRIATE DATA.FRAMES FOR PLOTTING
  # CAREFUL IF MAKING CHANGES
  
  # Default: all vars except "Temp" will be stored in the data.frame no. 1, "Temp" in data.frame no. 2
  othrInd <- 1
  tempInd <- 2
  
  if (sum(OATAttVar %in% c("Temp")) > 0) {
    # no. of data frames = 2 if atleast two types of var exist and atleast one of them is "Temp"
    if (length(unique(OATAttVar)) > 1) {
      nDf <- 2
    } else {
      # only "Temp" exists
      nDf <- 1
      tempInd <- 1
    }
  } else {
    # only othrVar exists
    nDf <- 1
  }
  
  
  # check if two figures are required for otherVar, assume that "Temp" would anyway fit in one fig
  # if otherVars exists
  if (othrInd != tempInd) {
    # no. of othrAtts
    nAtts <- length(OATAttVar) - sum(OATAttVar %in% c("Temp"))
    if (nAtts > 8) {
      attDiv <- ceiling(nAtts/2)
      attGrp1 <- 1:attDiv
      
      # where should the new Grp be saved
      # if already two exists
      if (nDf == 2) {
        othrInd2 <- 2
        tempInd <- 3
        nDf <- 3
      } else {
        othrInd2 <- 2
        nDf <- 2
      }
    } else {
      attGrp1 <- 1:nAtts
    }
  }
  
  OATdf <- vector(mode = "list", length = nDf)
  
  # combine data.frames for all the OATAtts
  # column names of the combined df are: perturbation, names of the corresponding performance stats, and attribute
  # initialise
  for ( i in 1:nDf) {
    OATdf[[i]] <- OATDataList[[1]]
    OATdf[[i]] <- OATdf[[i]][FALSE, ]
  }
  countOthr <- 0
  for (i in 1:length(OATDataList)) {
    if (OATAttVar[i] %in% c("Temp")) {
      OATdf[[tempInd]] <- rbind(OATdf[[tempInd]], OATDataList[[i]])
    } else {
      countOthr <- countOthr + 1
      if (countOthr %in% attGrp1) {
        OATdf[[othrInd]] <- rbind(OATdf[[othrInd]], OATDataList[[i]])
      } else {
        OATdf[[othrInd2]] <- rbind(OATdf[[othrInd2]], OATDataList[[i]])
      }
    }
  }
  
  #***************
  return(OATdf)
}

OATPlot <- function(plotData, col = NULL, ylimits = NULL) {
  
  dfNames <- colnames(plotData)
  # the second column is always performance (named according to the input data)
  perfName <- dfNames[2]
  
  varName <- strsplit(plotData[["attribute"]][1], "_")[[1]][1]
  varUnits <- getVarUnits(varName)
  
  attNames <- unique(plotData[["attribute"]])
  attFullNames <- mapply(tagBlender_noUnits, attNames)
  
  # aggregate data if required
  plotDataMean <- aggregate(.~perturbation+attribute, plotData, mean)
  xLabeltext <- paste0("Perturbation", " (", varUnits, ")")
  
  if (is.null(col)) {
    col = foreSIGHT.colmap(1)
  }
  
  p1 <- ggplot(data = plotDataMean, aes(x = .data$perturbation, y = .data[[perfName]]))
  
  # UNSMOOTHENED VERSION 
  #------------------------------------------------------------
  # add ribbon if data exists
  if (!is.null(plotDataMean[["pMin"]])) {
    p1 <- p1 + geom_ribbon(aes(ymin = .data$pMin, ymax = .data$pMax), fill = col, alpha = OATplot_fillAlpha)
  }

  p1 <- p1 + geom_line(colour = col, size = OATplot_lineSize) + #, method = "loess", se = FALSE) +
    facet_wrap(vars(.data$attribute), scales = "free_x", nrow = 1, strip.position = "top", labeller = as_labeller(attFullNames, default = label_wrap_gen(35))) +
    theme_heatPlot(OATplot_textSize) + theme(strip.placement = "inside", strip.background = element_rect(fill = NA)) + xlab(xLabeltext) +
    theme(strip.text.x = element_text(size = OATplot_textSize, angle = 0)) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(tag = tag_text)
          #axis.text.x = element_text(color = "black", size = textSize, face = "plain", vjust = 0))
  #--------------------------------------------------------------
  
  # Fix smoothing of geom ribbon later
  #==============================================================
  # pMean <- p1 + stat_smooth(colour = col, size = OATplot_lineSize, method = "loess", se = FALSE) + facet_wrap(vars(attribute))
  # ggPMean <- ggplot_build(pMean)
  # 
  # pMin <- p1 + stat_smooth(aes(x = perturbation, y = pMin), colour = col, size = OATplot_lineSize, method = "loess", se = FALSE) + facet_wrap(vars(attribute))
  # ggPMin <- ggplot_build(pMin)
  # 
  # pMax <- p1 + stat_smooth(aes(x = perturbation, y = pMax), colour = col, size = OATplot_lineSize, method = "loess", se = FALSE) + facet_wrap(vars(attribute))
  # ggPMax <- ggplot_build(pMax)
  # 
  # p1 <- p1 + stat_smooth(colour = col, size = OATplot_lineSize, method = "loess", se = FALSE) +
  #   facet_wrap(vars(attribute), scales = "free_x", nrow = 1, strip.position = "top", labeller = as_labeller(attFullNames, default = label_wrap_gen(20))) +
  #   theme_heatPlot(OATplot_textSize) + theme(strip.placement = "inside", strip.background = element_rect(fill = NA)) + xlab(xLabeltext) +
  #   theme(strip.text.x = element_text(size = OATplot_textSize, angle = 0)) +
  #   theme(panel.spacing = unit(1, "lines")) +
  #   labs(tag = tag_text)
  #         #axis.text.x = element_text(color = "black", size = textSize, face = "plain", vjust = 0))
  # 
  # p1Build <- ggplot_build(p1)
  # 
  # # add ribbon if data exists
  # if (!is.null(plotDataMean[["pMin"]])) {
  #   
  #   # pMin <- ggplot(data = plotDataMean, aes(x = perturbation, y = !!as.symbol(perfName))) + 
  #   # p1 <- p1 + geom_ribbon(aes(ymin = pMin, ymax = pMax), fill = col, alpha = OATplot_fillAlpha)
  # }
  #=================================================================
  
  if (!is.null(ylimits)) {
    p1 <- p1 + ylim(ylimits)
  }
  print(p1)
  return(p1)
  # return(list(p1, p1Build, ggPMean, ggPMin, ggPMax))
}





