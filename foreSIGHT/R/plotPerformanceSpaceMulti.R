
# example of performance thersholds
# perfThreshMin <- c(NA, 0.81)
# perfThreshMax <- c(29, NA)
# add attSlices to the args of this function as well

# is perturbation space the correct word to use?

#' Plots contours of the number of performance thresholds exceeded in the perturbation space
#'
#' \code{plotPerformanceSpaceMulti} uses multiple system model performances calculated using the function \code{runSystemModel} and 
#' the summary of the simulation generated using the functions \code{generateScenarios} & \code{getSimSummary} as input to plot filled contours showing the 
#' number of performance thresholds exceeded in the perturbation space.
#' The user may specify the attributes to be used as the axes of the perturbation space.
#' @inheritParams plotPerformanceSpace
#' @param performance a list; each element of the list should be a performance metric. May be calculated using the function \code{runSystemModel}
#' @param perfThreshMin a vector; the minimum threshold value of each performance metric. The length of the vector should be equal to \code{length(performance)}. 
#' If the metric does not have a minimum threshold, specify the corresponding element in \code{perfThreshMin} as \code{NA}.
#' @param perfThreshMax a vector; the maximum threshold value of each performance metric. The length of the vector should be equal to \code{length(performance)}.  
#' If the metric does not have a maximum threshold, specify the corresponding element in \code{perfThreshMax} as \code{NA}.
#' @param climData data.frame; the values of \code{attX} and \code{attY} from other sources like climate models. This data will be plotted as points in the perturbation space.
#' If the \code{Name} of the data is available in the data.frame, the points will be identified using the \code{Name}. 
#' Please refer data provided with the package that may be loaded using \code{data("egClimData")} for an example of the expected format of \code{climData}.
#' @param col a vector of colours; The length of the vector should atleast be sufficient to assign unique colours to all
#' the different values in the generated plot. If \code{NULL}, the default foreSIGHT colours is used.
#' @details If the space contains more than two perturbed attributes, the performance values are averaged across the perturbations in the attributes other than \code{attX} and \code{attY}.
#' The user may specify argument \code{attSlices} to slice the performance space at specific values of the other perturbed attributes. If \code{attSlices} are used to 
#' specify minimum-maximum values to subset other perturbed attributes, the performance values are averaged across the subsetted perturbations in these attributes. This function
#' cannot be used with \code{sim} perturbed on an "OAT" grid since contours of the number of performance thresholds exceeded cannot be calculated for an irregular perturbation space.
#' @return The plot showing the number of thresholds exceeded and the ggplot object.
#' @seealso \code{runSystemModel}, \code{generateScenarios}, \code{getSimSummary}, \code{plotPerformanceSpace}
#' @examples
#' # load example datasets
#' data("egSimPerformance")
#' data("egSimSummary")
#' data("egClimData")
#' 
#' plotPerformanceSpaceMulti(egSimPerformance, egSimSummary, 
#' perfThreshMin = c(NA, 0.80), perfThreshMax = c(30, NA))
#' 
#' # add alternate climate data and specify different colours for the plot
#' plotPerformanceSpaceMulti(egSimPerformance, egSimSummary, perfThreshMin = c(NA, 0.80), 
#' perfThreshMax = c(30, NA), climData = egClimData, col = viridisLite::magma(3))
#' 
#' # example using simple scaled simulations
#' data("egScalPerformance")
#' data("egScalSummary")
#' data("egClimData")
#' plotPerformanceSpaceMulti(egScalPerformance, egScalSummary, perfThreshMin = c(NA, 0.80), 
#' perfThreshMax = c(30, NA), climData = egClimData)
#' @export

plotPerformanceSpaceMulti <- function(performance,      # system model performance, list containing matrices of different performance metrics size targets x replicates 
                          sim,                      # summary of simulation containing metadata of all targets and replicates
                          perfThreshMin,            # vector; minimum performance thresholds (NA if thresholds don't apply)
                          perfThreshMax,            # vector; maximum performance thresholds (NA if thresholds don't apply)
                          attX = NULL,              # attribute to be plotted on x-axis
                          attY = NULL,              # attribute to be plotted on y-axis
                          attSlices = NULL,         # list containing the slices of attributes to use for plotting
                          topReps = NULL,           # number of topReps based on fitness (i.e., the objective function score which is -ve)
                          climData = NULL,          # data.frame containing alternate climate data
                          col = NULL                # vector of colours
) {
  
  # checks specific to this function
  if (!is.list(performance)) stop("performance should be a list containing all the performance metrics to be analysed")
  nMet <- length(performance)
  if (nMet != length(perfThreshMin)) {
    stop("perfThreshMin should be specified for all performance metrics. Use NA if minimum threshold does not exist for the metric.")
  }
  if (nMet != length(perfThreshMax)) {
    stop("perfThreshMax should be specified for all performance metrics. Use NA if maximum threshold does not exist for the metric.")
  }
  if (sum(perfThreshMax <= perfThreshMin, na.rm = TRUE) > 0) {
    stop("perfThreshMax should be greater than perfThreshMin.")
  }
  
  # unpacking sim metadata
  repNames <- names(sim[grep("Rep", names(sim))])
  tarNames <- names(sim[[repNames[1]]])
  nRep <- length(repNames)
  nTar <- length(tarNames)
  targetMat <- sim[["expSpace"]][["targetMat"]]
  
  if (sim[["expSpace"]][["attPerturbType"]] == "OAT") stop("Performance thresholds cannot be plotted for simulations using an OAT sampling.")
  
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
  
  if (is.list(sim[["controlFile"]])) {
    if (!is.null(topReps)) {
      if (topReps > nRep) {
        cat(paste0("sim does not contain ", topReps, " replicates. Using all replicates available in sim.\n"))
        topReps <- nRep
      }
    }
    # get the simulation fitness
    simFitness <- getSimFitness(sim)
  } else {
    if (sim[["controlFile"]] == "scaling") {
      simFitness <- NULL
      topReps <- NULL
    } else {
      stop("sim$controlFile is unrecognized.")
    }
  }
  
  attNames <- colnames(targetMat)
  colAttX <- which(attNames == attX)
  colAttY <- which(attNames == attY)
  tarAttX <- targetMat[ ,colAttX]
  tarAttY <- targetMat[ ,colAttY]
  
  # Find out how many unique targets there are (i.e., in case there are repeated targets)
  tempMat <- cbind(targetMat[colAttX], targetMat[colAttY])
  nTarUpdated <- nrow(tempMat[!duplicated(tempMat), ])
  
  # identify exceeded thresholds
  perfThresExceed <- matrix(0, nrow = nTarUpdated, ncol = nMet) 
  for (m in 1:nMet) {
    # check that the performance is calculated from this simulation
    if (!is.null(simFitness)) {
      if (!identical(dim(simFitness), dim(performance[[m]]))) {
        stop(paste0("The number of targets and replicates in sim and performance[[",  nMet, "]] should match. Is the performance metric calculated using sim?"))
      }
    }
    
    if (!is.null(tarInd)) {
      performance[[m]] <- performance[[m]][tarInd, ]
    }
    
    # get the average performance
    performanceAv <- getPerfStat(performance[[m]], simFitness, topReps, nRep, statFUN = mean)
    
    # Note: for plotPerformanceSpace, this calculation is done inside the heatPlot
    #-------
    # sometimes multiple targets exist with the same attX and attY combinations (eg: another perturbed attributes, or OAT sampling)
    # averaging the performance across these multiple targets before checking threshold exceedence
    tempMat <- cbind(targetMat[colAttX], targetMat[colAttY], performanceAv) 
    names(tempMat) <- c("x", "y", "z")
    performanceAvAgg <- aggregate(.~x+y, tempMat, mean)
    names(performanceAvAgg) <- c(attX, attY, "Performance")
    #-------
    
    # check performance against thresholds
    pMin <- perfThreshMin[m]
    pMax <- perfThreshMax[m]
    if (!is.na(pMin)) {
      ind <- which(performanceAvAgg[,3] < pMin)
      perfThresExceed[ind, m] <- 1
    }
    if (!is.na(pMax)) {
      ind <- which(performanceAvAgg[,3] > pMax)
      perfThresExceed[ind, m] <- 1
    }
  }
  # sum up exceeded thresholds for all metrics
  perfThreshComb <- rowSums(perfThresExceed)
  
  nx <- length(unique(tarAttX))
  ny <- length(unique(tarAttY))
  
  #                             attX,    attY, 
  threshPlotData <- cbind(performanceAvAgg[ ,1:2], perfThreshComb)
  
  p <- fillHeatPlot(threshPlotData, colMap = col, climData = climData)
  # plotThreshContour(threshPlotData, nx, ny)
  # return(threshPlotData)
  return(p)
  
}

plotThreshContour <- function(plotData, nx, ny){
  
  xyAtts <- colnames(plotData)[1:2]
  xyAttDefs <- mapply(tagBlender_noUnits, xyAtts, USE.NAMES = FALSE)
  
  varNames <- sapply(strsplit(xyAtts, "_"), `[[`, 1)
  varUnits <- getVarUnits(varNames)
  xyLabels <- paste0(xyAttDefs, " (", varUnits, ")")
  
  # max no. of thresholds
  nTh <-  max(plotData[,3])
  
  par(mar=c(5.5,3.8,1,1),oma=c(1,0.3,0.3,0.3), mgp = c(0,0.5,0))
  #plot(NA,NA,xlim = c(min(plotData[ ,1]), max(plotData[ ,1])), ylim = c(min(plotData[ ,2]), max(plotData[ ,2])),xaxs="i",yaxs='i',ylab="",xlab="")
  fields::quilt.plot(x = plotData[ ,1], y = plotData[ ,2], z = plotData[ ,3], nx = nx ,ny = ny,
                     add.legend = FALSE, nlevel = perfSpace_nlevel, col = foreSIGHT.colmap(perfSpace_nlevel), 
                     xlim = c(min(plotData[ ,1]), max(plotData[ ,1])), ylim = c(min(plotData[ ,2]), max(plotData[ ,2])),
                     zlim = c(0, nTh), xlab = "", ylab ="")
  box(col = "black")
  
  mtext(side=1,text=xyLabels[1],line=1.8)
  mtext(side=2,text=xyLabels[2],line=1.8)
  
  
  #get image for contours
  look <- fields::as.image(plotData[ ,3], ind = cbind(plotData[ ,1], plotData[ ,2]), nx = nx, ny = ny)
  #filled.contour(x = look$x, y = look$y, z = look$z, zlim = c(0, nTh), levels = c(0:nTh+1), col = foreSIGHT.colmap(nTh))
  
  #               method="edge", labcex = 1, nlevels = perfSpace_ncontour)
  contour(add = TRUE, x = look$x, y = look$y, z = look$z, method="edge", labcex = 1, levels = c(0:nTh+1), lwd=2)
  #contour(add = TRUE, x = look$x, y = look$y, z = look$z, nlevels = perfSpace_ncontour)
  # points(x=climAtts$P_ann_seasRatio_m[1:climLim],y=climAtts$P_ann_tot_m[1:climLim],pch=".",cex=2.0)
  
  #if 3 or 4 add ramps
  # if(colBar){
  #add legend
  # image.plot(legend.only = TRUE, zlim = colLim, col = foreSIGHT.colmap(perfSpLevel),
  #            horizontal = TRUE, legend.lab = metric.lab,smallplot=c(0.2,0.9,0.0001,0.015),legend.mar = 2)
  # fields::image.plot(legend.only = TRUE, zlim = c(0, nTh), col = foreSIGHT.colmap(perfSpace_nlevel),
  #                    horizontal = TRUE, smallplot=c(0.2,0.9,0.0001,0.02))
  #panel.lab <- colnames(plotData)[3]
  #add label
  #mtext(text = panel.lab,side = 1,line = 2.3,at =1.1,cex = 1.25,las=2)
  # }
  # legend at base
  #par(xpd=TRUE)
  legend("bottom",xpd = TRUE, inset = c(0, -0.25), legend=seq(0,nTh),bty="n",horiz = TRUE,
         fill = foreSIGHT.colmap(nTh+1), title = "No. of thresholds exceeded")
  #        #fill=c(makeTransparent(asc.col[1], alpha=40),makeTransparent(asc.col[nlevel], alpha=100),makeTransparent(asc.col[nlevel], alpha=200)))
  #mtext(side=1,text="No. of thresholds exceeded",adj=1.0, line=4,at=1.04)
  
  mtext(tag_text, side=1, line=0, adj=1.0, cex=0.8, col=tag_textCol, outer=TRUE)
  
}


# plotting code for a filled contour plot similar to heatPlot
fillHeatPlot <- function(plotData, colMap = NULL, climData = NULL) {
  
  xyAtts <- colnames(plotData)[1:2]
  xyAttDefs <- mapply(tagBlender_noUnits, xyAtts, USE.NAMES = FALSE)
  
  varNames <- sapply(strsplit(xyAtts, "_"), `[[`, 1)
  varUnits <- getVarUnits(varNames)
  xyLabels <- paste0(xyAttDefs, " (", varUnits, ")")
  xlimits <- c(min(plotData[ ,1]), max(plotData[ ,1]))
  ylimits <- c(min(plotData[ ,2]), max(plotData[ ,2]))
  
  nThresh <- sort(unique(plotData[,3]))
  # Add one more
  maxT <- max(nThresh) + 1
  breaks <- c(nThresh, maxT)  - 0.00001
  
  threshPlot_alpha <- 0.7
  
  p1 <- ggplot(data = plotData, aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]])) +
    stat_contour_filled(aes(z = .data$perfThreshComb), breaks = breaks) + #, alpha = perfSpace_alpha) +
    # to add outlines
    stat_contour(aes(z = .data$perfThreshComb), breaks = breaks, colour = "black") +
    labs(x = xyLabels[1], y = xyLabels[2]) +
    scale_x_continuous(expand=c(0, 0)) +    # no extra space on x and y axes
    scale_y_continuous(expand=c(0, 0)) +
    coord_cartesian(xlim=xlimits, ylim=ylimits) + 
    theme_heatPlot()
  
  # ******* set to NULL in case there is a column named "perfThreshComb" **** 
  # Remove if fill based on exceeded thresholds for climData points are to be implemented
  if (!is.null(climData[[colnames(plotData[3])]])) climData[[colnames(plotData[3])]] <- NULL
  
  p2List <- addClimData(p1, climData, colnames(plotData[3]), xyAtts, xlimits, ylimits)
  p2 <- p2List[[1]]

  
  if (is.null(colMap)) {
    coloursIn <- grDevices::adjustcolor(foreSIGHT.colmap(length(nThresh)), alpha.f = threshPlot_alpha)
  } else {
    if (length(colMap) < length(nThresh)) stop(paste0("col should contain atleast ", length(nThresh), " colours."))
    coloursIn <- colMap[1:length(nThresh)]
  }
  
  p2 <- p2 + scale_fill_manual(values = coloursIn, labels = nThresh) +
    guides(fill = guide_legend(title = "Number of Thresholds Exceeded", title.position = "bottom", title.hjust = 0.5, order = 1, override.aes = list(shape = NA))) + 
    # does not work for long legend titles: # https://stackoverflow.com/questions/48000292/center-align-legend-title-and-legend-keys-in-ggplot2-for-long-legend-titles
    # look for workaround if required
    theme(legend.title.align = 0.5) +
    #theme(legend.text = element_text(hjust = 0.5)) +
    labs(tag = tag_text)
    
  print(p2) 
  return(p2)
  
}


addClimData <- function(p1,                  # ggplot object to add the plot to
                        climData,            # climate data in a data.frame
                        perfName,            # name of the column to use for fill colour, if available
                        xyAtts,              # x & y attribute tags
                        xlimits,             # xlim (chop points outside this range)
                        ylimits,             # ylim ( " )
                        colLim = NULL,       # specified (colLimIn will be modified only if this is NULL)
                        colLimIn = NULL      # colour limits are checking with plotData
                        ) {
  # climate data
  if (!is.null(climData)) {
    
    # remove climate data points that fall outside the x-y limits
    outX <- which(climData[[xyAtts[1]]] < xlimits[1] | climData[[xyAtts[1]]] > xlimits[2])
    outY <- which(climData[[xyAtts[2]]] < ylimits[1] | climData[[xyAtts[2]]] > ylimits[2])
    if ((!identical(outX, integer(0))) | (!identical(outY, integer(0)))) {
      climData <- climData[-c(outX, outY), ]
    }
    
    # fix size based on the number of climData points
    climPoints <- nrow(climData)
    ptSize <- 4
    if (climPoints > 10) ptSize <- 3
    if (climPoints > 20) ptSize <- 2
    if (climPoints > 300) {
      ptSize <- 1.5
      perfSpace_climDataBg <- perfSpace_climDataBg2
    }
    
    if (is.null(climData$Name)) climData$Name <- "Climate Data"
    if (length(unique(climData$Name)) >= 2) {
      ncolLeg <- 2
    } else {
      ncolLeg <- length(unique(climData$Name))
    }
    
    if (sum(colnames(climData) %in% xyAtts) == 2) {
      if (climPoints == 0) {
        warning("climData is not plotted since the perturbations in the data are outside the ranges of the perturbations in sim.")
      } else {
        if (is.null(climData[[perfName]])) {
          p1 <-  p1 + geom_point(data = climData, mapping = aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], shape = .data$Name), show.legend = TRUE, size = ptSize, colour = perfSpace_climDataCol, fill = perfSpace_climDataBg) + 
            scale_shape_manual(name = NULL, values = rep(c(21, 22, 24, 25, 23, 16, 17, 18, 19, 20, c(0:14)), 20), guide = guide_legend(order = 2, ncol = ncolLeg))
        } else {
          # colour points based on performance value
          # currently this code is used only by plotPerformanceSpace - can implement for plotPerformanceThesh if required
          
          # expand the colLim to include the range of the climate data
          if (is.null(colLim)) {
            colLimIn[1] <- min(colLimIn[1], min(climData[[perfName]]))
            colLimIn[2] <- max(colLimIn[2], max(climData[[perfName]]))
          }
          
          # if there are more models than can be deficted using shape
          if (length(unique(climData$Name)) > 5) {
            climData$Name <- "Climate Data"
            message("Legend cannot differentiate between more than five unique climdata$Name. All data are plotted using the same shape.")
          }
          
          p1 <-  p1 + geom_point(data = climData, mapping = aes(x = .data[[xyAtts[1]]], y = .data[[xyAtts[2]]], fill = .data[[perfName]], shape = .data$Name), show.legend = TRUE, size = ptSize, colour = "black", alpha = perfSpace_alpha) + 
            scale_shape_manual(name = NULL, values = c(21, 22, 24, 25, 23), 
                               guide = guide_legend(overide.aes = list(fill = NA), order = 2, ncol = ncolLeg))
        }
      }
    } else {
      warning(paste0("climData is not plotted since it does not contain ", paste(xyAtts[(colnames(climData) %in% xyAtts)], sep = ","), "."))
    }
  }
  return(list(p1,
              colLimIn))
}
