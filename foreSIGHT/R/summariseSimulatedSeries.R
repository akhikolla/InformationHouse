#################################
##   SUMMARY PLOTTING WRAPPER  ##
#################################


#CONTAINS
  # plotSummary()
  # exposureSummary()
  
#NOT ALL VARS WILL NEED ALL PLOTS
# SOME JOINT PLOTTING WLL BE REQUIRED
#START WITH TRAFFIC LIGHT PLOTTING
#MORE DETAILED PLOTTING
#SOME TABLES?


#DAILY AMOUNTS TAKEN IN BRACKET
#metricTagList=c("dyAll","tot","totCorrel","dyCorrel","dyWet","nWet")
#plus - SOMETHING FOR EXTREMES OTHER THAN MAXIMA
#something for w/d-spells
  
# Anjana: this function does not appear to be used anywhere
# plotSummary<-function(obs=NULL,
#                       sim=NULL,
#                       simVar=NULL,
#                       datInd=NULL,
#                       attSel=NULL,
#                       attPrim=NULL,
#                       simTarget=NULL,
#                       target=NULL,
#                       targetType=NULL,
#                       modelTag=NULL,
#                       optimArgs=NULL
# ){
#   
# 
#   nTarget=nrow(target)  
#     
#   # DO EVERYTHING FOR OBSERVED SERIES ONCE
#   obsDat=list()
#   for(i in 1:length(simVar)){
#     plotVar=simVar[i]            #select variable to evaluate
#     obsDat[[plotVar]]=collateDat(TS=obs[[plotVar]],datInd=datInd[["obs"]],plotVar=plotVar)
#   }
#   
#   #LOOP OVER TARGETS
#   for(i in 1:nTarget){
#     simDat=list()
#     
#     #PLOT STUFF TO A PDF
#     #LABEL A PDF
#     fnam=nameMaker(attSel=attSel,target=target[i,])  #make filename
#     fnam=paste(fnam,".pdf",sep="")
#     pdf(file=fnam,height=8.27,width=11.69)   #landscape a4 page
#     par(mar=c(3,5,3,3),oma=c(2,2,2,2))
#     
#     #FRONT BOILERPLATE INFO
#     frontBoilerPlateInfo(modelTag=modelTag,
#                          targetLocn=target[i,],
#                          spot=i,
#                          nTarget=nTarget,
#                          attSel=attSel,
#                          attPrim=attPrim,
#                          optimArgs=optimArgs,
#                          sim=sim,
#                          simVar=simVar
#     )
#     
#     #TRAFFIC LIGHT PLOT HERE
#     if(modelTag[1] != "Simple-ann"){
#       trafficAttPlot(attSel=attSel,attPrim=attPrim,simPt=simPt[i,],target=target[i,],targetType=targetType)
#     }
#     
#     #ADD ADDITIONAL TRAFFIC LIGHT PANELS FOR REMAINING INFORMATION (SMALLER IN SIZE..., for eacg simVar set)
#     
#     
#     #SET LAYOUT - 2 ROWS, 1 COLUMN
#     par(mfrow=c(2,1),xaxs="i")
#     par(mar=c(3,5,3,3),oma=c(3,5,3,3),xpd=FALSE)
#     for(mod in 1:length(simVar)){
#       plotVar=simVar[mod] 
#       
#       switch(modelTag[1],
#         "Simple-ann" = {simTest=as.vector(unlist(sim[[i]][plotVar]))},
#                        {simTest=sim[[i]][[plotVar]]$sim}
#       )
#     
#       if(plotVar == "P"){mult=0.5}else{mult=1.05}
#       simTS.overlayMonthlyObsRange(obsDat=obsDat[[plotVar]],simTS=simTest,datInd=datInd[["obs"]],label=plotVar,range.mult=mult)
#     }
#   
#     for(mod in 1:length(simVar)){
#       
#       plotVar=simVar[mod]            #select variable to evaluate
#       switch(modelTag[1],
#              "Simple-ann" = {simTest=as.vector(unlist(sim[[i]][plotVar]))},
#                             {simTest=sim[[i]][[plotVar]]$sim}
#       )
#       
#       simDat[[plotVar]]=collateDat(TS=simTest,datInd=datInd[[modelTag[mod]]],plotVar=plotVar)
#       
#       #print(obsDat[[plotVar]][["ann_count_nWet"]])
#      # print(simDat[[plotVar]][["ann_count_nWet"]])
#       
#       #SET LAYOUT - 2 ROWS, 1 COLUMN
#       par(mfrow=c(2,1),xaxs="i",xpd=FALSE)
#       par(mar=c(3,5,3,3),oma=c(3,5,3,3))
#       #REGULAR PLOTS
#       #monthwise batch
#       runTag="mon_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
#       monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       runTag="mon_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
#       monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       
#       runTag="mon_sum_dyAll"; lab=paste(plotVar,": total", sep="")
#       monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       
#       if(plotVar == "P"){
#         runTag="mon_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
#         monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#         
#         runTag="mon_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
#         monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#         
#         # runTag="mon_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
#         # monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
# 
#         runTag="mon_count_nWet"; lab=paste(plotVar,": no. wet days", sep="")
#         monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       }
#   
#       
#       #seasonal batch
#       runTag="seas_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
#       seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       runTag="seas_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
#       seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       runTag="seas_sum_dyAll"; lab=paste(plotVar,": total", sep="")
#       seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       # # 
#       # if(plotVar == "P"){
#       # #   runTag="seas_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
#       # #   monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       # #   
#       # #   runTag="seas_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
#       # #   monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       # #   
#       #   runTag="seas_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
#       #   monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       # 
#       #   runTag="seas_sd_nWet"; lab=paste(plotVar,": no. wet days sd", sep="")
#       #   monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       # }
#       
#       #annual batch
#       #change 
#       par(mfrow=c(2,3),xaxs="i",xpd=FALSE)  #assuming a4 landscape layout
#       par(mar=c(3,5,3,3),oma=c(3,5,3,3))
#       
#       runTag="ann_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
#       annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       
#       runTag="ann_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
#       annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       
#       runTag="ann_sum_dyAll"; lab=paste(plotVar,": total", sep="")
#       annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#       
#       if(plotVar == "P"){
#         runTag="ann_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
#         annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#         runTag="ann_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
#         annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab) 
#         # runTag="ann_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
#         # annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#         runTag="ann_count_nWet"; lab=paste(plotVar,": no. wet days", sep="")
#         annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#         
#       }
#       
#     }
#     dev.off()  #STOP PLOTTING TO PDF
#   }
# 
# }

#--------------------------------------------------------
# Anjana: this function is not used anywhere
# exposureSummary<-function(targetMat=NULL,     # targetMat from the expSpace
#                           attSel=NULL,        
#                           paths=paths
#                           ){
#   fnam=paste0(paths$Plots,"/exposureSpace.pdf")
#   pdf(file=fnam,height=8.27,width=11.69)   #landscape a4 page
#   par(mar=c(3,5,3,3),oma=c(2,2,2,2))
#   exposureSlices(targetMat=targetMat,attSel=attSel)
#   dev.off()
# }

# Anjana: changed devPlotSummary to plotTargDiagnostics
#===========================================================================================
# Function that works on a single target - check if it is better to change input arguments
# There can be three functions
# 1. Function that works on the full simulation to create summary traffic light plots
# 2. Function that works on a specific target to create plots showing the span across replicates
# 3. Function that works on one single simulation (Rep-Targ) to create the same plots as function 2
# Functions 2&3 may be combined. This function would be similar to plotTargDiagnostics below.
# The different inputs to this function can be:
# A) The full simulation (i.e., list of reps - which is a list of targets): 
#                The function creates multiple PDF files, one for each target. The argument "writeToFile" has to be set to TRUE in this case, 
#                otherwise the user would end up with multiple plots in the console
# B) Simulation of a single target containing multiple replicates (include an example to show how to subset a target)
#                The function creates plots showing all the replicates; maybe use grey colour so that it shows up as a shading of replicates.
#                The plots are R plots by default. "writeToFile" can be set to TRUE to write the plots to a PDF with a front page information
# C) Simulation of a single target and replicate. The plots would contain only a single replicate
#                The function creates R plots by default. "writeToFile" can be set to TRUE to save plots in a PDF file.
#========================= NOTE =============================================================
# The base plotting code has to be the same for forward and inverse methods 
# So, ensure that the inverse elements are not propagated to the lowest level of the call tree
#============================================================================================

# Anjana: Creating heatmaps using "Diff" directly needs more work
# It would probably need to be split into two heatmaps plotted side-by-side with different cmap scales since the units of "Temp" are diff from the other variables
# It can be a model update for v1.1 (add another fn argument, performance = TRUE)
# This is done since it was decided to plot the biases directly instead of mapping them to performance

#' Creates summary plots of the biases in the scenarios
#'
#' \code{plotScenarios} uses a simulation performed using the function \code{generateScenarios} as input and creates heatmaps that show 
#' the biases in the simulated attributes with respect to the specified target values of the attributes. 
#' The plots show the magnitude (absolute value) of the mean biases, and the standard deviation of biases across replicates. The heatmaps can be used 
#' to evaluate how well the simulated attributes match the specified targets.
#' The biases are in units of percentage for attributes of variables like precipitation, and in units of degrees K for attributes of temperature. 
#' The function creates two heatmaps that show:
#' \itemize{
#' \item{magnitude of the mean biases across all the replicates}
#' \item{standard deviation of biases across all the replicates}
#' }
#' @param sim a list; contains a stochastic simulation or the summary of a stochastic simulation created using the function \code{generateScenarios}
#' @param simName a string; defaults to \code{NULL}). User-specified name of the simulation that will used as the heading in the 
#' saved pdf file to identify the simulation later. If \code{simName} is \code{NULL}, a random name will be assigned for the simulation. 
#' @param writeToFile logical; defaults to \code{FALSE}. Specifies whether the plots should be saved to a pdf file.
#' If set to true, the heatmaps will be saved to a pdf file that would also contain summary pages that show the attributes, models, and optimisation settings used to create \code{sim}.
#' @param fileName a string; defaults to \code{"plotScenarios.pdf"}. Specifies the name of the pdf file to be written, if the file exists it will be overwritten.
#' @param colMapRange a string; may be set to \code{"default"} or \code{"full"}. 
#' The argument specifies whether the colormap of the heatmap should span the full range of the data.
#' If set to \code{"default"}, the colourmap limits of attributes that are in units of percentage is set to 0\% to 10\%, 
#' and the colourmap limits of the attributes of temperature is set to 0 degrees K to 1 degrees K. 
#' If set to \code{"full"}, the colourmap limits are set to the minimum and maximum values in the data.
#' @details The argument \code{sim} may be a full stochastic simulation generated using the function \code{generateScenarrios} or the summary of the stochastic simulation
#' generated using \code{getSimSummary}
#' @return The function returns two R plots showing the biases in the targets of the scenarios generated using the function \code{generateScenarios}. 
#' The figures may be saved to a pdf file by setting the \code{writeToFile} argument to \code{TRUE}.
#' @seealso \code{createExpSpace}, \code{generateScenarions}, \code{getSimSummary}
#' @examples
#' \dontrun{
#' # the examples are nnot run since the run times are too long for CRAN
#' # create an exposure space
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m")
#' attHold <- c("P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m", 
#' "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(2, 1, 1)
#' attPerturbMin = c(0.9, 1, 1)
#' attPerturbMax = c(1.1, 1, 1)
#' expSpace <- createExpSpace(attPerturb = attPerturb, 
#'                            attPerturbSamp = attPerturbSamp, 
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType,
#'                            attHold = attHold,
#'                            attTargetsFile = NULL)
#' # load example data available in foreSIGHT
#' data(tankDat)
#' # perform stochastic simulation
#' sim <- generateScenarios(reference = tank_obs,
#'                          expSpace = expSpace,
#'                          simLengthNyrs = 30,
#'                          numReplicates = 2)
#' # plots heatmaps showing biases in simulated targets
#' plotScenarios(sim)
#' # to save the figures to a pdf file set writeToFile = TRUE
#' # using an example stochastic simulation summary provided with the package
#' data("egSimSummary")
#' plotScenarios(egSimSummary)
#' }
#' @export

plotScenarios <- function(sim,
                          simName = NULL,
                          writeToFile = FALSE,
                          fileName = "plotScenarios.pdf",
                          colMapRange = "default"
                          ) {
  
  if (is.null(sim[["controlFile"]])) {
      cat("controlFile is missing in the simulation. Are the scenarios generated using simple scaling?\n")
      stop("plotScenarios cannot be used on simple scaled data")
  }
  
  # get mean diff (set to diff not classified performance) and SD across replicates
  simTraffic <- getSimTraffic(sim)

  dataField <- c("mean", "SD")
  plotName <- c("Mean", "SD")
  
  plots <- list()
  # variable groups (Temp, other Vars)
  vGroups <- length(simTraffic[["mean"]])
   
  for (f in c(1,2)) {

    dataV <- NULL
    for (v in 1:vGroups) {
      dataIn <- simTraffic[[dataField[f]]][[v]]
      x <- colnames(dataIn)
      y <- rownames(dataIn)
      
      data <- expand.grid(Attribute=x, Target=y)
      data[[plotName[f]]] <- NA
      for (i in 1:nrow(data)) {
        colN <- as.character(data[i, 1])
        rowN <- as.character(data[i, 2])
        data[i, 3] <- abs(dataIn[rowN, colN])
      }
    dataV[[v]] <- data  
    }
    
    plots[[dataField[f]]] <- plotTrafficHeatmap(dataV, 
                       plotName[f],
                       attNameList = simTraffic[["attName"]],
                       markPrimList = simTraffic[["markPrim"]],
                       perturbedList = simTraffic[["perturbed"]],
                       writeToFile,
                       colMapRange)
  }

  # no of pages
  nPg <- length(plots[[1]])
  nTarg <- dim(simTraffic[["mean"]][[1]])[1]
  nAtt <- length(unlist(simTraffic[["attName"]])) # not completely accurate
  
  if (writeToFile) {
    if (file.exists(fileName)) {
      cat(paste0("\nThe file ", fileName, " will be overwritten.\n"))
      file.remove(fileName)
    }
    #pdf(file = fileName, width = 8.27, height = 11.69, paper = "a4")
    pdf(file = fileName, width = 11.69, height = 8.27, paper = "a4r")
    frontPageScenarios(sim, simName)
    advancedPageScenarios(sim)

    # *****
    # ADD function to add a table of targets here, i.e., what does Target1 mean 
    # *****
    
    for (i in 1:nPg) {
      # if (nTarg < nAtt) {
      #   nrow = NULL
      #   ncol = length(plots[["mean"]][[i]])
      # } else {
        nrow = length(plots[["mean"]][[i]])
        ncol = NULL
        if (nrow > 1) {
          prop <- (dim(simTraffic[["mean"]][[1]])[2])/(dim(simTraffic[["mean"]][[2]])[2])
          if (prop >= 4) {
            multiplier <- 0.5 - ((prop - 4)*0.06)
            if (multiplier < 0) multiplier <- 0.1
          } else if (prop == 1) {
            multiplier <- 0.975
          } else {
            multiplier <- 1 - (0.135*prop)
          }
          rel_heights <- c(prop*multiplier, 1)
        } else {
          rel_heights <- 1
        }
      # }
      print(cowplot::plot_grid(plotlist = plots[["mean"]][[i]],
                         align = "v",
                         nrow = nrow,
                         ncol = ncol,
                         rel_heights = rel_heights))
      # print(plots[["mean"]][[i]])
    }
    for (i in 1:nPg) {
      print(cowplot::plot_grid(plotlist = plots[["SD"]][[i]],
                         align = "v",
                         nrow = nrow,
                         ncol = ncol,
                         rel_heights = rel_heights))
                         #rel_heights = c(8,2)
      # print(plots[["SD"]][[i]])
    }
    dev.off()
    cat(paste0("\nFigures are saved to file: ", fileName, "."))
  } else {
    if (dim(simTraffic[["mean"]][[1]])[1] > 100) cat("The scenarios may contain too many targets to be examined in an R plot. Please call plotScenarios with writeToFile = TRUE to save the figures in a pdf file.")
    # print(plots[["mean"]][[1]])
    # print(plots[["SD"]][[1]])
    # if (nTarg < nAtt) {
    #   nrow = NULL
    #   ncol = length(plots[["mean"]][[1]])
    # } else {
      nrow = length(plots[["mean"]][[1]])
      ncol = NULL
    # }
      
      if (nrow > 1) {
        prop <- (dim(simTraffic[["mean"]][[1]])[2])/(dim(simTraffic[["mean"]][[2]])[2])
        if (prop >= 4) {
          multiplier <- 0.5 - ((prop - 4)*0.06)
          if (multiplier < 0) multiplier <- 0.1
        } else if (prop == 1) {
          multiplier <- 0.975
        } else {
          multiplier <- 1 - (0.135*prop)
        }
        rel_heights <- c(prop*multiplier, 1)
      } else {
        rel_heights <- 1
      }
      
    print(cowplot::plot_grid(plotlist = plots[["SD"]][[1]],
                       align = "v",
                       nrow = nrow,
                       ncol = ncol,
                       rel_heights = rel_heights))
    print(cowplot::plot_grid(plotlist = plots[["mean"]][[1]],
                       #axis = "l",
                       align = "v",
                       nrow = nrow,
                       ncol = ncol,
                       rel_heights = rel_heights))
  }
  #return(invisible())
  return(invisible(plots))
  }




# function to create the ggplot list objects using geom_tile - will be called by plotScenarios
#' @importFrom rlang .data
plotTrafficHeatmap <- function(targAttList,
                               field,
                               attNameList,
                               markPrimList,
                               perturbedList,
                               writeToFile = FALSE,
                               colMapRange = "default"
                               ) {
  # 1 or 2, tempeparture & other variables
  nMat <- length(targAttList)
  nTarg <- length(unique(targAttList[[1]]$Target))
  
  # xaxis label for penalty attributes
  if (any(unlist(markPrimList) == "*")) {
    penaltyLabel <- "*indicates penalty attributes"
  } else {
    penaltyLabel <- "There are no penalty attributes"
  }
  
  # if (field == "Variance") penaltyLabel <- paste0(penaltyLabel, "\ngrey areas have variances lower than 1")
  
  # sum of all attributes
  nAttTotal <- 0
  attDef <- NULL
  nPrim <- NULL
  attVarInd <- NULL
  nAtt <- NULL
  lTitle <- NULL
  colLimits <- NULL
  for (m in 1:nMat) {
    
    targAttMatrix <- targAttList[[m]]
    markPrim <- markPrimList[[m]]
    if (m ==1 ) {
      if (colMapRange == "full") {
        colLimits[[m]] <- c(min(targAttMatrix[,3]), max(targAttMatrix[,3]))
      #if (!is.na(max(targAttMatrix[,3])) & max(targAttMatrix[,3]) < 10) {
      } else {
        colLimits[[m]] <- c(0, trafficLim[["pc.lim"]][2])
      }
    } else {
      if (colMapRange == "full") {
      # if (!is.na(max(targAttMatrix[,3])) & max(targAttMatrix[,3]) < 0.1) {
        colLimits[[m]] <- c(min(targAttMatrix[,3]), max(targAttMatrix[,3]))
        # colLimits[[m]] <- c(0, 0.1)
      } else {
        colLimits[[m]] <- c(0, trafficLim[["diff.lim"]][2])
        #colLimits[[m]] <- c(min(targAttMatrix[,3]), max(targAttMatrix[,3]))
      }
    }
    
    heldFlag <- rep("H", length(perturbedList[[m]]))
    heldFlag[perturbedList[[m]]] <- "P"
    
    nAttTotal <- nAttTotal + length(unique(targAttMatrix$Attribute))
    nAtt[[m]] <- length(unique(targAttMatrix$Attribute))
  
    # Create labels using full attribute names
    attDef[[m]] <- paste0(mapply(tagBlender_noUnits, attNameList[[m]], USE.NAMES = FALSE), " (", heldFlag, ")")
    nPrim[[m]] <- length(which(markPrim=="*"))
  
    # Identifying points to add hlines betweeen variables
    attVarAll <- unlist(lapply(strsplit(attNameList[[m]], "_"), `[[`, 1))
    attVar <- unique(attVarAll)
    attVarInd_temp <- NULL
    if(length(attVar) > 1) {
      for (i in 2:length(attVar)) {
        tempInd <- which(attVarAll[(nPrim[[m]]+1):length(attVarAll)] == attVar[i])
        attVarInd_temp[i-1] <- (tempInd[1]+nPrim[[m]])
      }
    }
    attVarInd[[m]] <- attVarInd_temp
    if (attVar[1] == "Temp") {
      lTitle[[m]] <- legendTitleTraffic[["TempType"]]
    } else  {
      lTitle[[m]] <- legendTitleTraffic[["PType"]]
    }
    
  }
    
  # change size of the text based on the size of the matrix
  if (nTarg > 30 | nAttTotal > 40) {
    traffic_textSize = 11
    traffic_margins <- traffic_tightMargins
  }
  # if (nAtt > 80 | nTarg > 80) {
  #   traffic_margins <- traffic_tightMargins
  #   traffic_textSize = 7
  # }
  
  # if there are more than 50 targets - the plots need to be split between pages
  nTPg <- 50   # no. of targets per page
  targAttMatrixList <- list()
  aspectRatio <- list()
  if (writeToFile & nTarg > 50) {
    nPg <- ceiling(nTarg/nTPg)
    # recalculate targets per page for almost equal division
    nTPg <- ceiling(nTarg/nPg)
    allTarg <- unique(targAttList[[1]]$Target)
    aspectRatio <- replicate(nPg, vector("list", nMat), simplify = FALSE)
    targAttMatrixList <- replicate(nPg, vector("list", nMat), simplify = FALSE)
    for (i in 1:nPg) {
      iTargs <- allTarg[(1+(i-1)*nTPg):(i*nTPg)]
      # loop over Temp & other var
      for (m in 1:nMat) {
        indData <- which(targAttList[[m]]$Target %in% iTargs)
        targAttMatrixList[[i]][[m]] <- targAttList[[m]][indData, ]
        aspectRatio[[i]][[m]] = sum(!is.na(iTargs))/nAtt[[m]]
      }
    }
  } else {
    nPg <- 1
    aspectRatio <- replicate(nPg, vector("list", nMat), simplify = FALSE)
    targAttMatrixList[[1]] <- targAttList
    for (m in 1:nMat) {
      aspectRatio[[1]][[m]] = nTarg/(nAtt[[m]])
    }
  }

  p1 <- replicate(nPg, vector("list", nMat), simplify = FALSE)
  p2 <- replicate(nPg, vector("list", nMat), simplify = FALSE)
  # add vline to separate penalty attributes, and atts of different variables
  yInt <- NULL
  for (m in 1:nMat) {
    # yInt[[m]] <- factor(c(attNameList[[m]][nPrim[[m]]+1], attNameList[[m]][attVarInd[[m]]]),
    #               levels=levels(targAttList[[m]]$Attribute))
    # add lines only between atts of different variables
    yInt[[m]] <- factor(c(attNameList[[m]][attVarInd[[m]]]),
                        levels=levels(targAttList[[m]]$Attribute))
  }
  
  for (i in 1:nPg) {  
      for (m in 1:nMat){
        if(field == "Mean") {
          p1[[i]][[m]] <- ggplot(targAttMatrixList[[i]][[m]], aes(y = factor(.data$Attribute), x = .data$Target)) +
          geom_tile(aes(fill = .data$Mean), colour = traffic_tileOutline) + 
          scale_fill_gradientn(colours = traffic.col, limits = colLimits[[m]], 
                               oob = scales::squish,
                               guide = guide_colorbar(title = lTitle[[m]], title.position = "right", barwidth = 12, barheight = 0.5))
        } else {
          p1[[i]][[m]] <- ggplot(targAttMatrixList[[i]][[m]], aes(y = factor(.data$Attribute), x = .data$Target)) +
          geom_tile(aes(fill = .data$SD), colour = traffic_tileOutline) + 
          scale_fill_gradientn(colours = traffic.col, limits = colLimits[[m]],
                               oob = scales::squish,
                               guide = guide_colorbar(title = lTitle[[m]], title.position = "right", barwidth = 12, barheight = 0.5))
        } 
      #, limits = c(1, 3), breaks = c(1, 2, 3), labels = c("Good", "Fair", "Poor"))
      #scale_fill_gradientn(colours = traffic.col, limits = c(1, 3), breaks = c(1, 2, 3), labels = c("Good", "Fair", "Poor")) #+
      #geom_text(data = labeltext, aes(label = label), vjust = 0.5, hjust = textHjust, angle = 90, color = textCol, size = traffic_textSize*0.352777778)
    # } else {
    #   p1[[i]] <- ggplot(targAttMatrixList[[i]], aes(y = factor(Attribute), x = Target)) +
    #   geom_tile(aes(fill = Variance), colour = traffic_tileOutline) + 
    #   # Need to change this
    #   scale_fill_gradientn(colours = traffic.col, limits = c(1, 3), breaks = c(1, 2, 3), labels = c("1", "2", "3"))
    # }
  
    p1[[i]][[m]] <- p1[[i]][[m]] + geom_hline(yintercept = (as.integer(yInt[[m]])-0.5))
      }
    
    if (nMat == 2) {
      p2[[i]][[1]] <- p1[[i]][[1]] + theme_traffic_upper(traffic_textSize) +
        theme(plot.margin = unit(traffic_upperMarg,"cm")) +
        labs(title = plotTitleTraffic[[field]], y = "") + 
        scale_y_discrete(breaks = paste0(attNameList[[1]]), labels = paste0(attDef[[1]], markPrimList[[1]])) +
        scale_x_discrete(position = "top") + 
        theme(aspect.ratio = 1/aspectRatio[[i]][[1]])
      
      p2[[i]][[2]] <- p1[[i]][[2]] + theme_traffic_lower(traffic_textSize) +
        theme(plot.margin = unit(traffic_lowerMarg,"cm")) +
        labs(x = penaltyLabel, y = "") + 
        scale_y_discrete(breaks = paste0(attNameList[[2]]), labels = paste0(attDef[[2]], markPrimList[[2]])) +
        theme(aspect.ratio = 1/aspectRatio[[i]][[2]]) + labs(tag = tag_text)
    } else {
      p2[[i]][[1]] <- p1[[i]][[1]] + theme_traffic(traffic_textSize) +
        theme(plot.margin = unit(traffic_margins,"cm")) +
        labs(x = penaltyLabel, y = "", title = plotTitleTraffic[[field]]) + 
        scale_y_discrete(breaks = paste0(attNameList[[1]]), labels = paste0(attDef[[1]], markPrimList[[1]])) +
        theme(aspect.ratio = 1/aspectRatio[[i]][[1]]) + labs(tag = tag_text)
    }
   }
  
  return(p2)
}


# can have a targetNum (and RepNum) argument to point to the target to be plotted - set to "1" and "NULL" by default
# If targetNum is set to NULL all targets are to be plotted. In this case, write to file has to be true

plotTarget <- function(obs,                    # observations
                       sim,                    # the full simulation
                       targetNum,              # supply the targetNum, default set to Target1
                       repNum = 1,             # the default should be NULL & create plots for all replicates, setting to one for now so that the original function works
                       writeToFile = TRUE,
                       fileName = "plotTarget.pdf"
                       ) {
                            
  # Get required Information - Used by the original function
  #==================================================================
  # Unpack data from sim
  attSel <- c(sim[["expSpace"]][["attPerturb"]], sim[["expSpace"]][["attHold"]])
  attPrim <- sim[["controlFile"]][["penaltyAttributes"]]
  optimArgs <- sim[["controlFile"]][["optimisationArguments"]]
  optimArgs[["lambda.mult"]] <- sim[["controlFile"]][["penaltyWeights"]]
  nml <- sim[["controlFile"]]
  
  # Variable and target type
  attSel_varType <- vapply(attSel, FUN = get.attribute.varType, FUN.VALUE = character(1), USE.NAMES = FALSE)
  targetType <- vapply(attSel_varType, FUN = get.target.type, FUN.VALUE = character(1), USE.NAMES = FALSE)
  
  # modelTag and variables
  simVar <- names(nml[["modelType"]])
  modelTag <- NULL
  for (v in simVar) {
    modelTag <- c(modelTag, getModelTag(nml = nml, v))
  }
  
  # Date Indices
  datInd <- list()
  datInd[["obs"]] <- get.date.ind(dd = obs$day,
                                  mm = obs$month,
                                  yy = obs$year,
                                  nperiod = 12,
                                  southHemi = TRUE)
  datInd[["sim"]] <- get.date.ind(dd = sim$simDates$day,
                                  mm = sim$simDates$month,
                                  yy = sim$simDates$year,
                                  nperiod = 12,
                                  southHemi = TRUE)
  # Target to be plotted
  target <- sim[["expSpace"]][["targetMat"]][targetNum, ]
  
  tarName <- paste0("Target", targetNum)
  repName <- paste0("Rep", repNum)
  simPt <- sim[[repName]][[tarName]][["targetSim"]]
  simTarget <- sim[[repName]][[tarName]]
  
  #=======================================================================

  # DO EVERYTHING FOR OBSERVED SERIES ONCE
  obsDat=list()
  for(i in 1:length(simVar)){
    plotVar=simVar[i]            #select variable to evaluate
    obsDat[[plotVar]]=collateDat(TS=obs[[plotVar]],datInd=datInd[["obs"]],plotVar=plotVar)
  }
  
  #LOOP OVER TARGETS
  
    simDat=list()
    
    if (writeToFile) {
      #PLOT STUFF TO A PDF
      pdf(file=fileName,height=8.27,width=11.69)   #landscape a4 page
    }
    par(mar=c(3,5,3,3),oma=c(2,2,2,2))
    
    # Anjana: Need to modify this for the updated model settings
    #===========================================================
    #FRONT BOILERPLATE INFO
    # frontBoilerPlateInfo(modelTag=modelTag,
    #                      targetLocn=target,
    #                      spot=id,
    #                      nTarget=nTarget,
    #                      attSel=attSel,
    #                      attPrim=attPrim,
    #                      optimArgs=optimArgs,
    #                      sim=sim,
    #                      simVar=simVar
    # )
    #============================================================
    
    #TRAFFIC LIGHT PLOT HERE
    if(modelTag[1] != "Simple-ann"){
      trafficAttPlot(attSel=attSel,attPrim=attPrim,simPt=simPt,target=target,targetType=targetType)
    }
    
    
    #SET LAYOUT - 2 ROWS, 1 COLUMN
    par(mfrow=c(2,1),xaxs="i")
    par(mar=c(3,5,3,3),oma=c(3,5,3,3),xpd=FALSE)
    for(mod in 1:length(simVar)){
      plotVar=simVar[mod]

      switch(modelTag[1],
             "Simple-ann" = {simTest=as.vector(unlist(simTarget[plotVar]))},
             {simTest=simTarget[[plotVar]]$sim}
      )

      if(plotVar == "P"){mult=0.5}else{mult=1.05}
      simVobsTS(simTS=simTest[1:1460],obsTS=obs[[plotVar]][1:1460],datInd=datInd[["obs"]],varName=plotVar,asRollAv=TRUE) #datInd not used
      
      simTS.overlayMonthlyObsRange(obsDat=obsDat[[plotVar]],simTS=simTest,datInd=datInd[["obs"]],label=plotVar,range.mult=mult)
    }

    for(mod in 1:length(simVar)){

      plotVar=simVar[mod]            #select variable to evaluate
      switch(modelTag[1],
             "Simple-ann" = {simTest=as.vector(unlist(simTarget[plotVar]))},
             {simTest=simTarget[[plotVar]]$sim}
      )

      simDat[[plotVar]]=collateDat(TS=simTest,datInd=datInd[["sim"]],plotVar=plotVar)

      #print(obsDat[[plotVar]][["ann_count_nWet"]])
      # print(simDat[[plotVar]][["ann_count_nWet"]])

      #SET LAYOUT - 2 ROWS, 1 COLUMN
      par(mfrow=c(2,1),xaxs="i",xpd=FALSE)
      par(mar=c(3,5,3,3),oma=c(3,5,3,3))
      #REGULAR PLOTS
      #monthwise batch
      runTag="mon_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
      monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      runTag="mon_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
      monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      runTag="mon_sum_dyAll"; lab=paste(plotVar,": total", sep="")
      monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      if(plotVar == "P"){
        runTag="mon_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
        monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

        runTag="mon_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
        monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

        # runTag="mon_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
        # monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

        runTag="mon_count_nWet"; lab=paste(plotVar,": no. wet days", sep="")
        monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
      }


      #seasonal batch
      runTag="seas_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
      seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      runTag="seas_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
      seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      runTag="seas_sum_dyAll"; lab=paste(plotVar,": total", sep="")
      seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      #annual batch
      #change
      par(mfrow=c(2,3),xaxs="i",xpd=FALSE)  #assuming a4 landscape layout
      par(mar=c(3,5,3,3),oma=c(3,5,3,3))

      runTag="ann_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
      annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      runTag="ann_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
      annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      runTag="ann_sum_dyAll"; lab=paste(plotVar,": total", sep="")
      annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      if(plotVar == "P"){
        runTag="ann_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
        annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
        runTag="ann_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
        annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
        # runTag="ann_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
        # annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
        runTag="ann_count_nWet"; lab=paste(plotVar,": no. wet days", sep="")
        annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)

      }

    }
    if (writeToFile) dev.off()  #STOP PLOTTING TO PDF
  
  
}


