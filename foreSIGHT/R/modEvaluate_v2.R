modEvaluate <- function(data,                                         # observations in a data.frame
                        sim_path,                                     # path to simulation directory, path should end with '/'
                        sim_prefix,                                   # prefix of simulation filenames
                        nreps,                                        # number of replicates
                        write_pdf = FALSE,                            # write to pdf? T/F
                        out_path = ".",                               # output directory
                        out_filename = "modEvaluate_output.pdf"       # filename of the output pdf
) {
  
  # Read simulated data
  simData <- read_simData(path = sim_path, prefix = sim_prefix, nreps = nreps)
  
  ################################################
  # data AND simData CHECKs Go HERE - to be added
  ################################################
  
  # identify common variables
  obsVars <- names(data)[-c(1, 2, 3)]
  simVars <- names(simData[[1]])[-c(1, 2, 3)]
  commonVars <- intersect(obsVars, simVars)
  
  # process data & plot
  qqplot_list <- qqplot_master(data, simData, variables = commonVars)
  
  # Rplot or print plots to pdf
  if (write_pdf == TRUE) {
    pdf(paste0(out_path, "\\", out_filename), width = 12, height = 10)
    for (i in seq(1, length(qqplot_list))) {
      print(qqplot_list[[i]])
    }
    dev.off()
    return(invisible())
  } else {
    return(qqplot_list)
  }
}

# Functions
#-------------------------------------------------------------------------------------------------------------------

# Read sim from files and save in a list
#----------------------------------------------
read_simData <- function(path,
                         prefix,
                         nreps
){
  simData <- list()
  
  for (i_rep in 1:nreps) {
    
    file <- paste0(path, "\\", prefix, i_rep, ".csv")
    if (!file.exists(file)) { 
      stop( paste0("The Simulation data file \"", file, "\" does not exist."), call. = FALSE )
    }
    else {
      simRep <- read.csv(file, header = TRUE)
      simData[[i_rep]] <- simRep
      rm(simRep)
    }
    
  }
  return(simData)
  
}


# qqplot master function
# uses three main fns: qqplot_masterList
#                      qqplot_procData
#                      qqplot_createPlot
#---------------------------------------------

qqplot_master <- function(obs,
                          simData,
                          variables
) {
  
  # Read in the master lists
  qqplot_masterList <- qqplot_masterList()
  
  units <- qqplot_masterList$units
  calcFrom <- qqplot_masterList$calcFrom
  calc_FUN_list <- qqplot_masterList$calc_FUN_list
  agg_FUN_list <- qqplot_masterList$agg_FUN_list
  cdf_FUN_list <- qqplot_masterList$cdf_FUN_list
  qqplot_labels <- qqplot_masterList$qqplot_labels
  
  # add variables calculated from P/PET/Temp/Radn to variables
  for (i_var in c("P", "PET", "Temp", "Radn")) {
    if (i_var %in% variables) {
      variables <- c(variables, names(calcFrom[calcFrom == i_var]))
    }
  } 
  
  # loop to create plots for all variables
  qqplot_list <- list()
  plot_count <- 1
  
  for (i_var in variables) {
    
    # Variable to be used for aggregation & sort
    if (!is.null(calcFrom[[i_var]])){
      varname <- calcFrom[[i_var]]
    } else {
      varname <- i_var
    }
    
    # Process data
    procData <- qqplot_procData(varname = varname,
                                obs = obs,
                                simData = simData,
                                calc_FUN = calc_FUN_list[[i_var]],
                                agg_FUN = agg_FUN_list[[i_var]] )
    
    ylabel <- paste0(qqplot_labels[[i_var]], " (", units[[i_var]], ")" )
    
    # Create plot  
    qqplot_list[[plot_count]] <- qqplot_createPlot(procData = procData,
                                                   cdf_FUN = cdf_FUN_list[[i_var]],
                                                   ylabel = ylabel,
                                                   ylim = NULL )
    plot_count <- plot_count + 1
  }
  
  names(qqplot_list) <- variables
  return(qqplot_list)
}


# Process Data for qqplot
#------------------------------------------------
qqplot_procData <- function(varname,
                            obs,
                            simData,
                            calc_FUN,
                            agg_FUN
) {
  
  # identify columns in obs & sim
  col_obs <- which(names(obs) == varname)
  col_sim <- which(names(simData[[1]]) == varname)
  
  # subset obs
  if (is.null(calc_FUN)){
    obsDataByYr <- data.frame(as.numeric(obs$year), obs[col_obs])
  }
  else {
    obsDataByYr <- data.frame(as.numeric(obs$year), calc_FUN(obs[col_obs]))
  }
  names(obsDataByYr) <- c("year", "data") 
  
  # subset sim
  simDataByYr <- simData_getVar(simData = simData, col = col_sim, indexvar = "year", calc_FUN = calc_FUN)
  
  # Aggregate and sort data
  procData <- agg_sort_obsSim(obsDataByYr, simDataByYr, agg_FUN = agg_FUN)
  
  return(procData)
}



# Unpack data for plot, create x-axis based on the CDF
# call plotting function
#--------------------------------------------------
qqplot_createPlot <- function(procData,
                              cdf_FUN,
                              ylabel,
                              ylim
) {
  
  # Assemble data for plot
  # xaxis based on plotting position, CDF from cdf_FUN
  nyrs <- length((procData$obs))
  probs <- seq(1, nyrs) / (nyrs + 1)
  quants <- cdf_FUN(probs)
  
  plotdata <- data.frame(quants, procData$obs, procData$simLower, procData$simMedian, procData$simUpper)
  names(plotdata)<-c("quants", names(procData))
  
  xtick_labels <- c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99)
  xticks <- cdf_FUN(xtick_labels)
  
  p <- problim.plot(plotdata = plotdata, ylim = ylim, ylabel = ylabel, xticks = xticks, xtick_labels = xtick_labels)
  
}



# Subset and calculate variable from reps of sim
#------------------------------------------------
simData_getVar <- function(simData,    # list containing replicates. Each element is a data.frame containing all the variables
                           col,        # column index (in the data.frame) of the variable to be returned
                           indexvar,   # the name of the column to return as the indexing variable
                           calc_FUN    # calculation function to be applied on the variable
) {
  
  
  calc_data_bycol <- function(data, col, indexvar, calc_FUN) {
    if (is.null(calc_FUN)) {
      data_col <- data.frame(data[[indexvar]], data[col])
    } else {
      data_col <- data.frame(data[[indexvar]], calc_FUN(data[col]))
    }
    names(data_col) <- c(indexvar, "data")
    return(data_col)
  }
  
  simData_wIndex <- lapply(simData, FUN = calc_data_bycol, col = col, indexvar = indexvar, calc_FUN = calc_FUN)
  return(simData_wIndex)
}



# Aggregate and sort observations and simulations
#------------------------------------------------
agg_sort_obsSim <- function(obsDataByYr,
                            simDataByYr,
                            agg_FUN
) {
  
  nyrs <- length(unique(obsDataByYr$year))
  nreps <- length(simDataByYr)
  
  # Obs
  obsSort <- agg_sort(obsDataByYr, agg_FUN = agg_FUN)
  
  # All reps of Sim
  simSort_mat <- sapply(simDataByYr, agg_sort, agg_FUN = agg_FUN)
  simSort_vec <- as.vector(simSort_mat)
  
  quant <- rep(seq(1, nyrs), nreps)
  simSort <- data.frame(cbind(quant, simSort_vec))
  names(simSort) <- c("quant", "data")
  
  # Calc Stats of Sim
  simLower <- aggregate(x = simSort$data, by = list(simSort$quant), FUN = quantile, probs = 0.025)[[2]]
  simMedian <- aggregate(x = simSort$data, by = list(simSort$quant), FUN = median)[[2]]
  simUpper <- aggregate(x = simSort$data, by = list(simSort$quant), FUN = quantile, probs = 0.975)[[2]]
  
  return(list("obs" = obsSort,
              "simLower" = simLower,
              "simMedian" = simMedian,
              "simUpper" = simUpper)
  )
  
}




# Smaller utility functions
#-------------------------------------------------

agg_sort <- function(data, agg_FUN) {
  
  dataSummary <- aggregate(x = data$data, by = list(data$year), FUN = agg_FUN$FUN, na.rm = agg_FUN$na.rm)[2]
  return(sort(dataSummary[, 1]))
  
}

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
ma2 <- ma
formals(ma2)$n <- 2

wetdays <- function(x) {
  x[x > 0] <- 1
  return(x)
}


# Core Plotting Functions
#-----------------------------------------------------
#' @importFrom rlang .data
problim.plot <- function(plotdata, ylim, ylabel, xticks, xtick_labels) {
  
  theme_set(cowplot::theme_cowplot())
  # shapes
  shp_obs <- 19  # dot
  shp_sim <- 3   # +
  
  # colours
  col_obs <- "red"
  col_sim <- "black"
  
  p <- ggplot(data = plotdata)+
    
    # plots
    geom_point(aes(x = .data$quants, y = .data$obs, colour = "Observed"), shape = shp_obs)+
    geom_point(aes(x = .data$quants, y = .data$simMedian, colour = "Sim. Median"), shape = shp_sim)+
    geom_line(aes(x = .data$quants, y = .data$simUpper, colour = "Sim. 95% Limits"), linetype = "dashed")+
    geom_line(aes(x = .data$quants, y = .data$simLower, colour = "Sim. 95% Limits"), linetype = "dashed")+
    
    # xaxis & legend
    scale_x_continuous(breaks = xticks, labels = xtick_labels)+
    scale_colour_manual(name = "", labels = c("Observed", "Sim. Median", "Sim. 95% Limits"), values = c(col_obs, col_sim, col_sim),
                        guide = guide_legend(override.aes = list(linetype = c("blank", "blank", "dashed"), shape = c(shp_obs, shp_sim, NA))))+
    
    # layout
    theme(panel.border = element_rect(size = , colour = "black"),
          panel.grid.major.x = element_line(linetype = "dashed", colour = "grey"),
          panel.grid.major.y = element_line(linetype = "dashed", colour = "grey"),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          legend.position = "bottom", 
          legend.justification = "centre" )+
    
    # axis labels
    labs(x = "Cumulative Probability",
         y = ylabel) +
    
    # use ylim if specified
    if(!is.null(ylim)) { scale_y_continuous(limits = c(ylim[1], ylim[2]), breaks = c(seq(ylim[1], ylim[2], ylim[3]))) }
  
}

