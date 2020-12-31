#' Plots the location of points in a two-dimensional exposure space
#' 
#' The function uses an exposure space created using the function \code{createExpSpace} as input and creates a plot of the 
#' two dimensional (2D) exposure space. \code{plotExpSpace} plots only 2D spaces consisting of samples of 2 attributes.
#'
#' @param expSpace list; an exposure space created using the function \code{createExpSpace}
#' @param y a string; tag of a perturbed attribute to plot on the y-axis. Defaults to \code{expSpace[["attPerturb"]][1]}.
#' @param x a string; tag of a perturbed attribute to plot on the x-axis. Defaults to \code{expSpace[["attPerturb"]][2]}.
#' @details The number of dimensions of an exposure space is equal to the number of perturbed attributes in that space. If the
#' exposure space has more than 2 dimensions (perturbed attributes), this function can be used to plot 2D slices of the space.
#' Note that the default arguments of this function is defined to plot a slice showing the first two dimensions of the space, 
#' arguments \code{x} and \code{y} may be specified to plot alternate dimensions. 
#' @seealso \code{createExpSpace}
#' @export
#' @examples
#' # create an exposure space that has more than 2 dimensions
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_Feb_tot_m")
#' attHold <- c("P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m", "P_DJF_avgDSD_m", 
#' "Temp_ann_rng_m", "Temp_ann_avg_m")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(5, 5, 5)
#' attPerturbMin = c(0.8, 0.9, 0.85)
#' attPerturbMax = c(1, 1.1, 1.05)
#' expSpace <- createExpSpace(attPerturb = attPerturb, 
#'                            attPerturbSamp = attPerturbSamp, 
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType,
#'                            attHold = attHold,
#'                            attTargetsFile = NULL)
#' # plot the first two dimensions
#' plotExpSpace(expSpace)
#' # plot another slice
#' plotExpSpace(expSpace, y = "P_ann_tot_m", x = "P_Feb_tot_m")
plotExpSpace <- function(expSpace,
                         y = expSpace[["attPerturb"]][1],
                         x = expSpace[["attPerturb"]][2]
                         ){
  
  targetMat <- expSpace[["targetMat"]]
  
  # Find column numbers of x & y
  targetAtts <- colnames(expSpace[["targetMat"]])
  x.no <- which(targetAtts == x)
  y.no <- which(targetAtts == y)
  
  xUnits <- getVarUnits(strsplit(x, "_")[[1]][1])
  yUnits <- getVarUnits(strsplit(y, "_")[[1]][1])
  xyFullNames <- mapply(tagBlender_noUnits, c(x,y))
  y.lab <- paste0(xyFullNames[2], " (", yUnits, ")")
  x.lab <- paste0(xyFullNames[1], " (", xUnits, ")")
  
  out <- expSpace2dViz(targetMat[ ,x.no], targetMat[ ,y.no], x.lab = x.lab, y.lab = y.lab)
  
}