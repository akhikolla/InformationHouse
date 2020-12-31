#' @title Plot summary statistics for splits of mortality estimates
#'
#' @description The S3 \code{plot} method for \code{splitSummary} objects 
#'   constructs boxplots of the mortality estimates for all combinations of 
#'   splitting covariates summarized in the \code{splits} variable.
#'
#' For 1-covariate splits, box plots showing median, IQR, and confidence
#'   intervals (for the \code{CL} attribute for the splits object). For
#'   2-covariate splits, the box plots are in an array with levels of
#'   the temporal split (\code{split_SS} or \code{split_time}) arranged
#'   horizontally (if present) and the levels of the \code{split_CO} variable
#'   arranged vertically. If no temporal splits are present, then the box 
#'   plots along the levels of the first \code{split_CO} variable are arranged
#'   horizontally and the levels of the second variable are are arranged
#'   vertically.
#'
#' @param x A \code{splitSummary} object (result of \code{\link{calcSplits}})
#'  that includes summary statistics for simulated mortality estimates for all
#'  combinations of levels of 1 or 2 splitting covariates.
#'
#' @param rate \code{logical} scalar indicating whether the figures should be
#'  plotted as number of fatalities per split category (\code{rate = TRUE}) or
#'  fatality rates per unit time (\code{rate = TRUE}). If the splits do not
#'  include either a \code{split_SS} or \code{split_time} variable, the
#'  \code{rate} arg is ignored.
#'
#' @param commonScale boolean to indicate whether to plot separate splits panels
#'  with a common scale on their y-axes (or have y-axes scaled to fit each graph
#'  separately)
#'
#' @param ... additional arguments to be passed down
#'
#' @export
#'
plot.splitSummary <- function(x, rate = FALSE, commonScale = FALSE, ...){
  splits <- x
#  nvar <- ifelse(is.na(attr(splits, "vars")), 0, attr(splits, "vars"))
  nvar <- length(attr(splits, "vars"))
  vartype <- attr(splits, "type")
  if (nvar == 0 || vartype[1] == "CO") rate <- FALSE
  vars <- attr(splits, "vars")
  alpha <- 1 - attr(splits, "CL")
  times <- attr(splits, "times")
  probs <- c(alpha/2, 0.25, 0.5, 0.75, 1 - alpha/2)
  deltaT <- diff(times)
  if (nvar == 0){
    # split is for total...what kind of figure? Need to show Mhat
    simpleMplot(splits, ...)
    return(NULL)
  } else if (nvar == 1 & !is.list(splits)){
    splits <- list(splits)
    vnames <- NULL
  } else {
    vnames <- names(splits)
  }
  hnames <- rownames(splits[[1]])
  par(.par_splits)
  nlevel_h <- nrow(splits[[1]])
  nlevel_v <- length(splits)
  par(mfrow = c(nlevel_v, 1))
  cex.axis <- 1*(nlevel_v == 1) + (nlevel_v == 2)/0.83 + (nlevel_v > 2)/0.66
  if (!rate & commonScale & nvar == 2){
    ylim  <- c(0, 0)
    for (vi in 1:nlevel_v){
      ylim <- range(ylim, splits[[vi]][ , -1])
    }
  }
  for (vi in 1:nlevel_v){
    if ((vartype[1] %in% c("time", "SS")) & rate) {
      hwid <- deltaT/2
      xlim <- range(times)
      ylim <- range(splits[[vi]][ , -1]/deltaT)
    } else {
      hwid <- rep(0.15, nlevel_h) # half-width of boxes
      xlim <- c(1, nlevel_h) + 0.5 * c(-1, 1)
      if (!commonScale) ylim <- range(splits[[vi]][ , -1])
    }
    if (vi == 1 && !is.null(try(plot.new(), silent = TRUE))){
      par(mfrow = c(1, 1))
      return(1) # flags the plot as as an error
    } else {
      par(new = (vi == 1))
      plot(0, xlim = xlim, ylim = ylim, type = "n", axes = F, xlab = "", 
        ylab = "")
    }
    if (vartype[1] == "CO" | !rate){
      xx <- 1:nlevel_h 
    } else {
      xx <- times[-1] - hwid
    }
    for (hi in 1:nlevel_h){
      ratebars <- !(vartype[1] == "CO" || !rate)
      deno <- ifelse(ratebars, deltaT[hi], 1)
      qtls <- splits[[vi]][hi, -1]/deno
      polygon(xx[hi] + hwid[hi] * c(1, 1, -1, -1), qtls[c(2, 4, 4, 2)])
      lines(xx[hi] + hwid[hi] * c(1, -1), rep(qtls[3], 2), lwd = 3)
      if (alpha >= 0.5) yst <- c(3, 3) else yst <- c(2, 4)
      lines(rep(xx[hi], 2), qtls[c(1, yst[1])])
      lines(rep(xx[hi], 2), qtls[c(yst[2], 5)])
      lines(xx[hi] + hwid[hi]/(2 - ratebars) * c(1, -1), rep(qtls[1], 2))
      lines(xx[hi] + hwid[hi]/(2 - ratebars) * c(1, -1), rep(qtls[5], 2))
    }
    axis(2, las = 1)
    
    ymid <- mean(par("usr")[3:4])
    if (!rate | vartype[1] == "CO"){
      if (vartype[1] == "time"){
        at <- 0.5 + 0:nlevel_h
        lab <- c(0, rownames(splits[[vi]]))
      } else {
        at <- 1:nlevel_h
        lab <- rownames(splits[[vi]])
      }
    } else {
      at <- times
      lab <- times
    }
    if (vi == nlevel_v){
      axis(1, at = at, labels = lab, cex.axis = cex.axis)
    }
    if (nlevel_v > 1){
      axis(4, at = mean(par("usr")[c(3, 4)]), labels = vnames[vi],
        tck = 0, mgp = c(3, 0.5, 0), cex.axis = cex.axis)
    }
    box()
  }
  mtext(side = 1, vars[1], line = 4.6, cex = 1.2)
  if (vartype[1] == "SS" & rate) {
    mtext(side = 1, line = 3, text = hnames, at = times[-1] - diff(times)/2)
  }
  ylab <- ifelse(vartype[1] == "CO" | !rate,
    paste0("Estimated M"),
    paste0("Estimated M/day")
  )
  mtext(side = 2, text = ylab, outer = T, line = 3.5, cex = 1.2)

  title <- ifelse(vartype[1] == "CO" | !rate,
    paste0("Estimated mortality by ", vars[1]),
    paste0("Estimated daily mortality rate")
  )
  mtext(text = title, side = 3, line = 2, cex = 1.3, outer = T)
  mtext(side = 3, line = 0.5, cex = 0.9, outer = T,
    text = paste0("Median, IQR, and ", 100*(1 - alpha), 
             "% confidence intervals")
           )
  mtext(side = 4, text = vars[2], outer = T, line = 2.5, cex = 1.2)
  par(.par_default)
  return(NULL)
}

#' @title Plot summary statistics for splits of mortality estimates
#'
#' @description The S3 \code{plot} method for \code{splitFull} objects 
#'   constructs boxplots of the mortality estimates for all combinations of 
#'   splitting covariates summarized in the \code{splits} variable. This is a
#'   simple wrapper function for creating a \code{splitSummary} object by 
#'   calling \code{\link{summary.splitFull}} and plotting the result via
#'   \code{\link{plot.splitSummary}}.
#'
#' @param x A \code{splitSummary} object (result of \code{\link{calcSplits}})
#'   that includes summary statistics for simulated mortality estimates for
#'   all combinations of levels of 1 or 2 splitting covariates.
#'
#' @param rate \code{logical} scalar indicating whether the figures should be
#'   plotted as number of fatalities per split category (\code{rate = TRUE})
#'   or fatality rates per unit time (\code{rate = TRUE}). If the splits do
#'   not include either a \code{split_SS} or \code{split_time} variable, the
#'   \code{rate} arg is ignored.
#'
#' @param CL desired confidence level to show in box plots
#'
#' @param commonScale Boolean: Should panels share a common y-axis scale?
#'  Relevant only when there are two splitting variables.
#'
#' @param ... to be passed down
#'
#' @export
#'
plot.splitFull <- function(x, rate = FALSE, CL = 0.90, commonScale = FALSE,...){
  nvar <- length(attr(x, "vars"))
  if (nvar == 0){
    simpleMplot(x, ..., CL = CL)
  } else {
    splitSum <- summary(x, CL)
    if(!is.null(plot(splitSum, rate, CL = CL, commonScale = commonScale))){
      stop("Second split too fine for plotting. Consider transposing.")
    }
  }
}
