#' @title Plot cell-specific decay curve for carcass persistence
#'
#' @description Produce the figure panel for a specific cell (factor
#'   level combination) including the specific fitted decay curves.
#'
#' @param model model of class cpm
#'
#' @param specificCell name of the specific cell to plot
#'
#' @param col color to use
#'
#' @param axis_x logical of whether or not to plot the x axis
#'
#' @param axis_y logical of whether or not to plot the y axis
#'
#' @export
#'
cpmCPCellPlot <- function(model, specificCell, col, axis_y = TRUE,
                          axis_x = TRUE){
# Juniper's original function
  dist <- model$dist
  CL <- model$CL
  cellwise <- model$cell_ab
  cellNames <- model$cells[ , "CellNames"]

  whichCarcs <- which(model$carcCell == specificCell)
  observations <- model$observations[whichCarcs, ]
  ncarc <- nrow(observations)

  whichSpecificCell <- which(cellNames == specificCell)
  a <- cellwise[whichSpecificCell, "pda_median"]
  b <- cellwise[whichSpecificCell, "pdb_median"]
  abs <- rcp(n = 1000, model = model, type = "ppersist")
  as <- abs[[whichSpecificCell]][ , "pda"]
  bs <- abs[[whichSpecificCell]][ , "pdb"]

  xVals <- model$observations[model$observations!= Inf]
  max_x <- max(xVals, na.rm = TRUE)
  max_x <- ceiling(max_x / 10) * 10
  pred_x <- seq(0, max_x, length.out = max_x * 10)
  pts <- pred_x[2:length(pred_x)]
  pta0 <- rep(0, length(pts))
  pta1 <- rep(0.000001, length(pts))
  CP <- ppersist(pda = as, pdb = bs, dist = dist,
    t_arrive0 = pta0, t_arrive1 = pta1, t_search = pts)
  pred_y <- apply(CP, 1, median)
  pred_yl <- apply(CP, 1, quantile, probs = (1 - CL) / 2)
  pred_yu <- apply(CP, 1, quantile, probs = 1 - (1 - CL) / 2)

  t1 <- observations[ , 1]
  t2 <- observations[ , 2]
  event <- rep(3, length(t1))
  event[which(!is.finite(t2))] <- 0
  event[which(t1 == t2)] <- 1
  t1[which(t1 == 0)] <- 0.0001
  survobj <- survival::Surv(t1, t2, event, "interval")
  form <- formula("survobj ~ 1")
  smod <- tryCatch(
    survival::survfit(form, data = observations),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  plot(smod, ylim = c(0, 1), xlim = c(0, max_x), # this gives the persistence probabilities
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o", lwd = c(2, 1, 1)
  )
  axis(1, las = 1, cex.axis = 0.9, at = seq(0, max_x, by = 10), labels = axis_x)
  axis(2, las = 1, cex.axis = 0.9, at = seq(0, 1, 0.2), labels = axis_y)
  points(pts, pred_y, type = "l", col = col, lwd = 3)
  points(pts, pred_yl, type = "l", col = col, lwd = 2, lty = 3)
  points(pts, pred_yu, type = "l", col = col, lwd = 2, lty = 3)
  countText <- paste("N = ", ncarc, sep = "")
  text(x = max_x, y = 0.95, countText, cex = 0.75, xpd = TRUE, adj = 1)
  text(max_x, 1.02, specificCell, adj = 1, cex = 0.75, font = 2)
}

#' @title Plot results of a single CP model
#'
#' @description Plot a single \code{\link{cpm}} model
#'
#' @param x model of class cpm
#'
#' @param col color to use
#'
#' @param ... to be passed down
#'
#' @examples
#'   data(wind_RP)
#'   mod <- cpm(formula_l = l ~ Season, formula_s = s ~ Season,  
#'            data = wind_RP$CP, left = "LastPresent", right = "FirstAbsent")
#'  plot(mod)
#'
#' @export
#'
plot.cpm <- function(x, col = "black", ...){
# loop through cells
# create header
  model <- x
  name_l <- format(model$formula_l)
  name_s <- format(model$formula_s)
  modelName <- paste(name_l, "; ", name_s, sep = "")
  if (model$dist == "exponential"){
    modelName <- name_l 
  }

  par(mar = c(0, 0, 0, 0), fig = c(0, 1, 0.95, 1))
  plot(1, 1, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", 
    ylab = "", ylim = c(0, 1), xlim = c(0, 1))

  points(c(0.01, 0.06), c(0.25, 0.25), type = 'l', lwd = 2, col = col)
  text(x = 0.07, y = 0.3, "= Median", adj = 0, cex = 0.9)
  points(c(0.2, 0.25), c(0.25, 0.25), type = 'l', lwd = 2, lty = 3, col = col)
  text(x = 0.26, y = 0.3, "= Confidence Bounds", adj = 0)

  labelsText <- paste(model$predictors, collapse = ".")
  text_label <- paste("Labels: ", labelsText, sep = "")
  text_model <- paste("Model: ", modelName, sep = "")
  text(x = 0.58, y = 0.3, text_label, adj = 0, cex = 0.75)
  text(x = 0.58, y = 0.7, text_model, adj = 0, cex = 0.75)
  par(fig = c(0, 1, 0, 1), mar = c(1, 1, 1, 1), new = TRUE)
  plot(1,1, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", 
    ylab = ""
  )
  mtext(side = 1, "Time", line = -0.25, cex = 1.5)
  mtext(side = 2, "Carcass Persistence", line = -0.25, cex = 1.5)

  ncell <- model$ncell
  cellNames <- model$cells[ , "CellNames"]

  nmatrix_col <- min(3, ncell)
  nmatrix_row <- ceiling(ncell / nmatrix_col)
  figxspace <- 0.925 / nmatrix_col
  figyspace <- 0.875 / nmatrix_row
  x1 <- rep(figxspace * ((1:nmatrix_col) - 1), nmatrix_row) + 0.05
  x2 <- rep(figxspace * ((1:nmatrix_col)), nmatrix_row) + 0.05
  y1 <- rep(figyspace * ((nmatrix_row:1) - 1), each = nmatrix_col) + 0.04
  y2 <- rep(figyspace * ((nmatrix_row:1)), each = nmatrix_col) + 0.04
  bottomCells <- seq(ncell - (nmatrix_col - 1), ncell, 1)
  leftCells <- which(1:ncell %% nmatrix_col == 1)
  if (length(leftCells) == 0){
    leftCells <- 1
  }

  for (celli in 1:ncell){
    par(mar = c(2.5, 2, 0, 0))
    par(fig = c(x1[celli], x2[celli], y1[celli], y2[celli]), new = T)
    specificCell <- cellNames[celli]
    axis_x <- FALSE
    axis_y <- FALSE
    if (celli %in% bottomCells){
      axis_x <- TRUE
    }
    if (celli %in% leftCells){
      axis_y <- TRUE
    }
    cpmCPCellPlot(model, specificCell, col, axis_y, axis_x)
  }
  par(.par_default)
}

#' @title Plot results of a set of CP models
#'
#' @description Produce a set of figures for a set of CP models, as fit by
#'   \code{\link{cpmSet}}
#'
#' @param x pk model set of class pkmSet
#'
#' @param specificModel the name(s) or index number(s) of specific model(s) to 
#'   restrict the plot
#'
#' @param cols named vector of the colors to use for the distributions
#'
#' @param ... to be passed down
#'
#' @examples
#'   data(wind_RP)
#'   mod <- cpmSet(formula_l = l ~ Season, formula_s = s ~ Season,  
#'            data = wind_RP$CP, left = "LastPresent", right = "FirstAbsent")
#'  \donttest{plot(mod)}
#'
#' @export
#'
plot.cpmSet <- function(x, specificModel = NULL, cols = NULL, ...){
  modelSet <- tidyModelSetCP(x)
  cells_set <- modelSetCells(modelSet)
  ncell <- nrow(cells_set)
  preds_set <- modelSetPredictors(modelSet)
  if (!is.null(specificModel) &&
     (!"character" %in% class(specificModel) || anyNA(specificModel))){
    stop("specific model must be specified by name (character string) or number")
  }
  if (is.null(specificModel)){
    specMods <- names(modelSet)
  } else {
    for (spi in specificModel){
      if (is.null(modelSet[[spi]])){ # spi not in modelSet (maybe exponential?)
        if (!grepl("exponential", spi))
          stop("specific CP model not in model set")
        tmp_spc <- strsplit(spi, ";")[[1]]
        tmp_spc[3] <- " NULL"
        if (is.null(modelSet[paste(tmp_spc, collapse = ";")]))
          stop("specific CP model not in model set")
      }
    }
    specMods <- specificModel
  }

  dist_set <- NULL # which distributions are included?
  # the possibilities are limited to names in the .cols_CP vector (internal data)
  for (di in names(.cols_CP))
    if(any(grepl(di, names(modelSet)))) dist_set <- c(dist_set, di)
  if (is.null(cols)){
    cols_CP <- .cols_CP[dist_set]
  } else if (length(dist_set) == 1 & length(cols) == 1){
    cols_CP <- cols
    names(cols_CP) <- dist_set
  } else {
    if (!all(dist_set %in% names(cols))){
      stop("'cols' must be NULL or a named vector of colors to be used in ",
        "plotting the distributions. Names are the names of the distributions ",
        "to associate with each respective color in the vector.")
    }
    cols_CP <- .cols_CP[dist_set]
  }
  preds_set <- modelSetPredictors(modelSet)
  n_col  <- ifelse(length(preds_set) == 0, 1,
    length(unique(cells_set[ , preds_set[1]])))

  n_row <- nrow(cells_set)/n_col
  for (spi in specMods){
    if (spi != specMods[1]){ # more than one specific & more than one fig to draw
      devAskNewPage(TRUE)
    }
    par(oma = c(5, 7.5, ifelse(nrow(cells_set) > 1, 9, 6), ifelse(n_row > 1, 4, 0.5)),
      mar = c(0, 0, 0, 0))
    par(mfrow = c(n_row, n_col))
    CPFig(modelSet, spi, cells_set, preds_set, cols_CP)
  }
  devAskNewPage(FALSE)
  par(.par_default)
}
# draw a CP figure graph
# NOTE: This is an internal utility function called by plot.cpmSet
#' @keywords internal
CPFig <- function(modelSet, specificModel, cells_set, preds_set, cols_CP){
  n_col  <- ifelse(length(preds_set) == 0, 1,
    length(unique(cells_set[ , preds_set[1]])))
  n_row <- nrow(cells_set)/n_col
  distr_featured <- sub(";", "", strsplit(specificModel, " ")[[1]][2])
  if (distr_featured == "exponential" && grepl("NULL", specificModel)){
    modelSet <- modelSet[grepl("exponential", names(modelSet))]
    cols_CP <- cols_CP["exponential"]
  }
  if (distr_featured == "exponential" && length(cols_CP) == 1){
     dists <- "exponential"
     forms <- strsplit(specificModel, "; ")[[1]]
     exind <- grep(paste0(forms[2], ";"), names(modelSet), fixed = TRUE)
     model_spc <- modelSet[[exind]]
     names(modelSet)[[exind]] <- specificModel
   } else {
    # reorder distributions to plot featured distr last
    dists <- c(names(cols_CP), distr_featured)
    dists <- dists[-min(which(dists == distr_featured))]
    if (("exponential" %in% dists) & length(dists) > 1){
      # change "NULL" to "s ~ <whatever's needed>"
      forms <- strsplit(specificModel, "; ")[[1]]
      exind <- intersect(
        grep("exponential", names(modelSet)),
        grep(paste0(forms[2], ";"), names(modelSet), fixed = TRUE))
      names(modelSet)[exind] <-
        sub("NULL", strsplit(specificModel, "; ")[[1]][3], names(modelSet)[exind])
    }
    if (distr_featured == "exponential" & length(dists > 1)){
      specificModel <- sub("exponential", dists[1], specificModel)
    }
    model_spc <- modelSet[[specificModel]]
   }
  # draw the cells
  for (celli in 1:nrow(cells_set)){
    specificCell <- cells_set[celli, "CellNames"]
    if (nrow(cells_set) == 1){
      observations <- model_spc$observations
    } else {
      carcCells <- apply(data.frame(model_spc$data[ , preds_set]),
        MARGIN = 1, FUN = paste, collapse = ".")
      observations <- model_spc$observations[carcCells == specificCell, ]
    }
    ncarc <- nrow(observations)
    cind <- which(cells_set$CellNames == specificCell)
    cell_spc <- matchCells(model_spc, modelSet)[cind]
    if ("exponential" %in% dists)
      cell_exp <- matchCells(modelSet[[exind]], modelSet)[cind]
    xVals <- observations[observations != Inf]
    max_x <- max(xVals, na.rm = TRUE)
    max_x <- ceiling(max_x / 10) * 10
    t1 <- observations[ , 1]
    t2 <- observations[ , 2]
    event <- rep(3, length(t1))
    event[which(!is.finite(t2))] <- 0
    event[which(t1 == t2)] <- 1
    t1[which(t1 == 0)] <- 0.0001

    smod <- tryCatch(
      survival::survfit(survival::Surv(t1, t2, event, type = "interval") ~ 1),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    plot(smod, ylim = c(0, 1), xlim = c(0, max_x), xlab = "", ylab = "",
      axes = FALSE, lwd = c(2, 1, 1))
    box()
    if (celli <= n_col & length(preds_set) > 0)
      mtext(side = 3, line = 0.3, text = cells_set[celli, preds_set[1]])
    if (celli %% n_col == 0 && length(preds_set) == 2)
      mtext(side = 4, line = 0.3, text = cells_set[celli, preds_set[2]])
    if (celli %% n_col == 1 | n_col == 1){
      axis(2, las = 1, at = seq(0, 1, by = 0.2))
      axis(2, las = 1, at = seq(0, 1, by = 0.1), labels = FALSE, tcl = -0.3)
    }
    if (celli > n_col * (n_row - 1))
      axis(1)
    xvals <- seq(0.001, max_x, length = 500)
    text(x = max_x, y = 0.95, paste0("N = ", ncarc), adj = 1)
    for (distr in dists){
      lwd  <- 2 + 2 * (distr == distr_featured)
      model_line <- gsub(model_spc$dist, distr, specificModel)
      mod <- modelSet[[model_line]]
      if (!is.null(mod)){
        if (distr == "exponential"){
          crow <- which(mod$cell_ab$cell == cell_exp)
          lines(xvals, 1 - pexp(xvals, rate = 1/mod$cell_ab[crow, "pdb_median"]),
            col = cols_CP[distr], lwd = lwd)
        } else {
          cind <- which(mod$cell_ab$cell == cell_spc)
          pda <- mod$cell_ab[cind, "pda_median"]
          pdb <- mod$cell_ab[cind, "pdb_median"]
          if (distr == "weibull"){
            lines(xvals, 1 - pweibull(xvals, shape = pda, scale = pdb),
              lwd = lwd, col = cols_CP[distr])
          } else if (distr == "loglogistic"){
            lines(xvals, 1 - pllogis(xvals, pda = pda, pdb = pdb),
              lwd = lwd, col = cols_CP[distr])
          } else if (distr == "lognormal"){
            lines(xvals, 1 - plnorm(xvals, meanlog = pdb, sdlog = sqrt(pda)),
              lwd = lwd, col = cols_CP[distr])
          }
        }
      }
    }
  }
  mtext(side = 1, line = 3.3, text = "Time", outer = TRUE, cex = 1.4)
  mtext(side = 2, line = 5.2, text = "Carcass Persistence", outer = TRUE, cex = 1.4)
  mtext(side = 2, line = 3.3, text = "Pr(persist > t)", cex = 1.2, outer = TRUE,)
  if (length(preds_set) > 0)
    mtext(side = 3, line = 2.5, text = preds_set[1], cex = 1.2, outer = TRUE)
  if (length(preds_set) == 2)
    mtext(side = 4, line = 2.5, text = preds_set[2], cex = 1.2, outer = TRUE)

  par(new = TRUE, mar = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mfrow = c(1, 1) )
  plot(0, type = "n", axes = FALSE)
  if (length(dists) > 1)
    dists <- dists[c(length(dists), 1:(length(dists) - 1))]
  lwd <- rep(2, length(dists))
  lwd[1] <- 4
  if (grepl("exponential", specificModel) & length(dists) == 1){
    # exponential model only is plotted
    legend("topleft",
      xpd = TRUE,
      legend = "exponential", lty = 1, lwd = lwd, ncol = 2, bty = "n",
      cex = 0.85, col = cols_CP["exponential"],
      title = paste("model:",
        gsub("~ 1", "constant", paste(strsplit(specificModel, ";")[[1]][2],
        collapse = ";"))),
      title.adj = 0)
  } else {
    legend("topleft", xpd = TRUE, legend = dists, lty = 1, col = cols_CP[dists],
      lwd = lwd, ncol = 2, bty = "n", cex = 0.85,
      title = paste("model:",
        gsub("~ 1", "~ constant",
          paste(strsplit(specificModel, ";")[[1]][-1], collapse = ";"))),
      title.adj = 0)
  }
}

#' @title Plot results of a single CP model in a set
#'
#' @description Produce a figures for a specific CP model, as fit by
#'   \code{\link{cpmSet}}
#'
#' @param modelSet cp model set of class \code{cpmSet}
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @param cols named vector of the colors to use for the distributions
#'
#' @return a plot
#'
#' @export
#'
plotCPFigure <- function(modelSet, specificModel, cols = CPcols()){
  plotCPHeader(modelSet, specificModel, cols)
  plotCPCells(modelSet, specificModel, cols)
}

#' @title Plot the cellwise results of a single model in a set of CP models
#'
#' @description Produce a set of cellwise figures for a specific CP model, as 
#'   fit by \code{\link{cpmSet}}
#'
#' @param modelSet cp model set of class cpmSet
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @param cols named vector of the colors to use for the distributions
#'
#' @return a plot
#'
#' @export
#'
plotCPCells <- function(modelSet, specificModel, cols){
## need cells_set, preds_set
  preds_set <- modelSetPredictors(modelSet)
  cells_set <- modelSetCells(modelSet)
  par(fig = c(0, 1, 0, 0.925), new = TRUE, mar = c(0, 0, 0, 0),
    oma = c(4, 4, 4, 0.5))
  plot(1, 1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "",
      ylab = "")

  mtext(side = 1, "Time", line = -0.25, cex = 1.5)
  mtext(side = 2, "Carcass Persistence", line = -0.25, cex = 1.5)

  n_col <- length(unique(cells_set[ , preds_set[1]]))
  n_row <- nrow(cells_set)/n_col
  H <- .res * (.header + n_row * .panel_H)
  W <- min(1200, .res * .panel_W * max(2, n_col))
  cell_W <- round(W/n_col)
  box_W <- round(W/2)

  # cells
  if (nrow(cells_set) > 1) omd = c(.res/W, 1 - .res/W,
      .res * .footer/H, 1 - .res * (.header + .box_H + 0.75 * .buffer)/H)
  if (nrow(cells_set) == 1) omd = c(.res/W, 0.96,
      .res * .footer/H, 1 - .res * (.header + .box_H + 0.2 * .buffer)/H)

  par(omd = omd, mar = c(0, 0, 0, 0), xaxt = "n", yaxt = "n")
  par(mfrow = c(n_row, n_col))

  cells_set <- modelSetCells(modelSet)
  ncell <- nrow(cells_set)
  cellNames <- cells_set[ , "CellNames"]
  nmatrix_col <- min(3, ncell)
  nmatrix_row <- ceiling(ncell / nmatrix_col)
  figxspace <- 0.925 / nmatrix_col
  figyspace <- 0.875 / nmatrix_row
  x1 <- rep(figxspace * ((1:nmatrix_col) - 1), nmatrix_row) + 0.04
  x2 <- rep(figxspace * ((1:nmatrix_col)), nmatrix_row) + 0.04
  y1 <- rep(figyspace * ((nmatrix_row:1) - 1), each = nmatrix_col) + 0.03
  y2 <- rep(figyspace * ((nmatrix_row:1)), each = nmatrix_col) + 0.03
  bottomCells <- seq(ncell - (nmatrix_col - 1), ncell, 1)
  leftCells <- which(1:ncell %% nmatrix_col == 1)
  if (length(leftCells) == 0){
    leftCells <- 1
  }
  for (celli in 1:ncell){
    par(mar = c(2.5, 2, 0, 0))
    par(fig = c(x1[celli], x2[celli], y1[celli], y2[celli]), new = TRUE)
    specificCell <- cellNames[celli]
    axis_x <- FALSE
    axis_y <- FALSE
    if (celli %in% bottomCells){
      axis_x <- TRUE
    }
    if (celli %in% leftCells){
      axis_y <- TRUE
    }
    axes <- c("x" = axis_x, "y" = axis_y)
    cpmSetSpecCPCellPlot(modelSet, specificModel, specificCell, cols, axes)
  }
}

#' @title Plot cell-specific decay curve for carcass persistence
#'
#' @description Produce the figure panel for a specific cell (factor 
#'   level combination) including the specific fitted decay curves.
#'
#' @param modelSet modelSet of class cpmSet
#'
#' @param specificModel name of the specific submodel to plot
#'
#' @param specificCell name of the specific cell to plot
#'
#' @param cols named vector of the colors to use for the distributions
#'
#' @param axes named vector of logical values indicating whether or not to 
#'   plot the x axis and the y axis
#'
#' @return a specific cell plot panel
#'
#' @export
#'
cpmSetSpecCPCellPlot <- function(modelSet, specificModel, specificCell,
                                 cols, axes){

  model_spec <- modelSet[[specificModel]]
  dist <- model_spec$dist
  CL <- model_spec$CL
  cellNames_set <- modelSetCells(modelSet)[ , "CellNames"]
  preds_set <- modelSetPredictors(modelSet)
  carcCells <- apply(data.frame(model_spec$data[ , preds_set]),
    MARGIN = 1, FUN = paste, collapse = ".")
  whichCarcs <- which(carcCells == specificCell)
  if (specificCell == "all"){
    whichCarcs <- 1:length(carcCells)
  }
  observations <- model_spec$observations[whichCarcs, ]
  ncarc <- nrow(observations)
  cellMatch_spec <- matchCells(model_spec, modelSet)
  whichSpecificCell_spec <- cellMatch_spec[cellNames_set == specificCell]

  xVals <- model_spec$observations[model_spec$observations!=Inf]
  max_x <- max(xVals, na.rm = TRUE)
  max_x <- ceiling(max_x / 10) * 10
  t1 <- observations[ , 1]
  t2 <- observations[ , 2]
  event <- rep(3, length(t1))
  event[which(!is.finite(t2))] <- 0
  event[which(t1 == t2)] <- 1
  t1[which(t1 == 0)] <- 0.0001

  smod <- tryCatch(
    survival::survfit(survival::Surv(t1, t2, event, type = "interval") ~ 1),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  plot(smod, ylim = c(0, 1), xlim = c(0, max_x),
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o", lwd = c(2, 1, 1)
  )
  axis(1, las = 1, cex.axis = 0.9, at = seq(0, max_x, by = 10), 
    labels = axes["x"]
  )
  axis(2, las = 1, cex.axis = 0.9, at = seq(0, 1, 0.2), labels = axes["y"])
  countText <- paste("N = ", ncarc, sep = "")
  text(x = max_x, y = 0.95, countText, cex = 0.75, xpd = TRUE, adj = 1)
  text(max_x, 1.02, specificCell, adj = 1, cex = 0.75, font = 2)

  nmodelsInSet <- length(modelSet)
  modelSetNames <- names(modelSet)
  modelSetNames_components <- unlist(strsplit(modelSetNames, "; "))
  modelSetNames_matrix <- matrix(modelSetNames_components, ncol = 3, byrow = TRUE)

  colnames(modelSetNames_matrix) <- c("dist", "form_l", "form_s")

  modelMatches <- vector("list", length = nmodelsInSet)
  for (modelsInSeti in 1:nmodelsInSet){
    mod_l <- gsub(" ", "", modelSetNames_matrix[modelsInSeti, 2])
    mod_s <- gsub(" ", "", modelSetNames_matrix[modelsInSeti, 3])
    modComp_l <- gsub(" ", "", modelSetNames_matrix[ , 2])
    modComp_s <- gsub(" ", "", modelSetNames_matrix[ , 3])
    matches <- which(modComp_l == mod_l & modComp_s == mod_s)
    modelMatches[[modelsInSeti]] <- matches
  } 

  whichSpecificModel <- which(names(modelSet) == specificModel)
  modsToPlot <- modelMatches[[whichSpecificModel]]
  modsToPlot <- modsToPlot[-which(modsToPlot == whichSpecificModel)]
  modsToPlot <- c(modsToPlot, whichSpecificModel)
  nmodsToPlot <- length(modsToPlot)

  for (modsToPloti in 1:nmodsToPlot){
    model_i <- modelSet[[modsToPlot[modsToPloti]]]
    dist_i <- model_i$dist
    cellwise_i <- model_i$cell_ab
    cellNames_i <- cellwise_i[ , "cell"]
    cellMatch_i <- matchCells(model_i, modelSet)
    reducedCell_i <- cellMatch_i[cellNames_set == specificCell]
    whichSpecificCell <- which(cellNames_i == reducedCell_i)
    a <- cellwise_i[whichSpecificCell, "pda_median"]
    b <- cellwise_i[whichSpecificCell, "pdb_median"]
    pred_x <- seq(0, max_x, length.out = max_x * 10)
    pts <- pred_x[2:length(pred_x)]
    pta0 <- rep(0, length(pts)) 
    pta1 <- rep(0.000001, length(pts))
    pred_y <- t(ppersist(pda = a, pdb = b, dist = dist_i,
      t_arrive0 = pta0, t_arrive1 = pta1, t_search = pts))

    col_i <- cols[dist_i]
    lwd <- 2
    if (modsToPloti == nmodsToPlot){
      lwd <- 5
    }
    points(pts, pred_y, type = "l", col = col_i, lwd = lwd)
  }

}

#' @title The CP plot header
#'
#' @description Produce the header for a CP plot
#'
#' @param modelSet cp model set of class cpmSet
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @param cols named vector of the colors to use for the distributions
#'
#' @return a plot
#'
#' @export
#'
plotCPHeader <- function(modelSet, specificModel, cols = CPcols()){

  modelSetNames <- names(modelSet)
  nmodelsInSet <- length(modelSetNames)
  whichSpecificModel <- which(names(modelSet) == specificModel)

  modelSetNames_components <- unlist(strsplit(modelSetNames, "; "))
  modelSetNames_matrix <- matrix(modelSetNames_components, ncol = 3,
    byrow = TRUE)
  colnames(modelSetNames_matrix) <- c("dist", "form_l", "form_s")
  distsIncluded <- gsub("dist: ", "", unique(modelSetNames_matrix[ , "dist"]))

  par(fig = c(0, 1, 0.925, 1), mar = c(0, 0, 0, 0))
  plot(1,1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", 
    ylab = "", ylim = c(0, 1), xlim = c(0, 1))

  x1 <- c(0, 0.29, 0, 0.29)
  y1 <- c(0.4, 0.4, 0.0, 0.0)
  x2 <- c(0.04, 0.33, 0.04, 0.33)
  y2 <- c(0.6, 0.6, 0.2, 0.2)
  xt <- c(0.05, 0.34, 0.05, 0.34)
  yt <- c(0.5, 0.5, 0.1, 0.1)
  ndists <- length(distsIncluded)
 
  for (disti in 1:ndists){
    rect(x1[disti], y1[disti], x2[disti], y2[disti], border = NA,
      col = cols[distsIncluded[disti]]
    )
    distName <- distsIncluded[disti]
    distName <- paste(toupper(substring(distName, 1, 1)),
                      substring(distName, 2), sep = "", collapse = " ")
    distText <- paste0("= ", gsub("Log", "Log-", distName))
    text(x = xt[disti], y = yt[disti], adj = 0, cex = 0.8, distText)
  }

  labelsText <- paste(modelSetPredictors(modelSet), collapse = ".")
  text_label <- paste("Labels: ", labelsText, sep = "")
  forms <- modelSetNames_matrix[whichSpecificModel, c("form_l", "form_s")]
  modelsText <- paste(forms, collapse = "; ")
  modelsText <- gsub("~ 1", "~ constant", modelsText)
  text_model <- paste("Model: ", modelsText, sep = "")
  text(x = 0.59, y = 0.1, text_label, adj = 0, cex = 0.75)
  text(x = 0.59, y = 0.5, text_model, adj = 0, cex = 0.75)
}

#' @title Error check a specific model selection for a CP plot
#'
#' @description Make sure it's available and good, update the name for usage
#'
#' @param modelSet cp model set of class cpmSet
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @return updated name of the model to use
#'
#' @export
#'
checkSpecificModelCP <- function(modelSet, specificModel){
  if (!is.null(specificModel) && anyNA(specificModel)){
    stop(
      "specificModel must be NULL or a vector of model names or positions.",
      "\nNAs not allowed."
    )
  }
  if (length(specificModel) > 0){
    if (is.numeric(specificModel)){
      if (anyNA(specificModel)){
        warning("specificModel cannot be NA. NA models removed.")
        specificModel <- specificModel[!is.na(specificModel)]
        if (length(specificModel) == 0){
          stop("No valid specificModel")
        }
      }
      if (any(specificModel > length(modelSet))){
        stop(paste0("there are only ", length(modelSet), " model choices."))
      }
      specificModel <- names(modelSet)[specificModel]
    } else{
      specificModel <- gsub("NULL", "s ~ 1", specificModel)
    }
    if (any(specificModel %in% names(modelSet)) == FALSE){
      stop("Selected model not in set. To see options use names(modelSet).")
    }
    modNames <- specificModel
    for (modi in modNames){
      if (cpmFail(modelSet[[modi]])){
        stop("specificModel ", modi, " is not a well-fit model")
      }
    }
  } else{
    specificModel <- names(cpmSetFailRemove(modelSet))
  }
  return(specificModel)
}

#' @title Expand a CP model set for plotting
#'
#' @description Expand the exponential models across the other distributions
#'
#' @param modelSet cp model set of class cpmSet
#'
#' @return updated model set
#'
#' @export
#'
expandModelSetCP <- function(modelSet){

  modelSetNames <- names(modelSet)
  nmodelsInSet <- length(modelSetNames)

  modelSetNames_components <- unlist(strsplit(modelSetNames, "; "))
  modelSetNames_matrix <- matrix(modelSetNames_components, ncol = 3, 
                            byrow = TRUE
                          )
  colnames(modelSetNames_matrix) <- c("dist", "form_l", "form_s")
  whichExp <- which(modelSetNames_matrix[ , "dist"] == "dist: exponential")
  modelSetNames_matrix[whichExp, "form_s"] <- "s ~ 1"
  if (length(whichExp) > 0){
    uniqueForm_s <- unique(modelSetNames_matrix[ , "form_s"])
    names(modelSet)[whichExp] <- gsub("NULL", "s ~ 1", 
                                   names(modelSet)[whichExp]
                                 )
    if (length(uniqueForm_s) > 1){
      uniqueForm_s <- uniqueForm_s[-which(uniqueForm_s == "s ~ 1")]
      nuniqueForm_s <- length(uniqueForm_s)
      for (uniqueForm_si in 1:nuniqueForm_s){
        replacement <- uniqueForm_s[uniqueForm_si]
        newEntries <- modelSetNames_matrix[whichExp, ]
        newEntries[ , "form_s"] <- replacement
        modelSetNames_matrix <- rbind(modelSetNames_matrix, newEntries)
        newEntries <- modelSet[whichExp]
        names(newEntries) <- gsub("s ~ 1", replacement, names(newEntries))
        modelSet <- c(modelSet, newEntries)
      }
    }
  }
  return(modelSet)
}

#' @title Tidy a CP model set
#'
#' @description Remove bad fit models
#'
#' @param modelSet cp model set of class cpmSet
#'
#' @return a trimmed model set
#'
#' @export
#'
tidyModelSetCP <- function(modelSet){
  modelSet <- cpmSetFailRemove(modelSet)
  modelSet <- modelSet[order(sapply(modelSet, "[[", "AICc"))]
  class(modelSet) <- c("cpmSet", "list")
  return(modelSet)
}

#' @title Produce a named vector of standard CP plot colors
#' 
#' @description Produce a named vector of standard CP plot colors
#'
#' @export
#'
CPcols <- function(){
c(exponential = rgb(0.80, 0.38, 0.56), 
  weibull = rgb(1.00, 0.76, 0.15),
  loglogistic = rgb(0.00, 1.00, 1.00),
  lognormal = rgb(0.00, 0.41, 0.55))
}


