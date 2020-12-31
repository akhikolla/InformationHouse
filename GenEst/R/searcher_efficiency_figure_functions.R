#' @title Plot results of a single pk model
#'
#' @description Plot a single \code{\link{pkm}} model
#'
#' @param x model of class pkm
#'
#' @param col color to use
#'
#' @param CL confidence level to show in boxplots and confidence bounds
#'
#' @param ... arguments to be passed to sub functions
#'
#' @return a plot
#'
#' @examples
#'   data(wind_RP)
#'   mod <- pkm(formula_p = p ~ Season, formula_k = k ~ 1, data = wind_RP$SE)
#'   plot(mod)
#'
#' @export
#'
plot.pkm <- function(x, col = NULL, CL = NULL, ...){
  model <- x
  if (anyNA(model$varbeta) || sum(diag(model$varbeta) < 0) > 0){
    stop("Variance in pkm not well-defined. Cannot plot.")
  }
  name_p <- format(model$formula_p)
  name_k <- model$formula_k
  if (!is.null(model$pOnly) && model$pOnly){
    stop("k missing from pk model. Cannot plot.")
  }
  if (class(name_k) == "numeric"){
    name_k <- paste("k fixed at ", name_k, sep = "")
  } else if (class(name_k) == "character"){
    name_k <- "k not estimated"
   }else {
    name_k <- format(model$formula_k)
  }
  modelName <- paste0(name_p, "; ", name_k)

  if (is.null(CL)) CL <- model$CL
  if (is.null(col) || anyNA(col)){
    cols_SE <- .cols_SE
  }
  cells <- model$cells
  ncell <- nrow(cells)
  preds <- model$predictors
  ### summarize empical trial carcass observation data for drawing cell panel plots
  obsCol <- colnames(model$observations)
  nsearch <- length(obsCol)

  ### data structures for points in the SE decay panels
  # "found" and "available" give number of carcasses found and available on each
  # search (columns) for each cell (rows)
  found <- available <- matrix(nrow = ncell, ncol = nsearch,
    dimnames = list(cells$CellNames, obsCol))

  obs <- as.matrix(model$data0[ , obsCol])
  rownames(obs) <- switch(paste(length(preds)),
    "0" = rep("all", nrow(obs)),
    "1" = model$data0[ , preds],
    "2" = do.call(paste, c(model$data0[ , preds] , sep = "."))
  )
  for (cname in cells$CellNames){
    oind <- which(rownames(obs) == cname)
    found[cname, ] <- matrixStats::colCounts(
      matrix(obs[oind, ], ncol = ncol(obs)), value = 1, na.rm = TRUE)
    available[cname, ] <- matrixStats::colCounts(
      !is.na(matrix(obs[oind, ], ncol = ncol(obs))), value = TRUE)
  }

  # quantiles of p and k by cell (for plotting boxplots and fitted lines)
  lwr <- (1 - CL)/2
  upr <- 1 - (1 - CL)/2
  p <- c(bot = 0.005, lwr = lwr, qr1 = 0.25, med = 0.50,
          qr3 = 0.75, upr = upr, top = 0.995)

  pk <- qpk(p, model)
  gdat <- cbind(pk$p, pk$k)
  colnames(gdat) <-  c(paste0("p", names(p)), paste0("k", names(p)))
  specificModel <- NULL

  if (!is.null(model$pOnly) && model$pOnly){
    plot(0, 0, type = 'n', axes = F, xlab = '', ylab = '')
    text(0, .5, "k missing from pk model. Cannot plot.", cex = 2, col = 2)
  } else {
    SEfig(referenceModel = NULL, specificModel = modelName,
        gdat_spc = gdat, gdat_ref = NULL,
        found = found, available = available,
        cells_set = cells, cols_SE = cols_SE, CL = CL)
  }
  par(.par_default)
}

#' @title Plot results of a set of SE models
#'
#' @description Produce a set of figures for a set of SE models, as fit by
#'   \code{\link{pkmSet}}
#'
#' @param x pk model set of class pkmSet
#'
#' @param specificModel the name(s) or index number(s) of specific model(s)
#'   to restrict the plot
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @param CL confidence level
#'
#' @param ... to be sent to subfunctions
#'
#' @return a set of plots
#'
#' @examples
#'   data(wind_RP)
#'   mod <- pkmSet(formula_p = p ~ Season, formula_k = k ~ Season,
#'            data = wind_RP$SE
#'          )
#'   plot(mod)
#'
#' @export
#'
plot.pkmSet <- function(x, specificModel = NULL, cols = NULL, CL = NULL, ...){
  # data extraction and formatting
  modelSet <- x
  specMods <- checkSpecificModelSE(modelSet, specificModel)
  modelSet <- tidyModelSetSE(modelSet)
  if (is.null(CL)) CL <- modelSet[[1]]$CL
  if (is.null(cols) || anyNA(cols)){
    cols_SE <- .cols_SE
  } else if (length(cols) < 2 | !is.vector(cols)){
    stop ("cols must be a vector of length 2")
  } else {
    cols_SE <- cols[1:2]
    names(cols_SE) <- c("spc", "ref")
  }
  cells_set <- modelSetCells(modelSet)
  ncell <- nrow(cells_set)
  model_ref <- refMod(modelSet)
  referenceModel <- paste0(
    deparse(model_ref$formula_p), "; ",
    ifelse("formula" %in% class(model_ref$formula_k),
      deparse(model_ref$formula_k),
      paste("k fixed at ", model_ref$formula_k)
    )
  )
  preds_set <- modelSetPredictors(modelSet)
  preds_ref <- model_ref$predictors
  if (length(preds_ref) > 0 & !anyNA(preds_ref))
    key_ref <- match(preds_ref, preds_set)
  ### summarize empical trial carcass observation data for drawing cell panel plots
  obsCol <- colnames(model_ref$observations)
  nsearch <- length(obsCol)

  ### data structures for points in the SE decay panels
  # "found" and "available" give number of carcasses found and available on each
  # search (columns) for each cell (rows)
  found <- available <- matrix(nrow = ncell, ncol = nsearch,
    dimnames = list(cells_set$CellNames, obsCol))

  obs <- as.matrix(model_ref$data0[ , obsCol])
  rownames(obs) <- switch(paste(length(preds_set)),
    "0" = rep("all", nrow(obs)),
    "1" = model_ref$data0[ , preds_ref],
    "2" = do.call(paste, c(model_ref$data0[ , preds_ref[key_ref]] , sep = "."))
  )
  for (cname in cells_set$CellNames){
    oind <- which(rownames(obs) == cname)
    found[cname, ] <- matrixStats::colCounts(
      matrix(obs[oind, ], ncol = ncol(obs)), value = 1, na.rm = TRUE)
    available[cname, ] <- matrixStats::colCounts(
      !is.na(matrix(obs[oind, ], ncol = ncol(obs))), value = TRUE)
  }
  # quantiles of p and k by cell (for plotting boxplots and fitted lines)
  lwr <- (1 - CL)/2
  upr <- 1 - (1 - CL)/2
  p <- c(bot = 0.005, lwr = lwr, qr1 = 0.25, med = 0.50,
          qr3 = 0.75, upr = upr, top = 0.995)

  pk_ref <- qpk(p, model_ref)
  gdat_ref <- gdat_spc <- matrix(nrow = ncell, ncol = 2*length(p),
    dimnames = list(
      cells_set$CellNames,
      c(paste0("p", names(p)), paste0("k", names(p)))))
  attr(gdat_ref, "p") <- p
  gdat_spc <- gdat_ref
  if (length(preds_ref) == 0 || anyNA(preds_ref)){
    for (ci in rownames(gdat_ref)){
      gdat_ref[ci, paste0("p", names(p))] <- pk_ref$p
      gdat_ref[ci, paste0("k", names(p))] <- pk_ref$k
    }
    key_ref <- NULL
  }
  nmod <- length(specMods)
  for (modi in 1:nmod){
    model_spc <- modelSet[[specMods[modi]]]
    preds_spc <- model_spc$predictors
    pk_spc <- qpk(p, model_spc)
    specificModel <- specMods[modi]
    if (length(preds_spc) == 0 || anyNA(preds_spc)){
      for (ci in rownames(gdat_spc)){
        gdat_spc[ci, paste0("p", names(p))] <- pk_spc$p
        gdat_spc[ci, paste0("k", names(p))] <- pk_spc$k
      }
      key_spc <- NULL
    } else {
      key_spc <- match(preds_spc, preds_set)
    }
    for (celli in 1:ncell){
      if (length(preds_set) == 0) break # NOTE: not plotted! ...needs correcting
      pi1 <- cells_set[celli, preds_set[1]]
      if (length(preds_set) == 2){
        pi2 <- cells_set[celli, preds_set[2]]
        ci <- paste(pi1, pi2, sep = ".")
      } else if (length(preds_set) == 1){
        pi2 <- NULL
        ci <- pi1
      }
      if (length(preds_spc) > 0){
         cellname_spc <- paste(unlist(mget(paste0("pi", key_spc))), collapse = ".")
         gdat_spc[ci, ] <- c(pk_spc$p[cellname_spc, ], pk_spc$k[cellname_spc, ])
      }
      if (length(preds_ref) > 0){
        cellname_ref <- paste(unlist(mget(paste0("pi", key_ref))), collapse = ".")
        gdat_ref[ci, ] <- c(pk_ref$p[cellname_ref, ], pk_ref$k[cellname_ref, ])
      }
    }
    if (modi == 2){
      devAskNewPage(TRUE)
    }
    if (!is.null(modelSet[[modi]]$pOnly) && modelSet[[modi]]$pOnly){
      plot(0, 0, type = 'n', axes = F, xlab = '', ylab = '')
      text(0, .5, "k missing from pk model. Cannot plot.", cex = 2, col = 2)
    } else {
      SEfig(referenceModel = referenceModel, specificModel = specificModel,
          gdat_spc = gdat_spc, gdat_ref = gdat_ref,
          found = found, available = available,
          cells_set = cells_set, cols_SE = cols_SE, p = p, CL = CL)
    }
  }
  devAskNewPage(FALSE)
  par(.par_default)
}


#' @title Plot results of a single SE model in a set
#'
#' @description Produce a figures for a specific SE model, as fit by
#'   \code{\link{pkmSet}}
#'
#' @keywords internal
SEfig <- function(referenceModel, specificModel, gdat_spc, gdat_ref,
  found, available, cells_set, cols_SE, p, CL){
  # data prep
  preds_set <- names(cells_set)[-length(names(cells_set))]
  nsearch <- ncol(found)
# cells
  xends <- nsearch/8 # how much wider than range(x) should cell axes be?
  # plot configuration
  n_col <- length(unique(cells_set[ , preds_set[1]]))
  n_row <- nrow(cells_set)/n_col
  H <- .res * (.header + .box_H + .buffer +  n_row * .panel_H + .footer)
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

  x <- 1:nsearch # x-values (shared among all cell panels)
  if (nrow(cells_set) == 1){
    y_spc <- gdat_spc["all", "pmed"] * gdat_spc["all", "kmed"] ^ (0:(nsearch - 1))
        SEpanel(found = found, available = available,
      y_spc = y_spc, y_ref = y_spc, xends = xends, cols_SE)
    axis(1, at = 1:ncol(available), xpd = TRUE, xaxt = "s")
    axis(2, xpd = TRUE, yaxt = "s", at = 0.1 * c(0, 2, 4, 6, 8, 10), las = 1)
    axis(2, xpd = TRUE, yaxt = "s", at = 0.1 * 0:10,
      las = 1, labels = FALSE, tcl = -0.25)

  } else if (length(preds_set) == 1){
    lev1 <- unique(cells_set[ , preds_set])
    for (pr1 in lev1){
      y_spc <- gdat_spc[pr1, "pmed"] * gdat_spc[pr1, "kmed"] ^ (0:(nsearch - 1))
      if (!is.null(gdat_ref)){
        y_ref <- gdat_ref[pr1, "pmed"] * gdat_ref[pr1, "kmed"] ^ (0:(nsearch - 1))
      } else {
        y_ref <- NULL
      }
      SEpanel(found = found[pr1, ], available = available[pr1, ],
        y_spc = y_spc, y_ref = y_ref, xends = xends, cols_SE)
      axis(1, at = 1:ncol(available), xpd = TRUE, xaxt = "s")
      mtext(pr1, side = 3, xpd = TRUE)
      if (pr1 == lev1[1]){
        axis(2, xpd = TRUE, yaxt = "s", at = 0.1 * c(0, 2, 4, 6, 8, 10), las = 1)
        axis(2, xpd = TRUE, yaxt = "s", at = 0.1 * 0:10,
          las = 1, labels = FALSE, tcl = -0.25)
      }
    }
    mtext(preds_set, side = 3, line = 3, cex = 1.25, outer = TRUE)
  } else {
    lev1 <- unique(cells_set[ , preds_set[1]])
    lev2 <- unique(cells_set[ , preds_set[2]])
    for (pr1 in lev1){
      for (pr2 in lev2){
        ci <- paste0(pr1, ".", pr2)
        y_spc <- gdat_spc[ci, "pmed"] * gdat_spc[ci, "kmed"] ^ (0:(nsearch - 1))
        if (!is.null(gdat_ref)){
          y_ref <- gdat_ref[ci, "pmed"] * gdat_ref[ci, "kmed"] ^ (0:(nsearch - 1))
        } else {
          y_ref <- NULL
        }
        SEpanel(found = found[ci, ], available = available[ci, ],
          y_spc = y_spc, y_ref = y_ref, xends = xends, cols_SE)
        if (pr1 == lev1[1])
          mtext(pr2, side = 3, xpd = TRUE)
        if (pr1 == lev1[length(lev1)])
          axis(1, at = 1:ncol(available), xpd = TRUE, xaxt = "s")
        if (pr2 == lev2[length(lev2)]){
          mtext(pr1, side = 4, line = 1.1, xpd = TRUE, las = 0)
        }
        if (pr2 == lev2[1]){
          axis(2, xpd = TRUE, yaxt = "s", at = 0.1 * c(0, 2, 4, 6, 8, 10), las = 1)
          axis(2, xpd = TRUE, yaxt = "s", at = 0.1 * 0:10,
            las = 1, labels = FALSE, tcl = -0.25)
        }
      }
    }
    mtext(preds_set[2], side = 3, line = 2, outer = TRUE, cex = 1.3, las = 0)
    mtext(preds_set[1], side = 4, line = 3.6, outer = TRUE, cex = 1.3, las = 0)
  }
  mtext("Search", side = 1, line = 3.5, outer = TRUE, cex = 1.5)
  mtext("Searcher Efficiency", side = 2, line = 3, outer = TRUE, cex = 1.5)

  # boxes
  par(mfrow = c(1, 1))
  par(omd = c(0, 0.475, 1 - .res * (.header + .box_H)/H, 1 - .res * .header/H),
      mar = c(2.5, 4, 1, 0.5), fg = "black",
      xaxt = "s", yaxt = "s",
      new = T)
  SEboxes("p", cells_set = cells_set, gdat_spc = gdat_spc, gdat_ref = gdat_ref,
    cols_SE = cols_SE)

  par(omd = c(0.475, 0.95, 1 - .res * (.header + .box_H)/H, 1 - .res * .header/H),
      mar = c(2.5, 4, 1, 0.5), fg = "black",
      xaxt = "s", yaxt = "s",
      new = T)
  SEboxes("k", cells_set = cells_set, gdat_spc = gdat_spc, gdat_ref = gdat_ref,
    cols_SE = cols_SE)

  # header
  par(omd = c(0, 1, 1 - .res * .header/H, 1), mar = c(1, 0, 0.5, 0), new = TRUE)
  plot(0, type = 'n', bty = 'n', xlab = "", axes = FALSE,
    ylab = "", ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")

  rect(0.01, 0.7, 0.06, 0.9, lwd = 2, col = cols_SE["spc"], border = NA)
  text(x = 0.07, y = 0.8, adj = 0, cex = 0.7,
    paste0("= Selected Model: ", gsub("~ 1", "~ constant", specificModel)))
  if (!is.null(gdat_ref)){
    rect(0.01, 0.35, 0.06, 0.55, lwd = 2, col = cols_SE["ref"], border = NA)
    text(x = 0.07, y = 0.45, adj = 0, cex = 0.7,
      paste0("= Reference Model: ", gsub("~ 1", "~ constant", referenceModel)))
  }
  text(x = 0.07, y = 0.4 * is.null(gdat_ref), cex = 0.6, adj = 0, xpd = TRUE,
    paste0("Box plots show median, IQR, and 99% &  ", 100 * CL, "% CIs"))
}

#' @title Error check a specific model selection for an SE plot
#'
#' @description Make sure it's available and good, update the name for usage
#'
#' @param modelSet pk model set of class pkmSet
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @return updated name of the model to use
#'
#' @export
#'
checkSpecificModelSE <- function(modelSet, specificModel){
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
    }
    if (any(specificModel %in% names(modelSet)) == FALSE){
      stop("Selected model not in set. To see options use names(modelSet).")
    }
    modNames <- specificModel
    for (modi in modNames){
      if (pkmFail(modelSet[[modi]])){
        stop("specificModel ", modi, " is not a well-fit pk model")
      }
    }
  } else{
    specificModel <- names(pkmSetFailRemove(modelSet))
  }
  return(specificModel)
}

#' @title Tidy an SE model set
#'
#' @description Remove bad fit models
#'
#' @param modelSet pk model set of class pkmSet
#'
#' @return a trimmed model set
#'
#' @export
#'
tidyModelSetSE <- function(modelSet){
  modelSet <- pkmSetFailRemove(modelSet)
  modelSet <- modelSet[order(sapply(modelSet, "[[", "AICc"))]
  class(modelSet) <- c("pkmSet", "list")
  return(modelSet)
}
#' @title Produce a single panel in an SE summary/diagnostic plot
#'
#' @description Each call to \code{SEpanel} produces a single panel showing
#'  searcher efficiency as a function of number of searches. Includes raw data
#'  (\code{found} and \code{available}) and model fits for a specific model
#'  (\code{y_spc}) and for the reference model (\code{y_ref}) for the
#'  \code{pkmSet} object from the reference model was extracted. For interal use
#'  only, for producing figs for \code{plot.pkmSet}. 
#'
#' @param found vector of number carcasses found on the ith attempt
#'
#' @param available vector of number carcasses found on the ith attempt
#'
#' @param y_spc vector of model fits for the specific model
#'
#' @param y_ref vector of model fits for the reference model
#'
#' @param xends x-axis buffer (numeric scalar) on sides of figs
#'
#' @param cols_SE named vector of colors (character)
#'
#' @return NULL inserts a panel with no labels into a preformatted figure
#'
SEpanel <- function(found, available, y_spc, y_ref, xends, cols_SE){
  x <- 1:length(available)
  y <- found/available
  plot(0, ylim = c(0, 1.07), xlim = range(x) + xends * c(-1, 1),
    xlab = "", ylab = "", axes = FALSE, type = "n")
  if (!is.null(y_ref)) lines(x, y_ref, lwd = 3, col = cols_SE["ref"])
  lines(x, y_spc, lwd = 3, col = cols_SE["spc"])
  points(x, y)
  labels <- paste0(found, "/", available)
  pos <- ifelse(y > 0.85, 1, 3)
  pos[is.na(pos)] <- 1
  for(i in 2:length(y)){
    if(is.na(y[i])) break
    if(abs(y[i] - y[i - 1]) < 0.1 & available[i - 1] >= 10){
      if (y[i - 1] <= y[i] & y[i] <= 0.9) {
        pos[c(i - 1, i)] <- c(3, 1)
        next
      }
      if (y[i - 1] <= y[i] & y[i] > 0.9) {
        pos[c(i - 1, i)] <- c(3, 1)
        next
      }
    }
  }
  text(x, y, labels = labels, pos = pos)
  axis(1, at = x, labels = FALSE)
  axis(2, labels = FALSE)
  box()
}

#' @title Produce boxplots \code{p} and/or \code{k} for all cells for reference
#'  model and specific model
#'
#' @description Each call to \code{SEboxes} produces a series of boxplots for
#'  either \code{p} or \code{k}. For interal use only, for producing figs for
#'  \code{plot.pkmSet}. Function requires that ncell, preds_set, lev1, lev2,
#'  boxW, bsep, gdat_spc, and gdat_ref are defined prior to the function call.
#'
#' @param pk either \code{p} or \code{k}, depending on which type of boxplots
#'  need to be inserted
#'
#' @return NULL inserts a boxplot panel into pkmSet plot
#' @keywords internal
SEboxes <- function(pk, cells_set, gdat_spc, gdat_ref, cols_SE = cols_SE){
  if (!is.null(gdat_ref)){
    boxW <- 0.2
    bsep <- 0.05
    mrt <- 0
  } else {
    boxW <- 0.5
    bsep <- 0
    mrt <- boxW/2
  }
  ncell <- nrow(cells_set)
  preds_set <- colnames(cells_set)[-ncol(cells_set)]
  par(fg = "black")
  plot(0, xlim = c(0.5, ncell + 0.5), ylim = c(0, 1),  ylab = "", xlab = "",
    axes = FALSE, mgp = c(1.5, 1, 0),  type = "n", xaxs = "i")
  if (ncell > 1){
    text(1:ncell + 0.2, -0.09, srt = 35, labels = cells_set[ , 1],
      adj = 1, xpd = T, cex = 0.6 + 0.4 * (ncell <= 6))
    lev1 <- unique(cells_set[ , preds_set[1]])
  }
  if (length(preds_set) == 2){
    lev2 <- unique(cells_set[ , preds_set[2]])
    #mtext(preds_set[2], side = 3, line = 1.2, xpd = TRUE)
    mtext(lev2, side = 3, cex = 0.9,
      at = (1:length(lev2) * length(lev1)) + 0.5 - length(lev1)/2)
  }
  #  axis(1, at = 1:ncell, labels = cells_set[ , 1], cex.axis = 0.8, tck = 0)
  for (ci in 1:ncell){
    rect(ci - 0.5, par("usr")[3], ci + 0.5, par("usr")[4],
      col = ifelse(ci%%2 == 0, .cols_div , NA), border = NA)
    par(fg = colors()[125])
    rect(ci - boxW - bsep + mrt, gdat_spc[ci, paste0(pk, "qr1")],
      ci - bsep + mrt, gdat_spc[ci, paste0(pk, "qr3")],
      col = cols_SE["spc"])
    lines(c(ci - boxW - bsep, ci - bsep) + mrt, rep(gdat_spc[ci, paste0(pk, "med")], 2))
    lines(rep(ci, 2) - bsep - boxW/2 + mrt, gdat_spc[ci, paste0(pk, c("bot", "qr1"))])
    lines(rep(ci, 2) - bsep - boxW/2 + mrt, gdat_spc[ci, paste0(pk, c("top", "qr3"))])
    lines(ci - boxW/2 - bsep + boxW/4 * c(-1, 1) + mrt, rep(gdat_spc[ci, paste0(pk, "upr")], 2))
    lines(ci - boxW/2 - bsep + boxW/4 * c(-1, 1) + mrt, rep(gdat_spc[ci, paste0(pk, "lwr")], 2))
    if (!is.null(gdat_ref)){
      par(fg = colors()[218])
      rect(ci + bsep, gdat_ref[ci, paste0(pk, "qr1")],
        ci + boxW + bsep, gdat_ref[ci, paste0(pk, "qr3")],
        col = cols_SE["ref"])
      lines(c(ci + boxW + bsep, ci + bsep), rep(gdat_ref[ci, paste0(pk, "med")], 2))
      lines(rep(ci + boxW/2 + bsep, 2), gdat_ref[ci, paste0(pk, c("bot", "qr1"))])
      lines(rep(ci + boxW/2 + bsep, 2), gdat_ref[ci, paste0(pk, c("top", "qr3"))])
      lines(ci + boxW/2 + bsep + boxW/4 * c(-1, 1), rep(gdat_ref[ci, paste0(pk, "upr")], 2))
      lines(ci + boxW/2 + bsep + boxW/4 * c(-1, 1), rep(gdat_ref[ci, paste0(pk, "lwr")], 2))
    }
    par(fg = "black")
  }

  box()
  axis(2, cex.axis = 0.8, las = 1)
  axis(2, at = (1 + 0:4*2) * 0.1, labels = FALSE, tcl = -0.25)
  mtext(pk, side = 2, line = 2.5, las = 1)
  if (length(preds_set) == 2){
    do.call(abline, list(v = 0.5 + which((1:ncell) %% length(lev1) == 0)))
  }
}

#' Plot parameter box plots for each cell for either p or k
#'
#' @description Boxplot for pk model cells (soon to be deprecated)
#'
#' @param model model of class pkm
#'
#' @param pk character of "p" or "k" to delineate between parameter graphed
#'
#' @param col color to use
#'
#' @return a parameter plot panel
#'
#' @export
#'
pkmParamPlot <- function(model, pk = "p", col){
  ncell <- model$ncell
  cellNames <- model$cells[ , "CellNames"]
  predictors <- model$predictors
  CL <- model$CL
  probs <- c(0, (1 - CL) / 2, 0.25, 0.5, 0.75, 1 - (1 - CL) / 2, 1)
  pks <- rpk(n = 1000, model = model)
  pks_full <- rpk(n = 1000, model = model)

  if (pk == "p"){
    maxy <- 1
  } else if (pk == "k"){
    maxy <- 1
    for (celli in 1:ncell){
      maxcell <- max(pks[[celli]][ , "k"]) * 1.01
      maxy <- max(maxy, maxcell)
    }
  }
  maxy[is.na(maxy)] <- 1

  plot(1, type = "n", xlab = "", ylab = "", bty = "L", xaxt = 'n', yaxt = 'n',
    ylim = c(0, maxy), xlim = c(0.5, ncell + 0.5)
  )

  for (celli in 1:ncell){
    x <- celli
    y <- quantile(pks[[celli]][ , pk], probs, na.rm = TRUE)

    med <- c(-0.1, 0.1)
    tb <- c(-0.07, 0.07)

    rect(x - 0.1, y[3], x + 0.1, y[5], lwd = 2, border = col)
    points(x + med, rep(y[4], 2), type = 'l', lwd = 2, col = col)
    points(x + tb, rep(y[2], 2), type = 'l', lwd = 2, col = col)
    points(x + tb, rep(y[6], 2), type = 'l', lwd = 2, col = col)
    points(rep(x, 3), y[1:3], type = 'l', lwd = 2, col = col)
    points(rep(x, 3), y[5:7], type = 'l', lwd = 2, col = col)

  }

  axis(1, at = 1:ncell, cellNames, las = 1, cex.axis = 0.75)
  axis(2, at = seq(0, 1, 0.5), las = 1, cex.axis = 0.75)
  axis(2, at = seq(0, 1, 0.1), labels = FALSE, tck = -0.05)
  mtext(side = 2, pk, line = 2.75, cex = 1.1)
}

#' Plot cell-specific decay curve for searcher efficiency
#'
#' @param model model of class pkm
#'
#' @param specificCell name of the specific cell to plot (soon to be deprecated)
#'
#' @param col color to use
#'
#' @param axis_x logical of whether or not to plot the x axis
#'
#' @param axis_y logical of whether or not to plot the y axis
#'
#' @return a cell plot panel
#'
#' @export
#'
pkmSECellPlot <- function(model, specificCell, col, axis_y = TRUE,
                          axis_x = TRUE){
  CL <- model$CL
  cellwise <- model$cell_pk
  cellNames <- model$cells[ , "CellNames"]

  whichCarcs <- which(model$carcCell == specificCell)
  observations <- as.matrix(model$observations[whichCarcs, ],
                    nrow = length(whichCarcs), ncol = ncol(model$observations)
                  )
  nobs <- ncol(observations)
  ncarc <- nrow(observations)
  carcFound <- apply(observations, 2, sum, na.rm = TRUE)
  carcUnavail <- apply(apply(observations, 2, is.na), 2, sum)
  carcAvail <- ncarc - carcUnavail

  whichSpecificCell <- which(cellNames == specificCell)
  p <- cellwise[whichSpecificCell, "p_median"]
  k <- cellwise[whichSpecificCell, "k_median"]
  pks <- rpk(n = 1000, model = model)
  ps <- pks[[whichSpecificCell]][ , "p"]
  ks <- pks[[whichSpecificCell]][ , "k"]
  searchTab <- matrix(1:nobs, nrow = length(ps), ncol = nobs, byrow = TRUE)
  ktab <- ks^(searchTab - 1)
  SE <- ps * ktab
  y <- apply(SE, 2, median)
  y_l <- apply(SE, 2, quantile, probs = (1 - CL) / 2)
  y_u <- apply(SE, 2, quantile, probs = 1 - (1 - CL) / 2)

  x_pts <- 1:nobs
  y_pts <- carcFound / carcAvail
  plot(x_pts, y_pts, ylim = c(0, 1), xlim = c(0.5, nobs + 0.5), xlab = "",
    ylab = "", xaxt = "n", yaxt = "n", bty = "L", col = rgb(0.02, 0.02, 0.02),
    lwd = 2, pch = 1, cex = 1.5
  )

  points(x_pts, y, type = 'l', lwd = 3, col = col)
  points(x_pts, y_l, type = 'l', lwd = 2, lty = 3, col = col)
  points(x_pts, y_u, type = 'l', lwd = 2, lty = 3, col = col)

  for (obi in 1:nobs){
    x1 <- x_pts[obi] - 0.25
    y1 <- y_pts[obi] + 0.06
    x2 <- x_pts[obi] + 0.35
    y2 <- y_pts[obi] + 0.15
    rect(x1, y1, x2, y2, border = NA, col = "white")
  }
  obsLabels <- paste(carcFound, carcAvail, sep = "/")
  text(x_pts + 0.05, y_pts + 0.1, obsLabels, cex = 0.65)

  axis(1, at = x_pts, las = 1, cex.axis = 0.75, labels = axis_x)
  axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.75, labels = axis_y)

  text(0.5, 0.95, specificCell, adj = 0, cex = 0.75, font = 2)
}

#' Plot cell-specific decay curve for searcher efficiency for a specific model
#'   with comparison to the cellwise model
#'
#' @param modelSet modelSet of class pkmSet (soon to be deprecated)
#'
#' @param specificModel name of the specific submodel to plot
#'
#' @param specificCell name of the specific cell to plot
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @param axes named vector of logical values indicating whether or not to
#'   plot the x axis and the y axis
#'
#' @return a specific cell plot panel
#'
#' @export
#'
pkmSetSpecSECellPlot <- function(modelSet, specificModel, specificCell,
                                 cols, axes){
  model_spec <- modelSet[[specificModel]]
  model_ref <- refMod(modelSet)

  cellwise_spec <- model_spec$cell_pk
  cellwise_ref <- model_ref$cell_pk
  cellNames_spec <- model_spec$cells[ , "CellNames"]
  cellNames_ref <- model_ref$cells[ , "CellNames"]

  cells_set <- modelSetCells(modelSet)
  preds_set <- modelSetPredictors(modelSet)
  carcCells <- apply(data.frame(model_spec$data0[ , preds_set]),
                 1, paste, collapse = "."
               )
  whichCarcs <- which(carcCells == specificCell)
  if (specificCell == "all"){
    whichCarcs <- 1:length(carcCells)
  }

  observations <- as.matrix(model_spec$observations[whichCarcs, ],
                    nrow = length(whichCarcs),
                    ncol = ncol(model_spec$observations)
                  )
  nobs <- ncol(observations)
  ncarc <- nrow(observations)
  carcFound <- apply(observations, 2, sum, na.rm = TRUE)
  carcUnavail <- apply(apply(observations, 2, is.na), 2, sum)
  carcAvail <- ncarc - carcUnavail

  if (any(grepl("k not estimated", specificModel))){
    return(1)
  }
  pks_spec <- rpk(n = 1000, model = model_spec)
  pks_ref <- rpk(n = 1000, model = model_ref)

#  kIncluded <- !any(grepl("k not estimated", specificModel))
#  if (kIncluded){
#    pks_spec <- rpk(n = 1000, model = model_spec)
#    pks_ref <- rpk(n = 1000, model = model_ref)
#  } else{
#    pks_spec <- rpk(n = 1000, model = model_spec, kFill = 1)
#    pks_ref <- rpk(n = 1000, model = model_ref, kFill = 1)
#  }
  cellMatch_spec <- matchCells(model_spec, modelSet)
  cellMatch_ref <- matchCells(model_ref, modelSet)

  cells_set <- modelSetCells(modelSet)
  cellNames_set <- cells_set$CellNames

  whichSpecificCell_spec <- cellMatch_spec[cellNames_set == specificCell]
  whichSpecificCell_ref <- cellMatch_ref[cellNames_set == specificCell]

  ps_spec <- pks_spec[[whichSpecificCell_spec]][ , "p"]
  ks_spec <- pks_spec[[whichSpecificCell_spec]][ , "k"]
  ps_ref <- pks_ref[[whichSpecificCell_ref]][ , "p"]
  ks_ref <- pks_ref[[whichSpecificCell_ref]][ , "k"]
  searchTab <- matrix(1:nobs, nrow = 1000, ncol = nobs, byrow = TRUE)
  ktab_spec <- ks_spec^(searchTab - 1)
  ktab_ref <- ks_ref^(searchTab - 1)
  SE_spec <- ps_spec * ktab_spec
  SE_ref <- ps_ref * ktab_ref
  y_spec <- apply(SE_spec, 2, median)
  y_ref <- apply(SE_ref, 2, median)

  x_pts <- 1:nobs
  y_pts <- carcFound / carcAvail
  plot(x_pts, y_pts, ylim = c(0, 1.1), xlim = c(0.5, nobs + 0.5), xlab = "",
    ylab = "", xaxt = "n", yaxt = "n", bty = "L", col = rgb(0.02, 0.02, 0.02),
    lwd = 2, pch = 1, cex = 1.5
  )

  points(x_pts, y_ref, type = 'l', lwd = 3, col = cols["ref"])
  points(x_pts, y_spec, type = 'l', lwd = 3, col = cols["spec"])

  for (obi in 1:nobs){
    x1 <- x_pts[obi] - 0.25
    y1 <- y_pts[obi] + 0.035
    x2 <- x_pts[obi] + 0.35
    y2 <- y_pts[obi] + 0.11
    rect(x1, y1, x2, y2, border = NA, col = "white")
  }
  obsLabels <- paste(carcFound, carcAvail, sep = "/")
  text(x_pts + 0.05, y_pts + 0.075, obsLabels, cex = 0.65)

  axis(1, at = x_pts, las = 1, cex.axis = 0.75, labels = axes["x"])
  axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 0.75, labels = axes["y"])
  text(0.5, 1.1, specificCell, adj = 0, cex = 0.75, font = 2, xpd = TRUE)
}
 #' @title p or k box plots for an SE model set
#'
#' @description Plot parameter box plots for each cell within a model for
#'   either p or k with comparison to the cellwise model (soor to be deprecated)
#'
#' @param modelSet modelSet of class pkmSet
#'
#' @param specificModel name of the specific submodel to plot
#'
#' @param pk character of "p" or "k" to delineate between parameter graphed
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @return a specific parameter plot panel
#'
#' @export
#'
pkmSetSpecParamPlot <- function(modelSet, specificModel, pk = "p", cols){
  model_spec <- modelSet[[specificModel]]
  model_ref <- refMod(modelSet)

  CL <- model_ref$CL
  probs <- c(0, (1 - CL) / 2, 0.25, 0.5, 0.75, 1 - (1 - CL) / 2, 1)
  observations_spec <- model_spec$observations
  observations_ref <- model_ref$observations
  ncell_spec <- model_spec$ncell
  ncell_ref <- model_ref$ncell
  cellNames_ref <- model_ref$cells[ , "CellNames"]
  predictors_spec <- model_spec$predictors
  predictors_ref <- model_ref$predictors

  if (any(grepl("k not estimated", specificModel))){
    return(1)
  }
  pks_spec <- rpk(n = 1000, model = model_spec)
  pks_ref <- rpk(n = 1000, model = model_ref)

#  kIncluded <- !any(grepl("k not estimated", specificModel))
#  if (kIncluded){
#    pks_spec <- rpk(n = 1000, model = model_spec)
#    pks_ref <- rpk(n = 1000, model = model_ref)
#  } else{
#    pks_spec <- rpk(n = 1000, model = model_spec)
#    pks_ref <- rpk(n = 1000, model = model_ref)
#  }
  cellMatch_spec <- matchCells(model_spec, modelSet)
  cellMatch_ref <- matchCells(model_ref, modelSet)

  cells_set <- modelSetCells(modelSet)
  cellNames_set <- cells_set$CellNames
  ncell_set <- nrow(cells_set)

  if (pk == "p"){
    maxy <- 1
  } else if (pk == "k"){
    maxy <- 1
    for (celli in 1:ncell_set){
      maxcell_spec <- max(pks_spec[[cellMatch_spec[celli]]][ , "k"]) * 1.01
      maxcell_ref <- max(pks_ref[[cellMatch_ref[celli]]][ , "k"]) * 1.01
      maxy <- max(c(maxy, maxcell_spec, maxcell_ref))
    }
  }
  maxy[is.na(maxy)] <- 1

  par(mar = c(4, 3, 2, 1))
  plot(1, type = "n", xlab = "", ylab = "", bty = "L", xaxt = 'n', yaxt = 'n',
    ylim = c(0, maxy), xlim = c(0.5, ncell_set + 0.5)
  )

  for (celli in 1:ncell_set){
    cMi_s <- cellMatch_spec[celli]
    cMi_f <- cellMatch_ref[celli]
    x_s <- celli - 0.2
    y_s <- quantile(pks_spec[[cMi_s]][ , pk], probs, na.rm = TRUE)
    x_f <- celli + 0.2
    y_f <- quantile(pks_ref[[cMi_f]][ , pk], probs, na.rm = TRUE)

    med <- c(-0.1, 0.1)
    tb <- c(-0.07, 0.07)

    rect(x_s - 0.1, y_s[3], x_s + 0.1, y_s[5], lwd = 1, border = cols["spec"])
    points(x_s + med, rep(y_s[4], 2), type = 'l', lwd = 1, col = cols["spec"])
    points(x_s + tb, rep(y_s[2], 2), type = 'l', lwd = 1, col = cols["spec"])
    points(x_s + tb, rep(y_s[6], 2), type = 'l', lwd = 1, col = cols["spec"])
    points(rep(x_s, 3), y_s[1:3], type = 'l', lwd = 1, col = cols["spec"])
    points(rep(x_s, 3), y_s[5:7], type = 'l', lwd = 1, col = cols["spec"])

    rect(x_f - 0.1, y_f[3], x_f + 0.1, y_f[5], lwd = 1, border = cols["ref"])
    points(x_f + med, rep(y_f[4], 2), type = 'l', lwd = 1, col = cols["ref"])
    points(x_f + tb, rep(y_f[2], 2), type = 'l', lwd = 1, col = cols["ref"])
    points(x_f + tb, rep(y_f[6], 2), type = 'l', lwd = 1, col = cols["ref"])
    points(rep(x_f, 3), y_f[1:3], type = 'l', lwd = 1, col = cols["ref"])
    points(rep(x_f, 3), y_f[5:7], type = 'l', lwd = 1, col = cols["ref"])
  }

  axis(1, at = 1:ncell_set, labels = FALSE, tck = -0.05)
  ang <- 0
  offy <- -0.25
  offx <- NULL
  if (ncell_set > 3){
    ang <- 35
    offx <- 1
  }
  xcex <- 0.75
  if (ncell_set > 6){
    xcex <- 0.5
    offy <- -0.125
  }
  text(1:ncell_set, offy, srt = ang, adj = offx, labels = cellNames_set,
    xpd = TRUE, cex = xcex
  )

  axis(2, at = seq(0, 1, 0.5), las = 1, cex.axis = 0.7)
  axis(2, at = seq(0, 1, 0.1), labels = FALSE, tck = -0.015)
  mtext(side = 2, pk, line = 2.2, cex = 1.125)
}

#' @title p and k box plots for an SE model set
#'
#' @description Plot parameter box plots for each cell within a model for
#'   both p and k with comparison to the cellwise model (soon to be deprecated)
#'
#' @param modelSet modelSet of class pkmSet
#'
#' @param specificModel name of the specific submodel to plot
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @return a set of parameter plot panels
#'
#' @export
#'
plotSEBoxPlots <- function(modelSet, specificModel, cols){
  par(mar = c(0, 0, 0, 0))
  par(fig = c(0, 0.45, 0.7, 0.965), new = TRUE)
  pkmSetSpecParamPlot(modelSet, specificModel, "p", cols)
  par(fig = c(0.45, 0.9, 0.7, 0.965), new = TRUE)
  if (!grepl("k not estimated", specificModel)){
    pkmSetSpecParamPlot(modelSet, specificModel, "k", cols)
  }
}

#' @title template box plot
#'
#' @description Plot template box plot (soon to be deprecated)
#'
#' @param modelSet modelSet of class pkmSet
#'
#' @param specificModel name of the specific submodel to plot
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @return a template box plot
#'
#' @export
#'
plotSEBoxTemplate <- function(modelSet, specificModel, cols){
  model_spec <- modelSet[[specificModel]]
  col_spec <- cols["spec"]
  par(mar = c(0, 0, 0, 0))
  par(fig = c(0.92, 1, 0.8, 0.95), new = TRUE)
  plot(1,1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "",
    ylab = "", ylim = c(0, 1), xlim = c(0, 1)
  )
  x_s <- 0.1
  CL_split <- (1 - model_spec$CL) / 2
  probs_y <- c(0, CL_split, 0.25, 0.5, 0.75, 1 - CL_split, 1)
  set.seed(12)
  y_s <- quantile(rnorm(1000, 0.5, 0.15), probs = probs_y)
  med <- c(-0.1, 0.1)
  tb <- c(-0.07, 0.07)
  rect(x_s - 0.1, y_s[3], x_s + 0.1, y_s[5], lwd = 1, border = col_spec)
  points(x_s + med, rep(y_s[4], 2), type = 'l', lwd = 1, col = col_spec)
  points(x_s + tb, rep(y_s[2], 2), type = 'l', lwd = 1, col = col_spec)
  points(x_s + tb, rep(y_s[6], 2), type = 'l', lwd = 1, col = col_spec)
  points(rep(x_s, 3), y_s[1:3], type = 'l', lwd = 1, col = col_spec)
  points(rep(x_s, 3), y_s[5:7], type = 'l', lwd = 1, col = col_spec)
  num_CL <- c(CL_split, 1 - CL_split) * 100
  text_CL <- paste(num_CL, "%", sep = "")
  text_ex <- c("min", text_CL[1], "25%", "50%", "75%", text_CL[2], "max")
  text(x_s + 0.2, y_s, text_ex, cex = 0.6, adj = 0)
}

#' @title Plot the cellwise results of a single model in a set of SE models
#'
#' @description Produce a set of cellwise figures for a specific SE model, as
#'   fit by \code{\link{pkmSet}} (soon to be deprecated)
#'
#' @param modelSet pk model set of class pkmSet
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @return a plot
#'
#' @export
#'
plotSECells <- function(modelSet, specificModel, cols){
  model_ref <- refMod(modelSet)

  par(fig = c(0, 1, 0, 0.65), new = TRUE, mar = c(1, 1, 1, 1))
  plot(1,1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "",
    ylab = ""
  )
  mtext(side = 1, "Search", line = -0.25, cex = 1.5)
  mtext(side = 2, "Searcher Efficiency", line = -0.25, cex = 1.5)

  cells_set <- modelSetCells(modelSet)
  ncell <- nrow(cells_set)
  cellNames <- cells_set[ , "CellNames"]

  nmatrix_col <- min(3, ncell)
  nmatrix_row <- ceiling(ncell / nmatrix_col)
  figxspace <- 0.95 / nmatrix_col
  figyspace <- 0.65 / nmatrix_row
  x1 <- rep(figxspace * ((1:nmatrix_col) - 1), nmatrix_row) + 0.035
  x2 <- rep(figxspace * ((1:nmatrix_col)), nmatrix_row) + 0.035
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
    pkmSetSpecSECellPlot(modelSet, specificModel, specificCell, cols, axes)
  }
}
#' @title Plot results of a single SE model in a set
#'
#' @description Produce a figures for a specific SE model, as fit by
#'   \code{\link{pkmSet}} (soon to be deprecated)
#'
#' @param modelSet pk model set of class \code{pkmSet}
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @param app logical indicating if the plot is for the app
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @return a plot
#'
#' @export
#'
plotSEFigure <- function(modelSet, specificModel, app, cols){
  plotSEHeader(modelSet, specificModel, app, cols)
  plotSEBoxPlots(modelSet, specificModel, cols)
  plotSEBoxTemplate(modelSet, specificModel, cols)
  plotSECells(modelSet, specificModel, cols)
}

#' @title The SE plot header
#'
#' @description Produce the header for an SE plot (soon to be deprecated)
#'
#' @param modelSet pk model set of class pkmSet
#'
#' @param specificModel the name of the specific model for the plot
#'
#' @param app logical indicating if the plot is for the app
#'
#' @param cols named vector of colors to use for the specific and reference
#'   models
#'
#' @return a plot
#'
#' @export
#'
plotSEHeader <- function(modelSet, specificModel, app = FALSE, cols = SEcols()){
  par(mar = c(0, 0, 0, 0))
  par(fig = c(0, 1, 0.935, 1))
  plot(1, 1, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "",
    ylab = "", ylim = c(0, 1), xlim = c(0, 1)
  )

  LL <- sapply(modelSet, "[[", "loglik")
  referenceMod <- names(modelSet)[which(LL == max(LL))]

  if (app){
    specificModel <- gsub("~ 1", "~ constant", specificModel)
    referenceMod <- gsub("~ 1", "~ constant", referenceMod)
  }

  rect(0.01, 0.725, 0.06, 0.925, lwd = 2, col = cols["spec"], border = NA)
  text_model <- paste("= Selected Model: ", specificModel, sep = "")
  text(x = 0.07, y = 0.85, text_model, adj = 0, cex = 0.9)

  rect(0.01, 0.325, 0.06, 0.525, lwd = 2, col = cols["ref"], border = NA)
  text_model <- paste("= Reference Model: ", referenceMod, sep = "")
  text(x = 0.07, y = 0.45, text_model, adj = 0, cex = 0.9)

  labelsText <- paste(modelSetPredictors(modelSet), collapse = ".")
  labelsText[labelsText == ""] <- "all"
  text_label <- paste("Labels: ", labelsText, sep = "")
  text(x = 0.9, y = 0.8, text_label, adj = 1, cex = 0.75)
}
#' @title Produce a named vectory with standard SE plot colors
#'
#' @description Produce a named vectory with standard SE plot colors. soon to be
#'  deprecated.
#'
#' @export
#'
SEcols <- function(){
  c(spec = "black", ref = "grey")
}
