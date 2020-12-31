#' @title Create the pretty versions of model and summary tables
#'
#' @description Format reader-friendly versions of results summary tables for
#'  searcher efficiency (GUI display and download), carcass persistence, and
#'  splits.
#'
#' @param modTab model table
#'
#' @param CL Confidence level
#'
#' @return pretty version of the SE model table
#'
#' @export
#'
prettyModTabSE <- function(modTab, CL = 0.90){

  kFit <- any(grepl("k_median", colnames(modTab)))
  
  if (!kFit){
    modTab$k_median <- NA
    modTab$k_lwr <- NA
    modTab$k_upr <- NA
  }

  out <- modTab[ , c("cell", "n", "p_median", "k_median")]
  ncell <- nrow(out)

  for (celli in 1:ncell){
    p_m <- round(modTab[celli, "p_median"], 3)
    p_l <- round(modTab[celli, "p_lwr"], 3)
    p_u <- round(modTab[celli, "p_upr"], 3)
    k_m <- round(modTab[celli, "k_median"], 3)
    k_l <- round(modTab[celli, "k_lwr"], 3)
    k_u <- round(modTab[celli, "k_upr"], 3)
    out[celli, "p_median"] <- paste0(p_m, " [", p_l, ", ", p_u, "]")

    if (is.na(k_m)){
      out[celli, "k_median"] <- ""
    } else{
      out[celli, "k_median"] <- paste0(k_m, " [", k_l, ", ", k_u, "]")
    }
  }

  colnames(out) <- c("Cell", "n", "p", "k")
  return(out)
}

#' @title Create the download version of the Searcher Efficiency model table
#'
#' @description Format a user-friendly version of the parameter table from
#'   a Searcher Efficiency model, based on confidence level of interest
#'
#' @param modTab model table
#'
#' @param CL Confidence level
#'
#' @return download version of the SE model table
#'
#' @export
#'
dlModTabSE <- function(modTab, CL = 0.90){

  kFit <- any(grepl("k_median", colnames(modTab)))
  if (!kFit){
    modTab$k_median <- NA
    modTab$k_lwr <- NA
    modTab$k_upr <- NA
  }

  out <- modTab
  coltypes <- c("Median", (1 - CL)/2, 1 - (1 - CL)/2)
  colnames(out) <- c("Cell", "n", paste0("p_", coltypes),
                     paste0("k_", coltypes))
  return(out)
}

#' @title Create the pretty version of the Carcass Persistence model table
#'
#' @description Format a reader-friendly version of the parameter table from
#'  a carcass persistence model showing CIs for medianCP and for rI's for
#'  intervals of Ir
#'
#' @param modTab \code{descCP} object or NULL
#'
#' @return pretty version of the CP model table in a data frame with point
#'  and interval estimates for medianCP and rI statistics. Output table is
#'  ready for rendering in shiny and posting in the GUI
#'
#' @export
#'
prettyModTabCP <- function(modTab){
  table_CP <- modTab
  if(is.null(table_CP))
    return(data.frame(msg = "Selected model was not successfully fit."))
  if (!"descCP" %in% class(table_CP)) stop("table_CP must be a descCP object")
  # descCP objects are matrices with have 3n named columns and ncell named rows

  modTab <- round(table_CP, 3)
  rcols <- grep("^r\\d", colnames(table_CP)) # indices for columns for r statistics
  rcols_abb <- rcols[!rcols %in% grep("_", colnames(table_CP))] # r1, r3, etc.
  cpcols <- grep("CP", colnames(table_CP))
  out <- data.frame(array(dim = c(nrow(table_CP), 2 + length(rcols_abb))))
  names(out) <- c("n", "medianCP", colnames(table_CP)[rcols_abb])
  rownames(out) <- row.names(modTab)
  out$n <- table_CP[ , "n"]
  for (i in 2:ncol(out)){
    out[ , i] <- paste0(modTab[, 2 + 3*(i - 2)],
      "  [", modTab[ , 3 + 3*(i - 2)], ", ", modTab[ , 4 + 3*(i - 2)], "]")
  }
  return(out)
}


#' @title Create the pretty version of the split summary table
#'
#' @description Format a reader-friendly version of the split summary table
#'   a mortality estimation
#'
#' @param splitSummary a split summary
#'
#' @return split pretty table 
#'
#' @export
#'
prettySplitTab <- function(splitSummary){
  if (!("splitSummary" %in% class(splitSummary))){
    out <- as.matrix(splitSummary, nrow = 1)
  } else if (is.list(splitSummary)){
    out <- NULL
    for (i in 1:length(splitSummary)){
      out <- round(rbind(out, splitSummary[[i]]), 2)
    }
    out <- cbind(rep(names(splitSummary), each = dim(splitSummary[[1]])[1]),
      rownames(out), out)
    rownames(out) <- NULL
    colnames(out) <- c(attr(splitSummary, "vars")[2:1], 
                       colnames(splitSummary[[1]]))
  } else {
    if (!is.matrix(splitSummary)){
      splitSummary <- t(as.matrix(splitSummary))
    }
    out <- cbind(rownames(splitSummary), round(splitSummary, 2))
    rownames(out) <- NULL
    colnames(out) <- c(attr(splitSummary, "vars"), colnames(splitSummary))
  }
  out
}