#' @title app utilities for formatting text, tables, figs, etc. for display
#'
#' @param rv Reactive values list for the GenEst GUI.
#'
#' @param type Model type, either "SE" or "CP" or "g".
#'
#' @param output \code{output} list to have elements \code{dontSuspend}
#'   (re)set to having \code{suspendWhenHidden = FALSE}.
#'
#' @param dontSuspend Names of elements in \code{output} to (re)set to
#'   having \code{suspendWhenHidden = FALSE}.
#'
#' @param x \code{list} object to have elements \code{toNULL} reset to
#'   \code{NULL}.
#'
#' @param toNULL Names of elements in \code{x} to reset to \code{NULL}.
#'
#' @param modelSet Model set of class \code{cpmSet} or \code{pkmSet}.
#'
#' @param modType "SE" or "CP"
#'
#' @name app_output_utilities
NULL

#' @rdname app_output_utilities
#'
#'
classText <- function(rv, type = "SE"){
  out <- ""
  if (type == "SE"){
    if (length(rv$sizeclasses_SE) > 1){
      out <- paste0("Carcass class: ", rv$sizeclass_SE)
    }        
  }
  if (type == "CP"){
    if (length(rv$sizeclasses_CP) > 1){
      out <- paste0("Carcass class: ", rv$sizeclass_CP)
    }    
  }
  if (type == "g"){
    if (length(rv$sizeclasses_g) > 1){
      out <- paste0("Carcass class: ", rv$sizeclass_g,
        " ........ Search schedule: I = ",
        round(rv$SS[["I"]],1), ", span = ", rv$SS[["span"]])
    }
  }
  renderText(out)
}

#' @rdname app_output_utilities
#'
estText <- function(rv, type = "SE"){
  out <- NULL
  if (type == "SE"){
    out <- paste0("Table shows median estimates and ", 100 * rv$CL,  
             "% confidence intervals")
  }
  if (type == "CP"){
    out <- paste0("Table shows median estimates and ", 100 * rv$CL,  
             "% confidence intervals for median persistence and r statistics")
  }
  renderText(out)
}

#' @rdname app_output_utilities
setNotSuspending <- function(output, dontSuspend){
  for(i in 1:length(dontSuspend)){
    outputOptions(output, dontSuspend[i], suspendWhenHidden = FALSE)
  }
}

#' @rdname app_output_utilities
reNULL <- function(x, toNULL){
  for(i in 1:length(toNULL)){
    x[[toNULL[i]]] <- NULL
  }
  x
}

#' @rdname app_output_utilities
#'
initialOutput <- function(rv, output){
  output$kNeed <- renderText("no")
  outputOptions(output, "kNeed", suspendWhenHidden = FALSE)
  return(output)
}
#' @rdname app_output_utilities
#'
setFigW <- function(modelSet){
  if (!any(attr(modelSet, "class") %in% c("cpmSet", "pkmSet"))){
    stop("modelSet must be a cpmSet or pkmSet object")
  }
  ncell <- nrow(modelSetCells(modelSet))
  if (ncell > 6){
    return(1200)
  } else{
    return(800)
  }
}

#' @rdname app_output_utilities
#'
setFigH <- function(modelSet, modType = "SE"){
  if (!any(attr(modelSet, "class") %in% c("cpmSet", "pkmSet"))){
    stop("modelSet must be a cpmSet or pkmSet object")
  }

  ncell <- nrow(modelSetCells(modelSet))
  nRow <- ceiling(ncell / 3)
  mult <- 200
  if (ncell > 6){
    mult <- 300
  }
  out <- max(c(nRow * mult + 400, 800))
  if ("cpmSet" %in% class(modelSet)) out <- out - 100
  out
}

