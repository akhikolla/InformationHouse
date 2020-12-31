#' @title GenEst App Messages
#'
#' @description lists of messages used in the app
#'
#' @param msgs message list
#'
#' @param clear logical indicator if clearing should happen.
#'
#' @param modelType "SE", "CP", "g", or "M"
#'
#' @param rv reactive values list
#'
#' @param type "SE", "CP", "M", "split", or "g"
#'
#' @param mods Set Size list of models
#'
#' @param special indicator of a special type of message
#'
#' @param fracNote the note regarding the input
#'
#' @name app_msg_functions
NULL

#' @rdname app_msg_functions
#'
#'
msgList <- function(){
  list(ModSE = NULL, ModCP = NULL, ModM = NULL, SS = NULL, Modg = NULL)
}

#' @rdname app_msg_functions
#'
clearNotifications <- function(msgs = msgList(), clear = TRUE){
  if(clear){
    if (!is.null(msgs$ModSE)){
      removeNotification(msgs$ModSE)
    }
    if (!is.null(msgs$ModCP)){
      removeNotification(msgs$ModCP)
    }
    if (!is.null(msgs$ModM)){
      removeNotification(msgs$ModM)
    }
    if (!is.null(msgs$SS)){
      removeNotification(msgs$SS)
    }
    if (!is.null(msgs$Modg)){
      removeNotification(msgs$Modg)
    }
  }
}

#' @rdname app_msg_functions
#'
msgModRun <- function(msgs, modelType, clear = TRUE){
  clearNotifications(msgs, clear)
  msg <- NULL
  if (modelType == "SE"){
    msg <- ("Running Searcher Efficiency Model")
  }
  if (modelType == "CP"){
    msg <- ("Running Carcass Persistence Model")
  }
  if (modelType == "g"){
    msg <- ("Running Detection Probability Model")
  }
  if (modelType == "M"){
    msg <- ("Estimating Mortality")
  }
  if(!is.null(msg)){
    return(showNotification(msg, duration = NULL))
  }
}

#' @rdname app_msg_functions
#'
msgModDone <- function(msgs, rv, type = "SE", clear = TRUE){
  clearNotifications(msgs, clear)
  if (type == "SE"){
    if ("error" %in% class(rv$mods_SE)){
      return(msgModFail(rv$mods_SE, type = "SE", special = "error"))
    } else {
      kfix <- which(!is.na(rv$kFixed))
      if (length(kfix) > 0 && any(rv$kFixed[kfix] > 1)){
        return(msgModFail(rv$mods_SE_og, "SE", "NA_kFixed"))
      } else if (length(kfix) > 0 && any(rv$kFixed[kfix] < 0)){
        return(msgModFail(rv$mods_SE_og, "SE", "NA_kFixed"))
      } else if (all(unlist(pkmSetSizeFail(rv$mods_SE_og)))){
        return(msgModFail(rv$mods_SE_og, "SE"))
      } else if(any(unlist(lapply(rv$mods_SE_og, pkmSetAllFail)))){
        return(msgModFail(rv$mods_SE_og, "SE", "size_k"))
      } else {
        return(msgModWarning(rv$mods_SE_og, "SE", rv))
      }
    }
  }
  if (type == "CP"){
    if (all(unlist(cpmSetSizeFail(rv$mods_CP_og)))){
      return(msgModFail(rv$mods_CP_og, "CP"))
    } else{
      return(msgModWarning(rv$mods_CP_og, "CP"))
    }
  }  
  if (type == "M"){
    if (is.null(rv$M)){
      if (!is.null(rv$fracNote)){
        return(msgFracNote(rv$fracNote))
      } else {
        return(msgModFail(mods = rv$M, "M"))
      }
    } else if ("error" %in% class(rv$M)){
      if (!is.null(rv$fracNote)){
        return(msgFracNote(rv$fracNote))
      } else {
        return(msgModFail(mods = rv$M[1], type = "M"))
      }
    }
  }
  if (type == "split"){
    if (rv$nsplit_CO + rv$nsplit_SS > 2 | rv$nsplit_SS > 1){
      return(msgSplitFail("setup"))
    }
    if (is.null(rv$Msplit)){
      return(msgSplitFail("run"))
    }
  }
  if (type == "g"){
    if ((is.null(rv$gGeneric))){
      if (is.null(rv$gGeneric[[1]])){    
        return(msgModFail(rv$gGeneric, "g"))
      }
    }
  }
  NULL
}

#' @rdname app_msg_functions
msgModPartialFail <- function(mods, type = "SE"){

  anyFail <- FALSE
  if (type == "SE"){
    if (any(unlist(pkmSetSizeFail(mods)))){
      anyFail <- TRUE
    }
  }
  if (type == "CP"){
    if (any(unlist(cpmSetSizeFail(mods)))){
      anyFail <- TRUE
    }
  }
  if (!anyFail){
    return(NULL)
  }

  nsizeclass <- length(mods)
  uniquemsgs <- NULL
  for (sci in names(mods)){
    newmsgs <- NULL
    if (type == "SE"){
      failedmods <- which(pkmSetFail(mods[[sci]]))
    }
    if (type == "CP"){
      failedmods <- which(cpmSetFail(mods[[sci]]))
    }
    nfailedmods <- length(failedmods)
    if (nfailedmods > 0){
      for(fmodi in 1:nfailedmods){
        newmsg <- mods[[sci]][[names(failedmods)[fmodi]]]
        if (length(newmsg) == 1){
          newmsg <- gsub("Failed model fit: ", "", newmsg)
        } else{
          newmsg <- "Failed fit for k."
        }
        newmsgs <- unique(c(newmsgs, newmsg))
      }
    }
    uniquemsgs <- unique(c(uniquemsgs, newmsgs))
  }
  paste0(
    "Some models were not successfully fit. Failed models were removed. ",
    paste(uniquemsgs, collapse = " ")
  )
}

#' @rdname app_msg_functions
#'
msgSampleSize <- function(mods){
  cellCounts <- countCarcs(mods)
  minCellCount <- min(na.omit(cellCounts))
  if (minCellCount < 10){
    if (length(cellCounts) == 1){
      return(paste0("Caution: n = ", cellCounts, " trial carcasses may be too ",
        "few for reliably estimating searcher efficiency parameters."))
    }
    return(paste0("Small (< 10) sample sizes in some cells. ",
      "Consider simplifying the model; parameter estimates may be unstable."))
  }
  NULL
}

#' @rdname app_msg_functions
#'
msgModWarning <- function(mods, type = "SE", rv = NULL){
  msg <- paste(msgModPartialFail(mods, type), msgSampleSize(mods), sep = " ")
  if (type == "SE"){
    msg <- paste(msg, msgModSENobs(rv), sep = " ")
  }
  if (length(msg) > 0){
    return(showNotification(msg, type = "warning", duration = NULL))    
  }
  NULL
}

#' @rdname app_msg_functions
#'
msgModSENobs <- function(rv){
  if (length(rv$obsCols_SE) == 1){
    if(length(rv$formula_k) > 0 & any(is.na(rv$kFixed))){
      return("Only one observation column, k not estimated.")
    }
  }
  NULL
}


#' @rdname app_msg_functions
#'
msgModFail <- function(mods, type = "SE", special = NULL){
  if (type %in% c("SE", "CP")){
    if (is.null(special)){
      msg <- paste(
               "No models were successfully fit.", 
                gsub("Failed model fit: ", "", unique(unlist(mods))),
                sep = " "
             )
    } else if (special == "error") {
      msg <- mods
    } else if (special == "size_k"){
      msg <- "Some carcass classes had no successful models. Consider a fixed k."
    } else if (special == "NA_kFixed"){
      msg <- "invalid value entered for fixed k"
    }
  }
  if (type == "g"){
    msg <- "Cannot estimate detection probability"
  }
  if (type == "M"){
    msg <- mods
  }
  if(!is.null(msg)){
    return(showNotification(msg, type = "error", duration = NULL))
  }
}

#' @rdname app_msg_functions
#'
msgSSavgFail <- function(msgs, rv, clear = TRUE){
  if (clear){
    clearNotifications(msgs)
  }
  if (is.na(rv$SStemp[1])){
    msg <- "Search Schedule can't be averaged using date searched column"
    return(showNotification(msg, type = "warning", duration = NULL))
  }
  NULL
}

#' @rdname app_msg_functions
#'
msgSSinputFail <- function(msgs, rv, clear = TRUE){
  if (clear){
    clearNotifications(msgs)
  }
  if (is.na(rv$SStemp[1])){
    msg <- "Search Schedule can't be created using inputs"
    return(showNotification(msg, type = "warning", duration = NULL))
  }
  NULL
}

#' @rdname app_msg_functions
#'
msgSplitFail <- function(type = NULL){

  if (is.null(type)){
    return(NULL)
  }
  if (type == "setup"){
    msg <- paste0(
             "Improper splits setup. Maximum two splits, only one of which ",
             "can be associated with the search schedule."
           )
  }
  if (type == "run"){
    msg <- "Splits calculation failed. Check split selections."
  }
  return(showNotification(msg, type = "error", duration = NULL))
}

#' @rdname app_msg_functions
#'
msgFracNote <- function(fracNote){
  return(showNotification(fracNote, type = "warning", duration = NULL))
}
