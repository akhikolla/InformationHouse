#' @title Update the reactive value list when an event occurs
#'
#' @description When an event occurs in the GenEst GUI, the reactive values 
#'   need to be updated. This function contains all of the possible updates
#'   based on the event options.
#'
#' @param eventName Character name of the event. One of "clear_all",
#'   "file_SE", "file_SE_clear", "file_CP", "file_CP_clear", "file_SS",
#'   "file_SS_clear", "file_DWP", "file_DWP_clear", "file_CO", 
#'   "file_CO_clear", "class", "obsSE", "predsSE", "run_SE", "run_SE_clear",
#'   "outSEclass", "outSEp", "outSEk", "ltp", "fta", "predsCP", "run_CP",
#'   "run_CP_clear", "outCPclass", "outCPdist", "outCPl", "outCPs",
#'   "run_M", "run_M_clear", "split_M", "split_M_clear", "transpose_split",
#'   "run_g", "run_g_clear", "outgclass", "load_RP", "load_RPbat",
#'   "load_cleared", "load_PV", "load_trough", "load_powerTower", or "load_mock"
#'
#' @param rv Reactive values list for the GenEst GUI, created by 
#'   \code{\link{initialReactiveValues}}, which calls 
#'   \code{\link[shiny]{reactiveValues}}
#'
#' @param input \code{input} list for the GenEst GUI.
#'
#' @return Updated \code{rv} list.
#'
update_rv <- function(eventName, rv, input){
  eventOptions <- c("clear_all", "file_SE", "file_SE_clear", "file_CP",
                    "file_CP_clear", "file_SS", "file_SS_clear", "file_DWP",
                    "file_DWP_clear", "file_CO", "file_CO_clear", "class",
                    "obsSE", "predsSE", "run_SE", "run_SE_clear",
                    "outSEclass", "outSEp", "outSEk", "ltp", "fta", "predsCP",
                    "run_CP", "run_CP_clear", "outCPclass", "outCPdist",
                    "outCPl", "outCPs", "run_M", "run_M_clear", "split_M",
                    "split_M_clear", "transpose_split",
                    "run_g", "run_g_clear", "outgclass",
                    "load_RP", "load_RPbat", "load_cleared", "load_PV",
                    "load_trough", "load_powerTower", "load_mock", "cscale")

  if (missing(eventName) || (eventName %in% eventOptions) == FALSE){
    stop("eventName missing or not in list of available eventNames")
  }

  if (eventName == "clear_all"){
    toNULL <- c("data_SE", "filename_SE", "colNames_SE", "colNames_SE_preds",
                "colNames_SE_preds0", "colNames_SE_obs", "colNames_SE_obs0",
                "toRemove_SE_obs", "toRemove_SE_preds", "sizeclass_SE",
                "obsCols_SE", "preds_SE", "predictors_SE", "formula_p", 
                "formula_k", "kFixedChoice", "kFixed", "mods_SE", 
                "mods_SE_og", "sizeclasses_SE", "outSEpk", "AICcTab_SE", 
                "modOrder_SE", "modNames_SE", "modNames_SEp", "modNames_SEk", 
                "modSet_SE", "best_SE", "modTab_SE", "modTabPretty_SE",
                "modTabDL_SE", "SStemp", "avgSI", "splittable_SS",
                "M", "Msplit", "unitCol", "sizeCol_M",
                "split_CO", "split_SS", "SEmodToUse",
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "SEmodToUse_g",
                "data_CP", "filename_CP", "colNames_CP", "colNames_CP_preds",
                "colNames_CP_preds0", "colNames_fta", "colNames_fta0",
                "colNames_ltp", "colNames_ltp0", "colNames_CP_preds", 
                "toRemove_fta", "toRemove_ltp", "toRemove_CP_preds", 
                "sizeclass_CP", "ltp", "fta", "preds_CP", "dist", 
                "predictors_CP", "formula_l", "formula_s", "mods_CP", 
                "mods_CP_og", "sizeclasses_CP", "AICcTab_CP", "modOrder_CP", 
                "modNames_CP", "modNames_CPl", "modNames_CPs", 
                "modNames_CPdist", "modSet_CP", "best_CP", "modTab_CP", 
                "modTabPretty_CP", "modTabDL_CP", "CPmodToUse", 
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "CPmodToUse_g",
                "data_SS", "colNames_SS", "data_DWP", "colNames_DWP",
                "data_CO", "colNames_CO", "colNames_COdates", "filename_SS",
                "filename_DWP", "filename_CO", "sizeCol", "colNames_size",
                "colNames_size0")
    rv <- reNULL(rv, toNULL)

    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "figH_SE", "figW_SE", "SS",
                 "gSearchInterval", "gSearchMax", "figH_CP", "figW_CP")
    rv <- reVal(rv, toReVal)
  }

  if (eventName == "file_SE"){
    toNULL <- c("data_SE", "filename_SE", "colNames_SE", "colNames_SE_preds",
                "colNames_SE_preds0", "colNames_SE_obs", "colNames_SE_obs0",
                "toRemove_SE_obs", "toRemove_SE_preds", "sizeclass_SE",
                "obsCols_SE", "preds_SE", "predictors_SE", "formula_p",
                "formula_k", "kFixedChoice", "kFixed", "mods_SE",
                "mods_SE_og", "sizeclasses_SE", "outSEpk", "AICcTab_SE",
                "modOrder_SE", "modNames_SE", "modNames_SEp", "modNames_SEk",
                "modSet_SE", "best_SE", "modTab_SE", "modTabPretty_SE",
                "modTabDL_SE", "SStemp", "avgSI",
                "M", "Msplit", "unitCol", "sizeCol_M",
                "split_CO", "split_SS", "SEmodToUse", "sizeCol", "sizeclasses",
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)

    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "figH_SE", "figW_SE", "SS")
    rv <- reVal(rv, toReVal)

    rv$data_SE <- try(readCSV(input$file_SE$datapath))
    rv$filename_SE <- input$file_SE$name
    if ("try-error" %in% class(rv$data_SE)){
      rv$data_SE <- data.frame(matrix(rv$data_SE, nrow = 1))
      names(rv$data_SE) <- rv$filename_SE
#      rv$data_SE <- NULL
      rv$filename_SE <- NULL
      rv$colNames_SE <- NULL
      rv$colNames_SE_preds0 <- NULL
      rv$colNames_SE_obs0 <- NULL
      rv$colNames_size0 <- NULL
      rv$colNames_SE_preds <- NULL
      rv$colNames_SE_obs <- NULL
    } else {
      rv$colNames_SE <- colnames(rv$data_SE)
      rv$colNames_SE_preds0 <- predsCols(rv$data_SE)
      rv$colNames_SE_obs0 <- obsCols_SE(rv$data_SE)
      rv$colNames_size0 <- updateColNames_size(rv)
      rv$colNames_SE_preds <- rv$colNames_SE_preds0
      rv$colNames_SE_obs <- rv$colNames_SE_obs0
      rv$colNames_size <- removeCols(rv$colNames_size0,
        c(rv$ltp, rv$fta, rv$preds_CP))
    }
    if (!is.null(rv$sizeCol) && !(rv$sizeCol %in% rv$colNames_size))
      rv$sizeCol <- NULL
  }


  if (eventName == "file_SE_clear"){
    toNULL <- c("data_SE", "filename_SE", "colNames_SE", "colNames_SE_preds",
                "colNames_SE_preds0", "colNames_SE_obs", "colNames_SE_obs0",
                "toRemove_SE_obs", "toRemove_SE_preds", "sizeclass_SE",
                "obsCols_SE", "preds_SE", "predictors_SE", "formula_p", 
                "formula_k", "kFixedChoice", "kFixed", "mods_SE", 
                "mods_SE_og", "sizeclasses_SE", "outSEpk", "AICcTab_SE", 
                "modOrder_SE", "modNames_SE", "modNames_SEp", "modNames_SEk", 
                "modSet_SE", "best_SE", "modTab_SE", "modTabPretty_SE",
                "modTabDL_SE", "SStemp", "avgSI", "M", "Msplit", "unitCol",
                "sizeCol_M", "split_CO", "split_SS", "SEmodToUse",
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)

    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "figH_SE", "figW_SE")
    rv <- reVal(rv, toReVal)
    rv$colNames_size0 <- updateColNames_size(rv)
    rv$colNames_size <- rv$colNames_size0
  }

  if (eventName == "file_CP"){
    toNULL <- c("data_CP", "filename_CP", "colNames_CP", "colNames_CP_preds",
                "colNames_CP_preds0", "colNames_fta", "colNames_fta0",
                "colNames_ltp", "colNames_ltp0",
                "toRemove_fta", "toRemove_ltp", "toRemove_CP_preds",
                "sizeclass_CP", "ltp", "fta", "preds_CP", "dist",
                "predictors_CP", "formula_l", "formula_s", "mods_CP",
                "mods_CP_og", "sizeclasses_CP", "AICcTab_CP", "modOrder_CP",
                "modNames_CP", "modNames_CPl", "modNames_CPs",
                "modNames_CPdist", "modSet_CP", "best_CP", "modTab_CP",
                "modTabPretty_CP", "modTabDL_CP", "SStemp", "avgSI", "M",
                "Msplit", "unitCol", "sizeCol_M", "split_CO", "split_SS",
                "sizeclasses", "sizeCol", "CPmodToUse", "sizeclasses_g",
                "nsizeclasses_g", "gGeneric", "CPmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "figH_CP", "figW_CP")
    rv <- reVal(rv, toReVal)

    rv$data_CP <- try(readCSV(input$file_CP$datapath))
    rv$filename_CP <- input$file_CP$name
    if ("try-error" %in% class(rv$data_CP)){
      rv$data_CP <- data.frame(matrix(rv$data_CP, nrow = 1))
      names(rv$data_CP) <- rv$filename_CP
      rv$filename_CP <- NULL
      rv$colNames_CP <- NULL
      rv$colNames_CP_preds0 <- NULL
      rv$colNames_fta0 <- NULL
      rv$colNames_ltp0 <- NULL
      rv$colNames_size0 <- NULL
      rv$colNames_fta <- NULL
      rv$colNames_ltp <- NULL
      rv$colNames_size <- NULL
    } else {
      rv$colNames_CP <- colnames(rv$data_CP)
      rv$colNames_CP_preds0 <- predsCols(rv$data_CP)
      rv$colNames_fta0 <- obsCols_fta(rv$data_CP)
      rv$colNames_ltp0 <- obsCols_ltp(rv$data_CP)
      rv$colNames_size0 <- updateColNames_size(rv)
      rv$colNames_fta <- rv$colNames_fta0
      rv$colNames_ltp <- rv$colNames_ltp0
      rv$colNames_size <- removeCols(rv$colNames_size0,
        c(rv$obsCols_SE, rv$preds_SE))
      if (!is.null(rv$sizeCol) && !(rv$sizeCol %in% rv$colNames_size))
        rv$sizeCol <- NULL
    }
  }

  if (eventName == "file_CP_clear"){
    toNULL <- c("data_CP", "filename_CP", "colNames_CP", "colNames_CP_preds",
                "colNames_CP_preds0", "colNames_fta", "colNames_fta0",
                "colNames_ltp", "colNames_ltp0", "colNames_CP_preds", 
                "toRemove_fta", "toRemove_ltp", "toRemove_CP_preds", 
                "sizeclass_CP", "ltp", "fta", "preds_CP", "dist", 
                "predictors_CP", "formula_l", "formula_s", "mods_CP", 
                "mods_CP_og", "sizeclasses_CP", "AICcTab_CP", "modOrder_CP", 
                "modNames_CP", "modNames_CPl", "modNames_CPs", 
                "modNames_CPdist", "modSet_CP", "best_CP", "modTab_CP", 
                "modTabPretty_CP", "modTabDL_CP", "SStemp", "avgSI", "M",
                "Msplit", "unitCol", "sizeCol_M", "split_CO", "split_SS",
                "CPmodToUse", "sizeclasses_g", "nsizeclasses_g", "gGeneric", 
                "CPmodToUse_g")
    rv <- reNULL(rv, toNULL)

    rv$colNames_size <- updateColNames_size(rv)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "figH_CP", "figW_CP")
    rv <- reVal(rv, toReVal)
    rv$colNames_size0 <- updateColNames_size(rv)
    rv$colNames_size <- rv$colNames_size0
  }

  if (eventName == "file_SS"){
    toNULL <- c("data_SS", "colNames_SS", "SStemp", "avgSI", "splittable_SS",
                "M", "Msplit", "unitCol", "sizeCol_M", "split_CO", "split_SS",
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "SS", "gSearchInterval", "gSearchMax")
    rv <- reVal(rv, toReVal)
    rv$data_SS <- try(readCSV(input$file_SS$datapath))
    rv$filename_SS <- input$file_SS$name
    if ("try-error" %in% class(rv$data_SS)){
      rv$data_SS <- data.frame(matrix(rv$data_SS, nrow = 1))
      names(rv$data_SS) <- rv$filename_SS
      rv$filename_SS <- NULL
      rv$colNames_SS <- NULL
      rv$splittable_SS <- NULL
      rv$SS <- NULL
      rv$avgSI <- NULL
    } else {
      rv$colNames_SS <- colnames(rv$data_SS)
      rv$splittable_SS <- rv$colNames_SS   # candidates for splittable columns
      badind <- NULL
      for (ci in 1:length(rv$splittable_SS)){
        if (is.numeric(rv$data_SS[ , rv$splittable_SS[ci]])){
          badind <- c(badind, ci) # SS splits are categorical (or dates)
        } else {
          SSlev <- rv$data_SS[ , rv$splittable_SS[ci]]
          if (length(unique(SSlev)) == 1 || # no levels to split on
            min(diff(match(SSlev, unique(SSlev)))) < 0) # contiguous blocks
            badind <- c(badind, ci)
        }
      }
      if (length(badind) > 0) rv$splittable_SS <- rv$splittable_SS[-badind]
      rv$SStemp <- tryCatch(averageSS(rv$data_SS), error = function(x){NA})
      if (!is.na(rv$SStemp[1])){ # aveSS for default SS in estg (if possible)
        rv$SS <- list("span" = max(rv$SStemp), "I" = rv$SStemp[2])
        rv$avgSI <-  mean(diff(rv$SStemp[-length(rv$SStemp)]))
      } else {
        rv$SS <- NULL # no default
      }
    }
  }


  if (eventName == "file_SS_clear"){
    toNULL <- c("data_SS", "filename_SS", "colNames_SS", "SStemp", "avgSI", "M",
                "Msplit", "splittable_SS", "unitCol", "sizeCol_M", "split_CO",
                "split_SS", "SS", "sizeclasses_g", "nsizeclasses_g", "gGeneric",
                "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "gSearchInterval", "gSearchMax")
    rv <- reVal(rv, toReVal)
  }

  if (eventName == "file_DWP"){
    toNULL <- c("data_DWP", "filename_DWP", "colNames_DWP", "M", "Msplit", "unitCol",
                "sizeCol_M", "split_CO", "split_SS")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
    rv$data_DWP <- try(readCSV(input$file_DWP$datapath))
    rv$filename_DWP <- input$file_DWP$name
    if ("try-error" %in% class(rv$data_DWP)){
      rv$data_DWP <- data.frame(matrix(rv$data_DWP, nrow = 1))
      names(rv$data_DWP) <- rv$filename_DWP
      rv$filenames_DWP <- NULL
      rv$colNames_DWP <- NULL
    } else {
      rv$colNames_DWP <- DWPCols(rv$data_DWP)
    }
  }


  if (eventName == "file_DWP_clear"){
    toNULL <- c("data_DWP", "filename_DWP", "colNames_DWP", "M", "Msplit", "unitCol",
                "sizeCol_M", "split_CO", "split_SS")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
  }

  if (eventName == "file_CO"){
    toNULL <- c("data_CO", "filename_CO", "colNames_CO", "colNames_COdates",
      "M", "Msplit", "xID", "unitCol", "sizeCol_M", "SEmodToUse", "split_CO",
      "split_SS")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
    rv$data_CO <- try(readCSV(input$file_CO$datapath))
    rv$filename_CO <- input$file_CO$name
    if ("try-error" %in% class(rv$data_CO)){
      rv$data_CO <- data.frame(matrix(rv$data_CO, nrow = 1))
      names(rv$data_CO) <- rv$filename_CO
      rv$filename_CO <- NULL
      rv$colNames_xID <- NULL
      rv$colNames_CO <- NULL
      rv$colNames_COdates <- NULL
      rv$colNames_size0 <- NULL
    } else {
      rv$data_CO <- readCSV(input$file_CO$datapath)
      rv$filename_CO <- input$file_CO$name
      rv$colNames_xID <- names(which(
        apply(rv$data_CO, FUN = function(x) length(unique(x)), MARGIN = 2) ==
        apply(rv$data_CO, FUN = length, MARGIN = 2)))
      rv$colNames_CO <- colnames(rv$data_CO)
      rv$colNames_COdates <- dateCols(rv$data_CO)
      rv$colNames_size0 <- updateColNames_size(rv)
      rv$colNames_size <- removeCols(rv$colNames_size0,
        c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta, rv$preds_CP))
    }
    if (!is.null(rv$sizeCol) && !(rv$sizeCol %in% rv$colNames_size))
      rv$sizeCol <- NULL
  }

  if (eventName == "file_CO_clear"){
    toNULL <- c("data_CO", "filename_CO", "colNames_CO", "colNames_COdates", "M", "Msplit",
                "unitCol", "sizeCol_M", "SEmodToUse", "split_CO", "split_SS")
    rv <- reNULL(rv, toNULL)
    rv$colNames_size <- updateColNames_size(rv)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
    rv$colNames_size0 <- updateColNames_size(rv)
    rv$colNames_size <- rv$colNames_size0
  }

  if (grepl("load_", eventName)){
    toNULL <- c("data_SE", "filename_SE", "colNames_SE", "colNames_SE_preds",
                "colNames_SE_preds0", "colNames_SE_obs", "colNames_SE_obs0",
                "toRemove_SE_obs", "toRemove_SE_preds", "sizeclass_SE",
                "obsCols_SE", "preds_SE", "predictors_SE", "formula_p",
                "formula_k", "kFixedChoice", "kFixed", "mods_SE",
                "mods_SE_og", "sizeclasses_SE", "outSEpk", "AICcTab_SE",
                "modOrder_SE", "modNames_SE", "modNames_SEp", "modNames_SEk",
                "modSet_SE", "best_SE", "modTab_SE", "modTabPretty_SE",
                "modTabDL_SE", "SStemp", "avgSI", "M", "Msplit", "unitCol",
                "sizeCol_M", "split_CO", "split_SS", "SEmodToUse",
                "splittable_SS", "sizeclasses_g", "nsizeclasses_g", "gGeneric",
                "SEmodToUse_g", "data_CP", "filename_CP", "colNames_CP",
                "colNames_CP_preds", "colNames_CP_preds0", "colNames_fta",
                "colNames_fta0", "colNames_ltp", "colNames_ltp0",
                "colNames_CP_preds", "toRemove_fta", "toRemove_ltp",
                "toRemove_CP_preds", "sizeclass_CP", "ltp", "fta", "preds_CP",
                "dist", "predictors_CP", "formula_l", "formula_s", "mods_CP",
                "mods_CP_og", "sizeclasses_CP", "AICcTab_CP", "modOrder_CP",
                "modNames_CP", "modNames_CPl", "modNames_CPs",
                "modNames_CPdist", "modSet_CP", "best_CP", "modTab_CP",
                "modTabPretty_CP", "modTabDL_CP", "CPmodToUse",
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "CPmodToUse_g",
                "data_SS", "colNames_SS", "data_DWP", "colNames_DWP",
                "data_CO", "colNames_CO", "colNames_COdates", "filename_SS",
                "filename_DWP", "filename_CO")
    rv <- reNULL(rv, toNULL)

    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M", "figH_g",
                 "figW_g", "figH_SE", "figW_SE", "SS",
                 "gSearchInterval", "gSearchMax", "figH_CP", "figW_CP")
    rv <- reVal(rv, toReVal)

    if (grepl("mock", eventName)){
      dataset <- GenEst::mock
      dataName <- "mock"
    } else {
      dataset <- switch(eventName,
        "load_RP" = GenEst::wind_RP,
        "load_RPbat" = GenEst::wind_RPbat,
        "load_cleared" = GenEst::wind_cleared,
        "load_PV" = GenEst::solar_PV,
        "load_trough" = GenEst::solar_trough,
        "load_powerTower" = GenEst::solar_powerTower
      )
      dataName <- switch(eventName,
        "load_RP" = "wind_RP",
        "load_RPbat" = "wind_RPbat",
        "load_cleared" = "wind_cleared",
        "load_PV" = "solar_PV",
        "load_trough" = "solar_trough",
        "load_powerTower" = "solar_powerTower"
      )
    }

    rv$data_SE <- dataset$SE
    rv$filename_SE <- paste0(dataName, "$SE")

    rv$data_CP <- dataset$CP
    rv$filename_CP <- paste0(dataName, "$CP")

    rv$data_SS <- dataset$SS
    rv$filename_SS <- paste0(dataName, "$SS")
    rv$SStemp <- tryCatch(averageSS(rv$data_SS), error = function(x){NA})
    if (!is.na(rv$SStemp[1]) && length(rv$SStemp > 1)){
      rv$avgSI <-  mean(diff(rv$SStemp[-length(rv$SStemp)]))
      if (max(rv$SStemp)%%rv$SStemp[2] == 0){
        rv$SS <- list("span" = max(rv$SStemp), "I" = rv$SStemp[2])
      } else {
        rv$SS <- list("span" = rv$SStemp[2] * length(rv$SStemp), "I" = rv$SStemp[2])
      }
    } else {
      rv$SS <- NULL
    }
    rv$data_DWP <- dataset$DWP
    rv$filename_DWP <- paste0(dataName, "$DWP")

    rv$data_CO <- dataset$CO
    rv$filename_CO <- paste0(dataName, "$CO")

    rv$colNames_CO <- colnames(rv$data_CO)
    rv$colNames_COdates <- dateCols(rv$data_CO)

    rv$colNames_SS <- colnames(rv$data_SS)
    rv$splittable_SS <- rv$colNames_SS
    badind <- NULL
    for (ci in 1:length(rv$splittable_SS)){
      if (is.numeric(rv$data_SS[ , rv$splittable_SS[ci]])){
        badind <- c(badind, ci)
      } else {
        SSlev <- rv$data_SS[ , rv$splittable_SS[ci]]
        if (min(diff(match(SSlev, unique(SSlev)))) < 0)
          badind <- c(badind, ci)
      }
    }
    if (length(badind) > 0) rv$splittable_SS <-   rv$splittable_SS[-badind]
    rv$colNames_DWP <- DWPCols(rv$data_DWP)
    rv$colNames_CP <- colnames(rv$data_CP)
    rv$colNames_CP_preds0 <- predsCols(rv$data_CP)
    rv$colNames_fta0 <- obsCols_fta(rv$data_CP)
    rv$colNames_ltp0 <- obsCols_ltp(rv$data_CP)
    rv$colNames_CP_preds <- rv$colNames_CP_preds0
    rv$colNames_fta <- rv$colNames_fta0
    rv$colNames_ltp <- rv$colNames_ltp0

    rv$colNames_SE <- colnames(rv$data_SE)
    rv$colNames_SE_preds0 <- predsCols(rv$data_SE)
    rv$colNames_SE_obs0 <- obsCols_SE(rv$data_SE)
    rv$colNames_SE_preds <- rv$colNames_SE_preds0
    rv$colNames_SE_obs <- rv$colNames_SE_obs0

    rv$colNames_size0 <- updateColNames_size(rv)
    rv$colNames_size <- rv$colNames_size0
    rv$colNames_xID <- names(which(
      apply(rv$data_CO, FUN = function(x) length(unique(x)), MARGIN = 2) ==
      apply(rv$data_CO, FUN = length, MARGIN = 2)))
    rv$xIDcol <- rv$colNames_xID[1]
    rv$sizeCol <- NULL
  }
  if (eventName == "class"){
    rv$sizeCol <- input$class
    rv$obsCols_SE <- input$obsSE
    rv$preds_SE <- input$predsSE
    rv$ltp <- input$ltp
    rv$fta <- input$fta
    rv$preds_CP <- input$predsCP
    #remove sizeCol from list of possibilities for SE_preds and SE_obs
    rv$colNames_SE_preds <- removeCols(rv$colNames_SE_preds0,
      c(rv$obsCols_SE, rv$sizeCol))
    rv$colNames_SE_obs <- removeCols(rv$colNames_SE_obs0,
      c(rv$preds_SE, rv$sizeCol))

    #remove sizeCol from list of possibilities for CP_preds, ltp, and fta
    rv$colNames_CP_preds <- removeCols(rv$colNames_CP_preds0,
      c(rv$ltp, rv$fta, rv$sizeCol))
    rv$colNames_ltp <- removeCols(rv$colNames_ltp0,
      c(rv$preds_CP, rv$fta, rv$sizeCol))
    rv$colNames_fta <- removeCols(rv$colNames_fta0,
      c(rv$preds_CP, rv$ltp, rv$sizeCol))

    if (!is.null(rv$sizeCol)){
      rv$sizeclasses <- sort(unique(
        c(rv$data_SE[ , rv$sizeCol], rv$data_CP[ , rv$sizeCol])))
      rv$sizeclasses_k <- sort(unique(rv$data_SE[ , rv$sizeCol]))
      rv$sizeclasses_CP <- sort(unique(rv$data_CP[, rv$sizeCol]))
      rv$sizeclasses_SE <- sort(unique(rv$data_SE[, rv$sizeCol]))
    } else {
#      sizeclasses <- NULL # why no rv$?
      rv$sizeclasses <- NULL
      rv$sizeclasses_k <- NULL
#      rv$sizeclasses_CP <- NULL
#      rv$sizeclasses_SE <- NULL
    }
#    rv$nsizeclasses <- length(sizeclasses)
    rv$nsizeclasses <- length(rv$sizeclasses)
    rv$nsizeclasses_k <- length(rv$sizeclasses_k)
    tmp <- rv$DWPCol
    if (rv$nsizeclasses > 1 & (is.null(tmp) || !(tmp %in% rv$sizeclasses)))
      rv$DWPCol <- rv$sizeclasses[1]
    if (rv$nsizeclasses == 0){
      rv$sizeclasses_k <- ""
      rv$nsizeclasses_k <- 1
      rv$DWPCol <- NULL
    }
  }

  if (eventName == "obsSE"){
    rv$obsCols_SE <- input$obsSE
    rv$colNames_SE_preds <- removeCols(rv$colNames_SE_preds0,
      c(rv$obsCols_SE, rv$sizeCol))
    rv$colNames_size <- removeCols(rv$colNames_size0,
      c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta, rv$preds_CP))
  }

  if (eventName == "predsSE"){
    rv$sizeCol <- input$class
    rv$obsCols_SE <- input$obsSE
    rv$preds_SE <- input$predsSE
    rv$ltp <- input$ltp
    rv$fta <- input$fta
    rv$preds_CP <- input$predsCP
    rv$colNames_SE_obs <- removeCols(rv$colNames_SE_obs,
      c(rv$preds_SE, rv$sizeCol))
    rv$colNames_size <- removeCols(rv$colNames_size0,
      c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta, rv$preds_CP))
  }

  if (eventName == "ltp"){
    rv$ltp <- input$ltp
    rv$colNames_fta <- removeCols(rv$colNames_fta0,
      c(rv$preds_CP, rv$ltp, rv$sizeCol))
    rv$colNames_CP_preds <- removeCols(rv$colNames_CP_preds0,
      c(rv$ltp, rv$fta, rv$sizeCol))
    rv$colNames_size <- removeCols(rv$colNames_size0,
      c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta, rv$preds_CP))
  }

  if (eventName == "fta"){
    rv$fta <- input$fta
    rv$colNames_ltp <- removeCols(rv$colNames_ltp0,
      c(rv$preds_CP, rv$fta, rv$sizeCol))
    rv$colNames_CP_preds <- removeCols(rv$colNames_CP_preds0,
      c(rv$ltp, rv$fta, rv$sizeCol))
    rv$colNames_size <- removeCols(rv$colNames_size0,
      c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta, rv$preds_CP))
  }

  if (eventName == "predsCP"){
    rv$preds_CP <- input$predsCP
    rv$colNames_ltp <- removeCols(rv$colNames_ltp0,
      c(rv$preds_CP, rv$fta, rv$sizeCol))
    rv$colNames_fta <- removeCols(rv$colNames_fta0,
      c(rv$preds_CP, rv$ltp, rv$sizeCol))
    rv$colNames_size <- removeCols(rv$colNames_size0,
      c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta, rv$preds_CP))
  }


  if (eventName == "run_SE"){
    toNULL <- c("predictors_SE", "formula_p", "formula_k", "outSEpk",
                "mods_SE", "mods_SE_og", "sizeclasses_SE", "AICcTab_SE",
                "modOrder_SE", "modNames_SE", "modNames_SEp", "modNames_SEk",
                "modSet_SE", "best_SE", "modTab_SE",
                "modTabPretty_SE", "modTabDL_SE", "M", "Msplit", "unitCol",
                "sizeCol_M", "SEmodToUse", "split_CO", "split_SS",
                "SStemp", "avgSI", "sizeclasses_g", "nsizeclasses_g",
                "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_SE", "figW_SE", "figH_M",
                 "figW_M", "figH_g", "figW_g")
    rv <- reVal(rv, toReVal)
    rv$kFixed <- numeric(rv$nsizeclasses_k) + NA
    if (length(rv$kFixed) > 1) {
      names(rv$kFixed) <- rv$sizeclasses_k
      for (sci in rv$sizeclasses_k){
        rv$kFixed[sci] <- input[[paste0("kFixed_val_", sci)]]
      }
    } else {
      rv$kFixed <- input[["kFixed_val_"]]
    }
    rv$obsCols_SE <- input$obsSE
    rv$preds_SE <- input$predsSE
    rv$predictors_SE <- prepPredictors(rv$preds_SE)
    rv$formula_p <- formula(paste0("p~", rv$predictors_SE))
    rv$formula_k <- formula(paste0("k~", rv$predictors_SE))

    rv$CL <- input$CL
    rv$sizeCol <- input$class

    rv$mods_SE <- suppressWarnings(
                    pkmSize(formula_p = rv$formula_p,
                      formula_k = rv$formula_k, data = rv$data_SE,
                      obsCol = rv$obsCols_SE, sizeCol = rv$sizeCol,
                      kFixed = rv$kFixed, kInit = 0.7,
                      CL = rv$CL, quiet = TRUE, allCombos = TRUE)
                  )
    rv$mods_SE_og <- rv$mods_SE
    rv$mods_SE <- pkmSetSizeFailRemove(rv$mods_SE)
    if (!all(unlist(pkmSetSizeFail(rv$mods_SE))) &&
        !any(unlist(lapply(rv$mods_SE_og, pkmSetAllFail)))){
      rv$sizeclasses <- updateSizeclasses(rv$data_SE, rv$sizeCol)
      rv$sizeclasses_SE <- sort(rv$sizeclasses)
      rv$sizeclass <- pickSizeclass(rv$sizeclasses, input$outSEclass)
      rv$sizeclass_SE <- rv$sizeclass
      rv$AICcTab_SE <- aicc(rv$mods_SE[[rv$sizeclass_SE]], quiet = TRUE, app = TRUE)
      rv$modOrder_SE <- as.numeric(row.names(rv$AICcTab_SE))
      rv$modNames_SE <- names(rv$mods_SE[[rv$sizeclass_SE]])[rv$modOrder_SE]
      rv$modNames_SEp <- modNameSplit(rv$modNames_SE, 1)
      rv$modNames_SEk <- modNameSplit(rv$modNames_SE, 2)
      rv$modSet_SE <- rv$mods_SE[[rv$sizeclass_SE]]
      rv$best_SE <- (names(rv$modSet_SE)[rv$modOrder_SE])[1]
      rv$modTab_SE <- rv$mods_SE[[rv$sizeclass_SE]][[rv$best_SE]]$cell_pk
      rv$modTabPretty_SE <- prettyModTabSE(rv$modTab_SE, rv$CL)
      rv$modTabDL_SE <- dlModTabSE(rv$modTab_SE, rv$CL)
      cells_set <- modelSetCells(rv$modSet_SE)
      n_col <- length(unique(cells_set[ , 1]))
      n_row <- nrow(cells_set)/n_col
      rv$figH_SE <- .res * (.header + .box_H + .buffer +  n_row * .panel_H + .footer)
      rv$figW_SE <- min(1200, .res * .panel_W * max(2, n_col))
    }
    rv$outSEpk <- modNamePaste(c(input$outSEp, input$outSEk))
  }


  if (eventName == "run_SE_clear"){
    toNULL <- c("predictors_SE", "formula_p", "formula_k", "kFixed", "outSEpk",
                "mods_SE", "mods_SE_og", "sizeclasses_SE", "AICcTab_SE",
                "modOrder_SE", "modNames_SE", "modNames_SEp", "modNames_SEk",
                "modSet_SE", "best_SE", "modTab_SE", "modTabPretty_SE",
                "modTabDL_SE", "M", "Msplit", "unitCol", "sizeCol_M",
                "SEmodToUse", "split_CO", "split_SS", "SStemp", "avgSI",
                "sizeclasses_g", "nsizeclasses_g", "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_SE", "figW_SE", "figH_M",
                 "figW_M", "figH_g", "figW_g")
    rv <- reVal(rv, toReVal)

    rv$toRemove_sizeCol <- c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta,
                                  rv$preds_CP)
    rv$colNames_size <- removeCols(rv$colNames_size0, rv$toRemove_sizeCol)
  }

  if (eventName == "outSEclass"){
    if (length(rv$mods_SE) > 0){
      rv$sizeclass <- pickSizeclass(rv$sizeclasses, input$outSEclass)
      rv$sizeclass_SE <- rv$sizeclass
      rv$AICcTab_SE <- aicc(rv$mods_SE[[rv$sizeclass_SE]], quiet = TRUE, 
                                            app = TRUE)
      rv$modOrder_SE <- as.numeric(row.names(rv$AICcTab_SE))
      rv$modNames_SE <- names(rv$mods_SE[[rv$sizeclass_SE]])[rv$modOrder_SE]
      rv$modNames_SEp <- modNameSplit(rv$modNames_SE, 1)
      rv$modNames_SEk <- modNameSplit(rv$modNames_SE, 2)
      rv$modSet_SE <- rv$mods_SE[[rv$sizeclass_SE]]
      rv$best_SE <- (names(rv$modSet_SE)[rv$modOrder_SE])[1]
      rv$modTab_SE <- rv$mods_SE[[rv$sizeclass_SE]][[rv$best_SE]]$cell_pk
      rv$modTabPretty_SE <- prettyModTabSE(rv$modTab_SE, rv$CL)
      rv$modTabDL_SE <- dlModTabSE(rv$modTab_SE, rv$CL)
      cells_set <- modelSetCells(rv$modSet_SE)
      n_col <- length(unique(cells_set[ , 1]))
      n_row <- nrow(cells_set)/n_col
      rv$figH_SE <- .res * (.header + .box_H + .buffer +  n_row * .panel_H + .footer)
      rv$figW_SE <- min(1200, .res * .panel_W * max(2, n_col))
    }
  }

  if (eventName == "outSEp"| eventName == "outSEk"){
    if (length(rv$mods_SE) > 0){
      rv$outSEpk <- modNamePaste(c(input$outSEp, input$outSEk))
      rv$modSet_SE <- rv$mods_SE[[rv$sizeclass]]
      if (rv$outSEpk %in% names(rv$modSet_SE)){
        rv$modTab_SE <- rv$modSet_SE[[rv$outSEpk]]$cell_pk
        rv$modTabPretty_SE <- prettyModTabSE(rv$modTab_SE, rv$CL)
        rv$modTabDL_SE <- dlModTabSE(rv$modTab_SE, rv$CL)
      } else {
        rv$modTab_SE <- NULL
        holder <- data.frame(msg = "Selected model was not successfully fit.")
        rv$modTabPretty_SE <- holder
        rv$modTabDL_SE <- holder
      }
    }
  }

  if (eventName == "run_CP"){
    toNULL <- c("dist", "predictors_CP", "formula_l", "formula_s",
                "mods_CP", "mods_CP_og", "sizeclasses_CP", "AICcTab_CP",
                "modOrder_CP", "modNames_CP", "modNames_CPl", "modNames_CPs",
                "modNames_CPdist", "modSet_CP", "best_CP", "modTab_CP",
                "modTabPretty_CP", "modTabDL_CP", "M", "Msplit", "unitCol",
                "sizeCol_M", "SEmodToUse", "split_CO", "split_SS", "SStemp",
                "avgSI", "sizeclasses_g", "nsizeclasses_g", "gGeneric",
                "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_CP", "figW_CP", "figH_M",
                 "figW_M", "figH_g", "figW_g")
    rv <- reVal(rv, toReVal)
    rv$ltp <- input$ltp
    rv$fta <- input$fta
    rv$preds_CP <- input$predsCP
    rv$dist <- input$dist
    rv$nsim <- input$nsim
    rv$CL <- input$CL
    rv$sizeCol <- input$class
    rv$predictors_CP <- prepPredictors(rv$preds_CP)
    rv$formula_l <- formula(paste("l~", rv$predictors_CP, sep = ""))
    rv$formula_s <- formula(paste("s~", rv$predictors_CP, sep = ""))

    rv$mods_CP <- suppressWarnings(
                    cpmSize(formula_l = rv$formula_l,
                      formula_s = rv$formula_s, data = rv$data_CP,
                      left = rv$ltp, right = rv$fta, dist = rv$dist,
                      sizeCol = rv$sizeCol, CL = rv$CL, quiet = TRUE,
                      allCombos = TRUE
                    )
                  )
    rv$mods_CP_og <- rv$mods_CP
    rv$mods_CP <- cpmSetSizeFailRemove(rv$mods_CP)

    if (!all(unlist(cpmSetSizeFail(rv$mods_CP)))){
      rv$sizeclasses <- updateSizeclasses(rv$data_CP, rv$sizeCol)
      rv$sizeclasses_CP <- sort(rv$sizeclasses)
      rv$sizeclass <- pickSizeclass(rv$sizeclasses, input$outCPclass)
      rv$sizeclass_CP <- rv$sizeclass
      rv$AICcTab_CP <- aicc(rv$mods_CP[[rv$sizeclass_CP]],
        quiet = TRUE, app = TRUE)
      rv$AICcTab_CP[ , "Scale Formula"] <-
        gsub("NULL", "", rv$AICcTab_CP[ , "Scale Formula"])
      rv$modOrder_CP <- as.numeric(row.names(rv$AICcTab_CP))
      rv$modNames_CP <- names(rv$mods_CP[[rv$sizeclass_CP]])[rv$modOrder_CP]
      rv$modNames_CPdist <- modNameSplit(rv$modNames_CP, 1)
      rv$modNames_CPl <- modNameSplit(rv$modNames_CP, 2)
      rv$modNames_CPs <- modNameSplit(rv$modNames_CP, 3)
      rv$modSet_CP <- rv$mods_CP[[rv$sizeclass_CP]]
      rv$best_CP <- (names(rv$modSet_CP)[rv$modOrder_CP])[1]
      rv$modTab_CP <- desc(model_CP = rv$mods_CP[[rv$sizeclass_CP]][[rv$best_CP]],
        CL = rv$CL)
      #rv$best_CP <- gsub("NULL", "s ~ 1", rv$best_CP)
      # size of graph page
      modelSet <- tidyModelSetCP(rv$modSet_CP)
      cells_set <- modelSetCells(modelSet)
      preds_set <- modelSetPredictors(modelSet)
      n_col <- ifelse(length(preds_set) == 0, 1,
        length(unique(cells_set[ , preds_set[1]])))
      n_row <- nrow(cells_set)/n_col
      rv$figH_CP <- .res * (.header + n_row * .panel_H + .footer)
      rv$figW_CP <- min(1200, .res * .panel_W * 2 * n_col)

#      rv$figH_CP <- setFigH(rv$modSet_CP, "CP")
#      rv$figW_CP <- setFigW(rv$modSet_CP)
    }
  }

  if (eventName == "run_CP_clear"){
    toNULL <- c("dist", "predictors_CP", "formula_l", "formula_s",
                "mods_CP", "mods_CP_og", "sizeclasses_CP", "AICcTab_CP",
                "modOrder_CP", "modNames_CP", "modNames_CPl", "modNames_CPs",
                "modNames_CPdist", "modSet_CP", "best_CP", "modTab_CP",
                "modTabPretty_CP", "modTabDL_CP", "M", "Msplit", "unitCol",
                "sizeCol_M", "SEmodToUse", "split_CO", "split_SS",
                "SStemp", "avgSI", "sizeclasses_g", "nsizeclasses_g",
                "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_CP", "figW_CP", "figH_M",
                 "figW_M", "figH_g", "figW_g")
    rv <- reVal(rv, toReVal)

    rv$toRemove_sizeCol <- c(rv$obsCols_SE, rv$preds_SE, rv$ltp, rv$fta,
                                  rv$preds_CP)
    rv$colNames_size <- removeCols(rv$colNames_size0, rv$toRemove_sizeCol)
  }

  if (eventName == "outCPclass"){
    if (length(rv$mods_CP) > 0){
      rv$sizeclass <- pickSizeclass(rv$sizeclasses, input$outCPclass)
      rv$sizeclass_CP <- rv$sizeclass
      rv$AICcTab_CP <- aicc(rv$mods_CP[[rv$sizeclass_CP]], quiet = TRUE,
                                             app = TRUE)
      rv$modOrder_CP <- as.numeric(row.names(rv$AICcTab_CP))
      rv$modNames_CP <- names(rv$mods_CP[[rv$sizeclass_CP]])[rv$modOrder_CP]
      rv$modNames_CPdist <- modNameSplit(rv$modNames_CP, 1)
      rv$modNames_CPl <- modNameSplit(rv$modNames_CP, 2)
      rv$modNames_CPs <- modNameSplit(rv$modNames_CP, 3)
      rv$modSet_CP <- rv$mods_CP[[rv$sizeclass_CP]]
      rv$best_CP <- (names(rv$modSet_CP)[rv$modOrder_CP])[1]
      rv$modTab_CP <- desc(model_CP = rv$mods_CP[[rv$sizeclass_CP]][[rv$best_CP]],
        CL = rv$CL)
      modelSet <- tidyModelSetCP(rv$modSet_CP)
      cells_set <- modelSetCells(modelSet)
      preds_set <- modelSetPredictors(modelSet)
      n_col <- ifelse(length(preds_set) == 0, 1,
        length(unique(cells_set[ , preds_set[1]])))
      n_row <- nrow(cells_set)/n_col
      rv$figH_CP <- .res * (.header + n_row * .panel_H + .footer)
      rv$figW_CP <- min(1200, .res * .panel_W * 2* n_col)
#      rv$figH_CP <- setFigH(rv$modSet_CP, "CP")
#      rv$figW_CP <- setFigW(rv$modSet_CP)
      #rv$best_CP <- gsub("NULL", "s ~ 1", rv$best_CP)
    }
  }

  if (eventName %in% c("outCPdist", "outCPl", "outCPs")){
    if (length(rv$mods_CP) > 0){
      rv$CPdls <- c(input$outCPdist, input$outCPl, input$outCPs)
      rv$outCPdlsfig <- modNamePaste(rv$CPdls, "CP")
      rv$outCPdlstab <- modNamePaste(rv$CPdls, "CP", tab = TRUE)
      rv$modSet_CP <- rv$mods_CP[[rv$sizeclass]]
      if (rv$outCPdlstab %in% names(rv$modSet_CP)){
        rv$modTab_CP <- desc(model_CP = rv$modSet_CP[[rv$outCPdlstab]], CL = rv$CL)
      } else {
        rv$modTab_CP <- NULL
      }
    }
  }

  if (eventName == "run_g"){
    rv$CL <- input$CL
    toNULL <- c("sizeclasses_g", "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("figH_g", "figW_g")
    rv <- reVal(rv, toReVal)
    if (length(rv$obsCols_SE) == 1 & any(is.na(rv$kFixed))){
      rv$kCheck_g <- rep(NA, rv$nsizeclasses_k)
      names(rv$kCheck_g) <- rv$sizeclasses_k
      counter <- 1
      for (sci in rv$sizeclasses_k){
        rv$kCheck_g[sci] <- rv$kFixed[sci]
      }
      if (length(na.omit(rv$kCheck_g)) != length(rv$kCheck_g)){
        return(rv)  
      }
    }
     rv$sizeclasses_g <- rv$sizeclasses
    rv$nsizeclasses_g <- length(rv$sizeclasses_g)
    if (length(rv$nsizeclasses_g) == 1){
      if (is.null(rv$sizeclasses_g)){
        rv$sizeclasses_g <- "all"
        rv$nsizeclasses_g <- 1
      }
    }
    rv$nsim <- input$nsim
    rv$gGeneric <- list()
    for (sci in 1:length(rv$sizeclasses_g)){
      if (is.null(input[[paste0("modelChoices_SE", sci)]])){
        showNotification(paste0("No SE model selected for ",
          rv$sizeclasses[sci], "...error."), type = "error")
        return(rv)
      }

      if (is.null(input[[paste0("modelChoices_CP", sci)]])){
        showNotification(paste0("No CP model selected for ",
          rv$sizeclasses[sci], "...error."), type = "error")
        return(rv)
      }
      rv$SEmodToUse_g <- input[[paste0("modelChoices_SE", sci)]]
      rv$CPmodToUse_g <- input[[paste0("modelChoices_CP", sci)]]

      rv$SEmodToUse_g <- gsub("~ constant", "~ 1", rv$SEmodToUse_g)
      rv$CPmodToUse_g <- gsub("~ constant", "~ 1", rv$CPmodToUse_g)
      if (!grepl("s ~", rv$CPmodToUse_g)){
        rv$CPmodToUse_g <- paste(rv$CPmodToUse_g, "; NULL", sep = "")
      }
      rv$CPmodToUse_g <- paste("dist: ", rv$CPmodToUse_g, sep = "")
      rv$SS <- list("span" = input[["gSearchMax"]], "I" = input[["gSearchInterval"]])
      if (is.null(rv$SS[["span"]]) || is.null(rv$SS[["I"]]) ||
        is.na(rv$SS[["span"]]) || is.na(rv$SS[["span"]])){
        showNotification("Enter search interval and span before estimating g",
          type = "error")
          return(rv)
      }
      if (rv$SS[["span"]] < 2 * rv$SS[["I"]]){
        days <- c(0, rv$SS[["I"]])
      } else {
        days <- tryCatch(seq(0, rv$SS[["span"]], by  = rv$SS[["I"]]), error = function(x) NULL)
      }
      if (is.null(days)) return(rv)
      rv$gGeneric[[sci]] <- tryCatch(
                              estgGeneric(nsim = rv$nsim, days = days,
                              model_SE = rv$mods_SE[[sci]][[rv$SEmodToUse_g]],
                              model_CP = rv$mods_CP[[sci]][[rv$CPmodToUse_g]]
                              ),
                              error = function(x){NULL}
                            )
    }
    names(rv$gGeneric) <- rv$sizeclasses_g
    rv$sizeclass_g <- rv$sizeclasses_g[1]
   }

  if (eventName == "run_g_clear"){
    toNULL <- c("sizeclasses_g", "gGeneric", "SEmodToUse_g")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("figH_g", "figW_g")
    rv <- reVal(rv, toReVal)
  }

  if (eventName == "outgclass"){
    rv$sizeclass_g <- pickSizeclass(rv$sizeclasses_g, input$outgclass)
    rv$CL <- input$CL
  }

  if (eventName == "run_M"){
    toNULL <- c("M", "Msplit", "unitCol", "sizeCol_M", "SEmodToUse",
                "split_CO", "split_SS")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
    if (length(rv$obsCols_SE) == 1 & any(is.na(rv$kFixed))){
      rv$kCheck <- rep(NA, rv$nsizeclasses_k)
      names(rv$kCheck) <- rv$sizeclasses_k
      counter <- 1
      for (sci in rv$sizeclasses_k){
        rv$kCheck[sci] <- rv$kFixed[sci]
      }
      if (length(na.omit(rv$kCheck)) != length(rv$kCheck)){
        return(rv)
      }
    }
    rv$nsizeclasses <- length(rv$sizeclasses)
    if (length(rv$nsizeclasses) == 1){
      if (is.null(rv$sizeclasses)){
        rv$sizeclasses <- "all"
       }
    }
    rv$COdate <- input$COdate
    rv$nsim <- input$nsim
    rv$xIDcol <- input$xID
    rv$frac <- input$frac
    if (rv$frac < 0.01 | rv$frac > 1) return(rv)
    rv$SEmodToUse <- rep(NA, rv$nsizeclasses)
    rv$CPmodToUse <- rep(NA, rv$nsizeclasses)
    names(rv$SEmodToUse) <- rv$sizeclasses
    names(rv$CPmodToUse) <- rv$sizeclasses

    for (sci in 1:length(rv$sizeclasses)){
      if (is.null(input[[paste0("modelChoices_SE", sci)]])){
        showNotification(paste0("No SE model selected for ",
          rv$sizeclasses[sci], "...error."), type = "error")
        return(rv)
      }
      rv$SEmodToUse[sci] <- input[[paste0("modelChoices_SE", sci)]]
      if (is.null(input[[paste0("modelChoices_CP", sci)]])){
        showNotification(paste0("No CP model selected for ",
          rv$sizeclasses[sci], "...error."), type = "error")
        return(rv)
      }
      rv$CPmodToUse[sci] <- input[[paste0("modelChoices_CP", sci)]]
      if (!grepl("s ~", rv$CPmodToUse[rv$sizeclasses[sci]])){
        rv$CPmodToUse[sci] <- paste(rv$CPmodToUse[rv$sizeclasses[sci]], "; NULL", sep = "")
      }
      rv$CPmodToUse[sci] <- paste("dist: ", rv$CPmodToUse[rv$sizeclasses[sci]], sep = "")
    }

    rv$SEmodToUse <- gsub("~ constant", "~ 1", rv$SEmodToUse)
    rv$CPmodToUse <- gsub("~ constant", "~ 1", rv$CPmodToUse)

    rv$models_SE <- tryCatch(
                      trimSetSize(rv$mods_SE, rv$SEmodToUse),
                      error = function(x){NULL}
                    )
    rv$models_CP <- tryCatch(
                      trimSetSize(rv$mods_CP, rv$CPmodToUse),
                      error = function(x){NULL}
                    )
    if(any(c(is.null(rv$models_SE), is.null(rv$models_CP)))){
      rv$M <- NULL
      return(rv)
    }
    if (rv$nsizeclasses > 1){
      dwpcol <- NULL
      rv$sizeCol_M <- rv$sizeCol
    } else {
      rv$DWPCol <- input$DWPCol
      dwpcol <- rv$DWPCol
      rv$sizeCol_M <- NULL
      rv$models_SE <- rv$models_SE[[1]]
      rv$models_CP <- rv$models_CP[[1]]
    }
    rv$CL <- input$CL
    rv$M <- tryCatch(
              estM(data_CO = rv$data_CO, data_SS = rv$data_SS, data_DWP = rv$data_DWP,
                frac = rv$frac, model_SE = rv$models_SE,
                model_CP = rv$models_CP,
                COdate = rv$COdate, DWPCol = dwpcol,
                sizeCol = rv$sizeCol_M, nsim = rv$nsim, IDcol = rv$xIDcol,
                max_intervals = 8
              ), error = function(e) e)
    if (!("error" %in% class(rv$M))){
      rv$Msplit <- tryCatch(
                     calcSplits(M = rv$M,
                       split_SS = NULL, split_CO = NULL,
                       data_SS = rv$data_SS, data_CO = rv$data_CO
                     ), error = function(x){NULL}, warning = function(x){NULL}
                   )
      rv$unitCol <- intersect(rv$colNames_CO, rv$colNames_DWP)
      rv$splittable_SS <- rv$colNames_SS
      badind <- NULL
      for (ci in 1:length(rv$splittable_SS)){
        if (is.numeric(rv$data_SS[ , rv$splittable_SS[ci]])){
          badind <- c(badind, ci)
        } else {
          SSlev <- rv$data_SS[ , rv$splittable_SS[ci]]
          if (min(diff(match(SSlev, unique(SSlev)))) < 0)
            badind <- c(badind, ci)
        }
      }
     if (length(badind) > 0) rv$splittable_SS <-   rv$splittable_SS[-badind]
    }
  }

  if (eventName == "run_M_clear"){
    toNULL <- c("M", "Msplit", "unitCol", "sizeCol_M", "SEmodToUse", 
                "split_CO", "split_SS")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
  }

  if (eventName == "split_M"){
    toNULL <- c("split_CO", "split_SS", "Msplit")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv$split_CO <- input$split_CO
    rv$split_SS <- input$split_SS
    rv$nsplit_CO <- length(rv$split_CO)
    rv$nsplit_SS <- length(rv$split_SS)
    rv$COdate <- input$COdate

    rv$Msplit <- tryCatch(
                   calcSplits(M = rv$M,
                     split_SS = rv$split_SS, split_CO = rv$split_CO,
                     data_SS = rv$data_SS, data_CO = rv$data_CO
                   ), error = function(x){NULL}, warning = function(x){NULL}
                 )
    if (!is.null(rv$Msplit)){
      rv$figH_M <- 600
      if (length(attr(rv$Msplit, "vars")) > 1){
        rv$figH_M <- max(600, 300 * length(rv$Msplit))
      }
    }
  }

  if (eventName == "split_M_clear"){
    toNULL <- c("split_CO", "split_SS", "Msplit")
    rv <- reNULL(rv, toNULL)
    toReVal <- c("nsplit_CO", "nsplit_SS", "figH_M", "figW_M")
    rv <- reVal(rv, toReVal)
    if (!is.null(rv$M)){
      rv$Msplit <- tryCatch(
                     calcSplits(M = rv$M,
                       split_SS = NULL, split_CO = NULL,
                       data_SS = rv$data_SS, data_CO = rv$data_CO
                     ), error = function(x){NULL}, warning = function(x){NULL}
                   )
      rv$unitCol <- intersect(rv$colNames_CO, rv$colNames_DWP)
      rv$splittable_SS <- rv$colNames_SS
      badind <- NULL
      for (ci in 1:length(rv$splittable_SS)){
        if (is.numeric(rv$data_SS[ , rv$splittable_SS[ci]])){
          badind <- c(badind, ci)
        } else {
          SSlev <- rv$data_SS[ , rv$splittable_SS[ci]]
          if (min(diff(match(SSlev, unique(SSlev)))) < 0)
            badind <- c(badind, ci)
        }
      }
     if (length(badind) > 0) rv$splittable_SS <-   rv$splittable_SS[-badind]
   }
  }

  if (eventName == "transpose_split"){
    if (rv$nsplit_CO + rv$nsplit_SS == 2){
      rv$Msplit <- transposeSplits(rv$Msplit)
    }
  }
  return(rv)
}