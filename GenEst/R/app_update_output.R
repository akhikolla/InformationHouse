#' @title Update the outputs when an event occurs
#'
#' @description When an event occurs in the GenEst GUI, the output values may
#'   need to be updated. This function contains all of the possible updates
#'   based on the event options (or lacks any updates if the event doesn't
#'   require any).
#'
#' @param eventName Character name of the event. One of "clear_all",
#'   "file_SE", "file_SE_clear", "file_CP", "file_CP_clear", "file_SS",
#'   "file_SS_clear", "file_DWP", "file_DWP_clear", "file_CO", 
#'   "file_CO_clear", "class", "obsSE", "predsSE", "run_SE", "run_SE_clear",
#'   "outSEclass", "outSEp", "outSEk", "ltp", "fta", "predsCP", "run_CP",
#'   "run_CP_clear", "outCPclass", "outCPdist", "outCPl", "outCPs",
#'   "run_M", "run_M_clear", "split_M", "split_M_clear", "transpose_split",
#'   "run_g", "run_g_clear", or "outgclass".
#'
#' @param rv Reactive values list for the GenEst GUI.
#'
#' @param output \code{output} list for the GenEst GUI.
#'
#' @param input \code{input} lisst for the GenEst GUI
#'
#' @return Updated \code{output} list.
#'
update_output <- function(eventName, rv, output, input){

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

  if (eventName == "clear_all" | grepl("load_", eventName)){
    toNULL <- c("data_SE", "filename_SE", "selected_SE", "SEModDone",
                "DWPNeed", "AICcTab_SE", "modTab_SE", "fig_SE",
                "sizeclasses_SE", "modelMenu_SE", "sizeclass_SE1",
                "sizeclass_SE2", "sizeclass_SE3", "sizeclass_SEyn",
                "text_SE_est", "data_CP", "filename_CP", "selected_CP",
                "CPModDone", "AICcTab_CP", "modTab_CP", "fig_CP",
                "sizeclasses_CP", "modelMenu_CP", "sizeclass_CP1",
                "sizeclass_CP2", "sizeclass_CP3", "sizeclass_CPyn",
                "text_CP_est", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2", "data_SS", "data_DWP", "data_CO", "filename_SS",
                "filename_CO", "filename_DWP")
    output <- reNULL(output, toNULL)
    output$kNeed <- setkNeed(rv)
    dontSuspend <- c("SEModDone", "sizeclasses_SE", "text_SE_est",
                     "kNeed", "MModDone", "gModDone", "sizeclass_gyn",
                     "CPModDone", "sizeclasses_CP", "text_CP_est",
      "filename_SE", "filename_CP", "filename_SS", "filename_DWP", "filename_CO")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "file_SE" | eventName == "file_SE_clear"){
    toNULL <- c("data_SE", "filename_SE", "selected_SE", "SEModDone",
                "DWPNeed", "AICcTab_SE", "modTab_SE", "fig_SE",
                "sizeclasses_SE", "modelMenu_SE", "sizeclass_SE1",
                "sizeclass_SE2", "sizeclass_SE3", "sizeclass_SEyn",
                "text_SE_est", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")
    output <- reNULL(output, toNULL)
    output$kNeed <- setkNeed(rv)
    if (eventName == "file_SE"){
      output$data_SE <- DT::renderDataTable(datatable(rv$data_SE,
        caption = paste0("File: ", rv$filename_SE)), server = FALSE)
      output$filename_SE <- renderText(paste0("File: ", rv$filename_SE))
    }
    dontSuspend <- c("SEModDone", "sizeclasses_SE", "filename_SE",
                     "text_SE_est", "kNeed", "MModDone", "gModDone",
                     "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "file_CP" | eventName == "file_CP_clear"){
    toNULL <- c("data_CP", "filename_CP", "selected_CP", "CPModDone",
                "AICcTab_CP", "modTab_CP", "fig_CP",
                "sizeclasses_CP", "modelMenu_CP", "sizeclass_CP1",
                "sizeclass_CP2", "sizeclass_CP3", "sizeclass_CPyn",
                "text_CP_est", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")
    output <- reNULL(output, toNULL)
    if (eventName == "file_CP"){
      output$data_CP <- DT::renderDataTable(datatable(rv$data_CP,
          caption = paste0("File: ", rv$filename_CP)), server = FALSE)
      output$filename_CP <- renderText(paste0("File: ", rv$filename_CP))
    }
    dontSuspend <- c("CPModDone", "sizeclasses_CP", "filename_CP",
                     "text_CP_est", "MModDone", "gModDone",
                     "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "file_SS" | eventName == "file_SS_clear"){
    toNULL <- c("data_SS", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2", "filename_SS")
    output <- reNULL(output, toNULL)
    if (eventName == "file_SS"){
      output$data_SS <- DT::renderDataTable(datatable(rv$data_SS, caption = paste0(
        "File: ", rv$filename_SS)), server = FALSE)
      output$filename_SS <- renderText(paste0("File: ", rv$filename_SS))
    }
    dontSuspend <- c("filename_SS", "MModDone", "gModDone", "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "file_DWP" | eventName == "file_DWP_clear"){
    toNULL <- c("filename_DWP", "data_DWP", "fig_M", "table_M", "MModDone")
    output <- reNULL(output, toNULL)
    if (eventName == "file_DWP"){
      output$data_DWP <- DT::renderDataTable(datatable(rv$data_DWP,
        caption = paste0("File: ", rv$filename_DWP)))
      output$filename_DWP <- renderText(paste0("File: ", rv$filename_DWP))
    }
    dontSuspend <- c("filename_DWP", "MModDone")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "file_CO" | eventName == "file_CO_clear"){
    toNULL <- c("data_CO", "filename_CO", "fig_M", "table_M", "MModDone")
    output <- reNULL(output, toNULL)
    if (eventName == "file_CO"){
      output$data_CO <- DT::renderDataTable(datatable(rv$data_CO,
        caption = paste0("File: ", rv$filename_CO)), server = FALSE)
      output$filename_CO <- renderText(paste0("File: ", rv$filename_CO))
    }
    dontSuspend <- c("filename_CO", "MModDone")
    setNotSuspending(output, dontSuspend)
  }


  if (grepl("load_", eventName)){
    toNULL <- c("data_SE", "data_CP", "data_SS", "data_DWP", "data_CO",
                "filename_SE", "filename_CP", "filename_SS", "filename_CO",
                "filename_DWP", "selected_SE", "SEModDone",
                "DWPNeed", "AICcTab_SE", "modTab_SE", "fig_SE",
                "sizeclasses_SE", "modelMenu_SE", "sizeclass_SE1",
                "sizeclass_SE2", "sizeclass_SE3", "sizeclass_SEyn",
                "text_SE_est", "selected_CP",
                "CPModDone", "AICcTab_CP", "modTab_CP", "fig_CP",
                "sizeclasses_CP", "modelMenu_CP", "sizeclass_CP1",
                "sizeclass_CP2", "sizeclass_CP3", "sizeclass_CPyn",
                "text_CP_est", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")

    output <- reNULL(output, toNULL)
    output$kNeed <- setkNeed(rv)

    output$data_SE <- DT::renderDataTable({datatable(rv$data_SE,
      caption = paste0("File: ", rv$filename_SE))}, server = FALSE)
    output$filename_SE <- renderText(paste0("File: ", rv$filename_SE))

    output$data_CP <- DT::renderDataTable({datatable(rv$data_CP,
      caption = paste0("File: ", rv$filename_CP))}, server = FALSE)
    output$filename_CP <- renderText(paste0("File: ", rv$filename_CP))

    output$data_SS <- DT::renderDataTable({datatable(rv$data_SS,
      caption = paste0("File: ", rv$filename_SS))}, server = FALSE)
    output$filename_SS <- renderText(paste0("File: ", rv$filename_SS))

    output$data_DWP <- DT::renderDataTable({datatable(rv$data_DWP,
      caption = paste0("File: ", rv$filename_DWP))}, server = FALSE)
    output$filename_DWP <- renderText(paste0("File: ", rv$filename_DWP))

    output$data_CO <- DT::renderDataTable({datatable(rv$data_CO,
      caption = paste0("File: ", rv$filename_CO))}, server = FALSE)
    output$filename_CO <- renderText(paste0("File: ", rv$filename_CO))

    dontSuspend <- c("SEModDone", "sizeclasses_SE", "text_SE_est",
                     "kNeed", "MModDone", "gModDone", "sizeclass_gyn",
                     "CPModDone", "sizeclasses_CP", "text_CP_est",
                     "filename_SE", "filename_CP", "filename_SS",
                     "filename_DWP", "filename_CO")
    setNotSuspending(output, dontSuspend)

  }

  if (eventName == "class"){
    if (!is.null(rv$obsCols_SE)){
      selectedCols <- c(rv$obsCols_SE, rv$sizeCol, rv$preds_SE)
      selectedData <- selectData(rv$data_SE, selectedCols)
      output$selected_SE <- DT::renderDataTable(datatable(selectedData), server = FALSE)
    }
    if (!is.null(c(rv$ltp, rv$fta))){
      obsColsSelected <- c(rv$ltp, rv$fta)
      selectedCols <- c(obsColsSelected, rv$sizeCol, rv$preds_CP)
      selectedData <- selectData(rv$data_CP, selectedCols)
      output$selected_CP <- DT::renderDataTable(datatable(selectedData), server = FALSE)
    }
    output$sizeclasses_SE <- renderText(paste(rv$sizeclasses_SE, collapse = " "))
    output$sizeclasses_CP <- renderText(paste(rv$sizeclasses_CP, collapse = " "))
    output$kFixedInput <- kFixedWidget(rv$sizeclasses_k)
  }

  if (eventName == "obsSE"){
    selectedCols <- c(rv$obsCols_SE, rv$sizeCol, rv$preds_SE)
    if (!is.null(rv$data_SE)){
      selectedData <- selectData(rv$data_SE, selectedCols)
      output$selected_SE <- DT::renderDataTable(datatable(selectedData), server = FALSE)
    }
  }

  if (eventName == "predsSE"){
    selectedCols <- c(rv$obsCols_SE, rv$sizeCol, rv$preds_SE)
    if (!is.null(rv$data_SE)){
      selectedData <- selectData(rv$data_SE, selectedCols)
      output$selected_SE <- DT::renderDataTable(datatable(selectedData), server = FALSE)
    }
  }

  if (eventName == "fta"){
    selectedCols <- c(rv$ltp, rv$fta, rv$sizeCol, rv$preds_CP)
    selectedData <- selectData(rv$data_CP, selectedCols)
    output$selected_CP <- DT::renderDataTable(datatable(selectedData), server = FALSE)
  }

  if (eventName == "ltp"){
    selectedCols <- c(rv$ltp, rv$fta, rv$sizeCol, rv$preds_CP)
    selectedData <- selectData(rv$data_CP, selectedCols)
    output$selected_CP <- DT::renderDataTable(datatable(selectedData), server = FALSE)
  }

  if (eventName == "predsCP"){
    selectedCols <- c(rv$ltp, rv$fta, rv$sizeCol, rv$preds_CP)
    selectedData <- selectData(rv$data_CP, selectedCols)
    output$selected_CP <- DT::renderDataTable(datatable(selectedData), server = FALSE)
  }

  if (eventName == "run_SE"){
    output$kNeed <- setkNeed(rv)
    output$DWPNeed <- renderText("yes")
    toNULL <- c("SEModDone", "AICcTab_SE", "modTab_SE", "fig_SE",
                "sizeclasses_SE", "modelMenu_SE", "sizeclass_SE1",
                "sizeclass_SE2", "sizeclass_SE3", "sizeclass_SEyn",
                "text_SE_est", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")
    output <- reNULL(output, toNULL)
    if (!all(unlist(pkmSetSizeFail(rv$mods_SE))) &&
        !any(unlist(lapply(rv$mods_SE_og, pkmSetAllFail)))){
      output$SEModDone <- renderText("OK")
      output$kNeed <- setkNeed(rv)
      if (length(rv$sizeclasses) == 1){
        output$DWPNeed <- renderText("yes")
      } else{
        output$DWPNeed <- renderText("no")
      }
      output$AICcTab_SE <- renderDataTable({rv$AICcTab_SE})
      output$modTab_SE <- renderDataTable({rv$modTabPretty_SE})
      output$fig_SE <- renderPlot(
        plot(rv$modSet_SE, specificModel = rv$best_SE, CL = rv$CL),
          height = rv$figH_SE, width = rv$figW_SE,
          pointsize = .pointsize, res = .res)

      output$sizeclasses_SE <- renderText(paste(rv$sizeclasses_SE, collapse = " "))
      output$modelMenu_SE <- modelSelectionWidget(rv$mods_SE, "SE")

      output$sizeclass_SE1 <- classText(rv, "SE")
      output$sizeclass_SE2 <- classText(rv, "SE")
      output$sizeclass_SE3 <- classText(rv, "SE")
      output$text_SE_est <- estText(rv, "SE")

      if (length(rv$sizeclasses_SE) == 1){
        output$sizeclass_SEyn <- renderText("NO")
      } else{
        output$sizeclass_SEyn <- renderText("YES")
      }
      output$dlSEest <- downloadTable("SE_estimates.csv", rv$modTabDL_SE,
                                            rv$csvformat)
      output$dlSEAICc <- downloadTable("SE_AICc.csv", rv$AICcTab_SE,
                                            rv$csvformat)
      output$dlSEfig <- downloadSEFig(rv)
      output$dlSEmod <- downloadSEmod(rv, input)
    }
    dontSuspend <- c("text_SE_est", "MModDone", "gModDone", "sizeclass_gyn",
                     "SEModDone", "kNeed", "DWPNeed", "sizeclasses_SE",
                     "sizeclass_SEyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "run_SE_clear"){
    output$kNeed <- setkNeed(rv)
    output$DWPNeed <- renderText("yes")
    toNULL <- c("SEModDone", "AICcTab_SE", "modTab_SE", "fig_SE",
                "sizeclasses_SE", "modelMenu_SE", "sizeclass_SE1",
                "sizeclass_SE2", "sizeclass_SE3", "sizeclass_SEyn",
                "text_SE_est", "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")
    output <- reNULL(output, toNULL)
    dontSuspend <- c("SEModDone", "sizeclasses_SE", "text_SE_est",
                     "kNeed", "MModDone", "gModDone",
                     "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "outSEclass"){
    if (length(rv$mods_SE) > 0){
      output$AICcTab_SE <- renderDataTable({rv$AICcTab_SE})
      output$modTab_SE <- renderDataTable({rv$modTabPretty_SE})
      output$fig_SE <- renderPlot(
        plot(rv$modSet_SE, specificModel = rv$best_SE, CL = rv$CL),
        height = rv$figH_SE, width = rv$figW_SE,
        pointsize = .pointsize, res = .res)

      output$sizeclass_SE1 <- classText(rv, "SE")
      output$sizeclass_SE2 <- classText(rv, "SE")
      output$sizeclass_SE3 <- classText(rv, "SE")

      output$dlSEest <- downloadTable("SE_estimates.csv", rv$modTabDL_SE,
                                            rv$csvformat)
      output$dlSEAICc <- downloadTable("SE_AICc.csv", rv$AICcTab_SE,
                                            rv$csvformat)
      output$dlSEfig <- downloadSEFig(rv)
      output$dlSEmod <- downloadSEmod(rv, input)
    }
  }

  if (eventName == "outSEp" | eventName == "outSEk"){
    if (length(rv$mods_SE) > 0){
        output$fig_SE <- renderPlot({
          tryCatch(
            plot(rv$modSet_SE, specificModel = rv$outSEpk, CL = rv$CL),
            error = function(x){plotNA()}
          )
        },
        height = rv$figH_SE,
        width = rv$figW_SE,
        res = .res,
        pointsize = .pointsize
      )
      output$dlSEfig <- downloadSEFig(rv)
      output$dlSEmod <- downloadSEmod(rv, input)
      if (!is.null(rv$modTab_SE)){
        output$modTab_SE <- renderDataTable({rv$modTabPretty_SE})
        output$dlSEest <- downloadTable("SE_estimates.csv", rv$modTabDL_SE,
                                            rv$csvformat)
      }
    }
  }

  if (eventName == "run_CP"){
    toNULL <- c("CPModDone", "AICcTab_CP", "modTab_CP", "fig_CP",
                "sizeclasses_CP", "modelMenu_CP", "sizeclass_CP1",
                "sizeclass_CP2", "sizeclass_CP3", "sizeclass_CPyn",
                "text_CP_est", "dlCPest", "dlCPAICc", "dlCPfig",
                "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")
    output <- reNULL(output, toNULL)
    if (!all(unlist(cpmSetSizeFail(rv$mods_CP)))){
      output$CPModDone <- renderText("OK")
      output$AICcTab_CP <- renderDataTable({rv$AICcTab_CP})
      output$modTab_CP <- renderDataTable({prettyModTabCP(rv$modTab_CP)})
      output$fig_CP <- renderPlot(
        plot(rv$modSet_CP, specificModel = rv$best_CP, CL = rv$CL),
        height = rv$figH_CP, width = rv$figW_CP,
        pointsize = .pointsize, res = .res
      )

      output$sizeclasses_CP <- renderText(paste(rv$sizeclasses_CP, collapse = " "))
      output$modelMenu_CP <- modelSelectionWidget(rv$mods_CP, "CP")

      output$text_CP_est <- estText(rv, "CP")
      output$sizeclass_CP1 <- classText(rv, "CP")
      output$sizeclass_CP2 <- classText(rv, "CP")
      output$sizeclass_CP3 <- classText(rv, "CP")

      if (length(rv$sizeclasses_CP) == 1){
        output$sizeclass_CPyn <- renderText("NO")
      } else{
        output$sizeclass_CPyn <- renderText("YES")
      }
      output$dlCPest <- downloadTable("CP_estimates.csv", rv$modTab_CP,
                                            rv$csvformat)
      output$dlCPAICc <- downloadTable("CP_AICc.csv", rv$AICcTab_CP, rv$csvformat)
      output$dlCPfig <- downloadCPFig(rv)
      output$dlCPmod <- downloadCPmod(rv, input)
    }
    dontSuspend <- c("CPModDone", "sizeclasses_CP", "sizeclass_CPyn",
                     "text_CP_est", "MModDone", "gModDone", "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "run_CP_clear"){
    toNULL <- c("CPModDone", "AICcTab_CP", "modTab_CP", "fig_CP",
                "sizeclasses_CP", "modelMenu_CP", "sizeclass_CP1",
                "sizeclass_CP2", "sizeclass_CP3", "sizeclass_CPyn",
                "text_CP_est", "dlCPest", "dlCPAICc", "dlCPfig",
                "fig_M", "table_M", "MModDone", "table_g",
                "fig_g", "gModDone", "sizeclass_gyn", "sizeclass_g1",
                "sizeclass_g2")
    output <- reNULL(output, toNULL)
    dontSuspend <- c("CPModDone", "sizeclasses_CP", "sizeclass_CPyn",
                     "text_CP_est", "MModDone", "gModDone", "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "outCPclass"){
    if (length(rv$mods_CP) > 0){
      output$modTab_CP <- DT::renderDataTable(datatable(prettyModTabCP(rv$modTab_CP)))
      output$fig_CP <- renderPlot(
        plot(rv$modSet_CP, specificModel = rv$best_CP, CL = rv$CL),
        height = rv$figH_CP, width = rv$figW_CP,
        pointsize = .pointsize, res = .res
      )

      preText <- paste0("Carcass class: ", rv$sizeclass_CP)
      output$sizeclass_CP1 <- classText(rv, "CP")
      output$sizeclass_CP2 <- classText(rv, "CP")
      output$sizeclass_CP3 <- classText(rv, "CP")
      output$dlCPest <- downloadTable("CP_estimates.csv", rv$modTab_CP,
                                            rv$csvformat)
      output$dlCPAICc <- downloadTable("CP_AICc.csv", rv$AICcTab_CP,
                                            rv$csvformat)
      output$dlCPfig <- downloadCPFig(rv)
      output$dlCPmod <- downloadCPmod(rv, input)
    }
  }

  if (eventName %in% c("outCPdist", "outCPl", "outCPs")){
    if (length(rv$mods_CP) > 0){
      output$modTab_CP <- DT::renderDataTable(datatable(prettyModTabCP(rv$modTab_CP)))
      specMod <- ifelse(!grepl("exponential", rv$outCPdlsfig), rv$outCPdlsfig,
        sub("s ~ 1", "NULL", rv$outCPdlsfig))
      output$fig_CP <- renderPlot(tryCatch(
        plot(rv$modSet_CP, specificModel = specMod, CL = rv$CL),
          error = function(x) plotNA()),
        height = rv$figH_CP, width = rv$figW_CP,
        pointsize = .pointsize, res = .res)
      output$dlCPest <- downloadTable("CP_estimates.csv", rv$modTab_CP,
                                             rv$csvformat)
      output$dlCPfig <- downloadCPFig(rv)
      output$dlCPmod <- downloadCPmod(rv, input)
      if (!is.null(rv$modTab_CP)){
        output$modTab_CP <- renderDataTable({prettyModTabCP(rv$modTab_CP)})
        output$dlCPest <- downloadTable("CP_estimates.csv", rv$modTab_CP,
                                             rv$csvformat)
      }
    }
  }


  if (eventName == "run_g"){
    toNULL <- c("table_g", "fig_g", "gModDone", "sizeclass_gyn",
                "sizeclass_g1", "sizeclass_g2")
    output <- reNULL(output, toNULL)
    if (length(rv$gGeneric) > 0 && !is.null(rv$gGeneric[[1]])){
      summaryTab <- summary(rv$gGeneric[[1]], CL = rv$CL)
      output$table_g <- renderDataTable(summaryTab,
        caption = paste0("I = ", rv$SS[["I"]], ", span = ", rv$SS[["span"]]))
      output$fig_g <- renderPlot({
                        tryCatch(
                          plot(rv$gGeneric[[1]], CL = rv$CL),
                          error = function(x){plot(1,1)},
                          warning = function(x){plot(1,1)}
                        )
                      }, height = rv$figH_g, width = rv$figW_g,
                      pointsize = .pointsize, res = .res)
      output$gModDone <- renderText("OK")
      if (length(rv$sizeclasses_SE) == 1){
        output$sizeclass_gyn <- renderText("NO")
      } else{
        output$sizeclass_gyn <- renderText("YES")
      }

      output$sizeclass_g1 <- classText(rv, "g")
      output$sizeclass_g2 <- classText(rv, "g")

#      output$dlgtab <- downloadTable("g_estimates.csv", summaryTab, rv$csvformat)
      output$dlgfig <- downloadgFig(rv, 1)
      output$dlgtab <- downloadgres(rv, input)
    }
    dontSuspend <- c("gModDone", "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "run_g_clear"){
    toNULL <- c("table_g", "fig_g", "gModDone", "sizeclass_gyn",
                "sizeclass_g1", "sizeclass_g2")
    output <- reNULL(output, toNULL)
    dontSuspend <- c("gModDone", "sizeclass_gyn")
    setNotSuspending(output, dontSuspend)
  }

  if (eventName == "outgclass"){
    if (class(rv$gGeneric[[rv$sizeclass_g]])[1] == "gGeneric"){
      summaryTab <- summary(rv$gGeneric[[rv$sizeclass_g]], CL = rv$CL)
      output$table_g <- renderDataTable(summaryTab)
      output$fig_g <- renderPlot({
                        tryCatch(
                          plot(rv$gGeneric[[rv$sizeclass_g]], CL = rv$CL),
                          error = function(x){plot(1,1)},
                          warning = function(x){plot(1,1)}
                        )
                      }, height = rv$figH_g, width = rv$figW_g,
                      pointsize = .pointsize, res = .res)
      output$gModDone <- renderText("OK")
      outputOptions(output, "gModDone", suspendWhenHidden = FALSE)

      preText <- paste0(
        "Carcass class: ", rv$sizeclass_g, " ........ ",
        "Search Schedule: I = ", round(rv$SS$I, 1), ", span = ", rv$SS$span)
      if (length(rv$sizeclasses_g) == 1){
        preText <- paste0("Search Schedule: I = ", round(rv$SS$I, 1), ", span = ", rv$SS$span)
      }
      scText <- renderText(preText)
      output$sizeclass_g1 <- scText
      output$sizeclass_g2 <- scText

#      output$dlgtab <- downloadTable("g_estimates.csv", summaryTab, rv$csvformat)
      output$dlgfig <- downloadgFig(rv, rv$sizeclass_g)
      output$dlgtab <- downloadgres(rv, input)
    }
  }

  if (eventName == "run_M"){
    toNULL <- c("fig_M", "table_M", "MModDone")
    output <- reNULL(output, toNULL)
    if (!is.null(rv$Msplit)){
      output$MModDone <- renderText("OK")
      output$fig_M <- renderPlot({
        plot(rv$Msplit, CL = rv$CL)}, height = rv$figH_M, width = rv$figW_M,
        pointsize = .pointsize, res = .res)
      summaryTab <-  prettySplitTab(summary(rv$Msplit, CL = rv$CL))
      output$table_M <- renderDataTable(datatable(summaryTab))
      #output$dlMtab <- downloadTable("M_table.csv", summaryTab, rv$csvformat)
      output$dlMfig <- downloadMFig(rv)
      output$dlMres <- downloadMres(rv, input)
    }
    outputOptions(output, "MModDone", suspendWhenHidden = FALSE)
  }

  if (eventName == "run_M_clear"){
    toNULL <- c("fig_M", "table_M", "MModDone")
    output <- reNULL(output, toNULL)
    outputOptions(output, "MModDone", suspendWhenHidden = FALSE)
  }

  if (eventName == "split_M"){
    if (is.null(rv$Msplit)){
      output$fig_M <- renderPlot({
                        tryCatch(plot(rv$M, CL = rv$CL),
                          error = function(x){plotNA()}
                        )
                      }, height = rv$figH_M, width = rv$figW_M,
                      pointsize = .pointsize, res = .res)
      output$dlMfig <- downloadMFig(rv, split = FALSE)

    } else {
      output$fig_M <- renderPlot({
        tryCatch(plot(rv$Msplit, CL = rv$CL, commonScale = input$cscale == "Yes"),
          error = function(x){plotNA("split")}
        )
      }, height = rv$figH_M, width = rv$figW_M,
      pointsize = .pointsize, res = .res)
      tmp <-  prettySplitTab(summary(rv$Msplit, CL = rv$CL))
      summaryTab <- data.frame(tmp, stringsAsFactors = FALSE)
      names(summaryTab) <- colnames(tmp)
      for (sti in 1:dim(tmp)[2]) {
        testcol <- suppressWarnings(as.numeric(summaryTab[ , sti]))
        if (!anyNA(testcol)) summaryTab[, sti] <- testcol
      }
      output$table_M <- renderDataTable(datatable(summaryTab))
      #output$dlMtab <- downloadTable("M_table.csv", summaryTab, rv$csvformat)
      output$dlMfig <- downloadMFig(rv)
      output$dlMres <- downloadMres(rv, input)
    }
    output$MSplitDone <- renderText("OK")
    output$nMSplits <- renderText(
                         as.character(rv$nsplit_CO + rv$nsplit_SS))
    outputOptions(output, "nMSplits", suspendWhenHidden = FALSE)
    outputOptions(output, "MSplitDone", suspendWhenHidden = FALSE)
  }
  if (eventName == "cscale"){
    commonScale <- input$cscale == "Yes"
    output$fig_M <- renderPlot({
      plot(rv$Msplit, CL = rv$CL, commonScale = commonScale)},
        height = rv$figH_M, width = rv$figW_M,
        pointsize = .pointsize, res = .res)
  }
  if (eventName == "split_M_clear"){
    if (!is.null(rv$Msplit)){
      output$MModDone <- renderText("OK")
      outputOptions(output, "MModDone", suspendWhenHidden = FALSE)

      output$fig_M <- renderPlot({plot(rv$Msplit, CL = rv$CL)},
                        height = rv$figH_M, width = rv$figW_M,
                        pointsize = .pointsize, res = .res
                      )
      summaryTab <-  prettySplitTab(summary(rv$Msplit, CL = rv$CL))
      output$table_M <- renderDataTable(datatable(summaryTab))
      #output$dlMtab <- downloadTable("M_table.csv", summaryTab, rv$csvformat)
      output$dlMfig <- downloadMFig(rv)
      output$dlMres <- downloadMres(rv, input)
    }

    output$MSplitDone <- NULL
    output$nMSplits <- renderText(as.character(0))
    outputOptions(output, "MSplitDone", suspendWhenHidden = FALSE)
    outputOptions(output, "nMSplits", suspendWhenHidden = FALSE)
  }

  if (eventName == "transpose_split"){
    if (!is.null(rv$Msplit)){
      output$fig_M <- renderPlot({
        tryCatch(
          plot(rv$Msplit, CL = rv$CL, commonScale = input$cscale == "Yes"),
          error = function(x){plotNA("split")}
        )
      },
      height = rv$figH_M, width = rv$figW_M, pointsize = .pointsize, res = .res)
      output$dlMfig <- downloadMFig(rv, TRUE)
    }
  }

  return(output)
}