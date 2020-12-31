#' @title Update the inputs when an event occurs
#'
#' @description When an event occurs in the GenEst GUI, the input values may
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
#' @param input \code{input} list for the GenEst GUI.
#'
#' @param session Environment for the GenEst GUI.
#'
update_input <- function(eventName, rv, input, session){

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

  if(eventName == "clear_all" | grepl("load_", eventName)){
    toReset <- c("file_SE", "predsSE", "obsSE", "outSEp", "outSEk",
                 "outSEclass", "DWPCol", "split_SS", "split_CO",
                 "modelChoices_SE1", "outgclass","file_CP", "predsCP", "ltp",
                 "fta", "outCPl", "outCPs", "outCPdist", "outCPclass", 
                 "modelChoices_CP1", "file_SS", "file_DWP", "file_CO", "COdate",
                 "sizeCol")
    lapply(toReset, reset)
    updateSelectizeInput(session, "predsSE", choices = "")
    updateSelectizeInput(session, "obsSE", choices = "")
#    updateSelectizeInput(session, "class", choices = scc, selected = scs)
# NOTE: the commented-out line glued the previous sizes onto the size menu
    updateSelectizeInput(session, "class", choices = "")
    updateSelectizeInput(session, "modelChoices_SE1", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outSEp", choices = "")
    updateSelectizeInput(session, "outSEk", choices = "")
    updateSelectizeInput(session, "outSEclass", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
    updateSelectizeInput(session, "predsCP", choices = "")
    updateSelectizeInput(session, "ltp", choices = "")
    updateSelectizeInput(session, "fta", choices = "")
    updateSelectizeInput(session, "modelChoices_CP1", choices = "")
    updateSelectizeInput(session, "outCPl", choices = "")
    updateSelectizeInput(session, "outCPs", choices = "")
    updateSelectizeInput(session, "outCPdist", choices = "")
    updateSelectizeInput(session, "outCPclass", choices = "")
    updateSelectizeInput(session, "xID", choices = "")
    updateSelectizeInput(session, "COdate", choices = "")
    updateNumericInput(session, "gSearchInterval",  value = NULL)
    updateNumericInput(session, "gSearchMax", value = NULL)
}

  if (eventName == "file_SE"){
    updateSelectizeInput(session, "predsSE", choices = rv$colNames_SE_preds)
    updateSelectizeInput(session, "obsSE", choices = rv$colNames_SE_obs)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol)
    updateTabsetPanel(session, "LoadedDataViz", "Searcher Efficiency")

    if (rv$nsizeclasses > 1){
      updateSelectizeInput(session, "DWPCol", selected = " ")
    }
  }

  if (eventName == "file_SE_clear" ||
      (eventName == "file_SE" & is.null(rv$filename_SE))){
    toReset <- c("file_SE", "predsSE", "obsSE", "outSEp", "outSEk",
                 "outSEclass", "DWPCol", "split_SS", "split_CO",
                 "modelChoices_SE1", "outgclass")
    lapply(toReset, reset)

    scc <- rv$colNames_size
    scs <- rv$sizeCol
    if (is.null(scc)){
      scs <- ""
      scc <- ""
    }
    updateSelectizeInput(session, "predsSE", choices = "")
    updateSelectizeInput(session, "obsSE", choices = "")
    updateSelectizeInput(session, "class", choices = scc, selected = scs)
    updateSelectizeInput(session, "modelChoices_SE1", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outSEp", choices = "")
    updateSelectizeInput(session, "outSEk", choices = "")
    updateSelectizeInput(session, "outSEclass", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
    updateSelectizeInput(session, "DWPCol", choices = "")
    updateTabsetPanel(session, "LoadedDataViz", "Searcher Efficiency")
  }

  if (eventName == "file_CP"){
    updateSelectizeInput(session, "predsCP", choices = rv$colNames_CP_preds)
    updateSelectizeInput(session, "ltp", choices = rv$colNames_ltp)
    updateSelectizeInput(session, "fta", choices = rv$colNames_fta)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )
    updateTabsetPanel(session, "LoadedDataViz", "Carcass Persistence")
  }


  if (eventName == "file_CP_clear"){

    toReset <- c("file_CP", "predsCP", "ltp", "fta", "outCPl", "outCPs",
                 "outCPdist", "outCPclass", "modelChoices_CP1",
                 "split_SS", "split_CO", "outgclass")
    lapply(toReset, reset)

    scc <- rv$colNames_size
    scs <- rv$sizeCol
    if (is.null(scc)){
      scs <- ""
      scc <- ""
    }
    updateSelectizeInput(session, "predsCP", choices = "")
    updateSelectizeInput(session, "ltp", choices = "")
    updateSelectizeInput(session, "fta", choices = "")
    updateSelectizeInput(session, "class", choices = scc, selected = scs)
    updateSelectizeInput(session, "modelChoices_CP1", choices = "")
    updateSelectizeInput(session, "outCPl", choices = "")
    updateSelectizeInput(session, "outCPs", choices = "")
    updateSelectizeInput(session, "outCPdist", choices = "")
    updateSelectizeInput(session, "outCPclass", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
    updateTabsetPanel(session, "LoadedDataViz", "Carcass Persistence")
  }

  if (eventName == "file_SS"){
    updateNumericInput(session, "gSearchInterval", value = rv$SS[["I"]])
    updateNumericInput(session, "gSearchMax", value = rv$SS[["span"]])
    updateSelectizeInput(session, "split_SS", choices = rv$splittable_SS)
    updateTabsetPanel(session, "LoadedDataViz", "Search Schedule")
  }

  if (eventName == "file_SS_clear"){

    toReset <- c("file_SS", "gSearchInterval", "gSearchMax",
      "split_SS", "split_CO", "outgclass")
    lapply(toReset, reset)

    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
    updateTabsetPanel(session, "LoadedDataViz", "Search Schedule")
  }

  if (eventName == "file_DWP"){
    updateSelectizeInput(session, "DWPCol", choices = rv$colNames_DWP)
    if (length(rv$colNames_DWP) == 1){
      updateSelectizeInput(session, "DWPCol", selected = rv$colNames_DWP)
    }
    if (rv$nsizeclasses > 1){
      updateSelectizeInput(session, "DWPCol", selected = rv$colNames_DWP[1])
    }
    updateTabsetPanel(session, "LoadedDataViz", "Density Weighted Proportion")
  }

  if (eventName == "file_DWP_clear"){

    toReset <- c("file_DWP", "DWPCol", "split_SS", "split_CO")
    lapply(toReset, reset)

    updateSelectizeInput(session, "DWPCol", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateTabsetPanel(session, "LoadedDataViz", "Density Weighted Proportion")
  }

  if (eventName == "file_CO"){
    updateSelectizeInput(session, "xID", choices = rv$colNames_xID,
      selected = rv$colNames_xID[1])
    updateSelectizeInput(session, "COdate", choices = rv$colNames_COdates)
    if (length(rv$colNames_COdates) == 1){
      updateSelectizeInput(session, "COdate", choices = rv$colNames_COdates,
        selected = rv$colNames_COdates)
    }
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol)
    updateTabsetPanel(session, "LoadedDataViz", "Carcass Observations")
  }

  if (eventName == "file_CO_clear"){

    toReset <- c("file_CO", "COdate", "split_SS", "split_CO")
    lapply(toReset, reset)

    scc <- rv$colNames_size
    scs <- rv$sizeCol
    if (is.null(scc)){
      scs <- ""
      scc <- ""
    }
    updateSelectizeInput(session, "COdate", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateTabsetPanel(session, "LoadedDataViz", "Carcass Observations")
  }

  if (grepl("load_", eventName)){
    updateSelectizeInput(session, "predsSE", choices = rv$colNames_SE_preds)
    updateSelectizeInput(session, "obsSE", choices = rv$colNames_SE_obs)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol)
    updateNumericInput(session, "gSearchInterval", value = rv$SS[["I"]])
    updateNumericInput(session, "gSearchMax", value = rv$SS[["span"]])
    updateSelectizeInput(session, "predsCP", choices = rv$colNames_CP_preds)
    updateSelectizeInput(session, "ltp", choices = rv$colNames_ltp)
    updateSelectizeInput(session, "fta", choices = rv$colNames_fta)
    updateSelectizeInput(session, "class", choices = rv$colNames_size)
    updateSelectizeInput(session, "xID", choices = rv$colNames_xID, selected = rv$xIDcol)
    updateSelectizeInput(session, "DWPCol", choices = rv$colNames_DWP)
    if (length(rv$colNames_DWP) == 1){
      updateSelectizeInput(session, "DWPCol", selected = rv$colNames_DWP)
    }
    if (rv$nsizeclasses > 1){
      updateSelectizeInput(session, "DWPCol", selected = rv$colNames_DWP[1])
    }
    updateSelectizeInput(session, "COdate", choices = rv$colNames_COdates)
    if (length(rv$colNames_COdates) == 1){
      updateSelectizeInput(session, "COdate",
        choices = rv$colNames_COdates, selected = rv$colNames_COdates
      )
    }
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )

    if (rv$nsizeclasses > 1){
      updateSelectizeInput(session, "DWPCol", selected = " ")
    }
    updateNavbarPage(session, "GenEstApp", selected = "Data Input")
    updateTabsetPanel(session, "LoadedDataViz", "Carcass Persistence")
    updateTabsetPanel(session, "LoadedDataViz", "Search Schedule")
    updateTabsetPanel(session, "LoadedDataViz", "Density Weighted Proportion")
    updateTabsetPanel(session, "LoadedDataViz", "Carcass Observations")
    updateTabsetPanel(session, "LoadedDataViz", "Searcher Efficiency")
  }

  if (eventName == "class"){
    toReset <- c(
      "outCPl", "outCPs", "outCPdist", "outsizeclassCP", "modelChoices_CP1",
      "outSEp", "outSEk", "outsizeclassSE", "modelChoices_SE1",
      "DWPCol", "split_SS", "split_CO", "outgclass")

    lapply(toReset, reset)

    updateSelectizeInput(session, "predsSE", choices = rv$colNames_SE_preds,
      selected = input$predsSE)
    updateSelectizeInput(session, "obsSE", choices = rv$colNames_SE_obs,
      selected = input$obsSE)
    updateSelectizeInput(session, "ltp", choices = rv$colNames_CP_nosel,
      selected = input$ltp)
    updateSelectizeInput(session, "fta", choices = rv$colNames_CP_nosel,
      selected = input$fta)
    updateSelectizeInput(session, "predsCP", choices = rv$colNames_CP_preds,
      selected = input$predsCP)
    updateSelectizeInput(session, "DWPCol", choices = rv$colNames_DWP,
      selected = rv$DWPCol)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol)
  }

  if (eventName == "obsSE"){
    updateSelectizeInput(session, "predsSE", choices = rv$colNames_SE_preds,
        selected = rv$preds_SE)
    updateSelectizeInput(session, "obsSE", choices = rv$colNames_SE_obs,
        selected = rv$obsCols_SE)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )
  }

  if (eventName == "predsSE"){
    updateSelectizeInput(session, "obsSE", choices = rv$colNames_SE_obs,
        selected = rv$obsCols_SE)
    updateSelectizeInput(session, "predsSE", choices = rv$colNames_SE_preds,
        selected = rv$preds_SE)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )
  }

  if (eventName == "ltp"){
    updateSelectizeInput(session, "fta", choices = rv$colNames_fta,
        selected = rv$fta)
    updateSelectizeInput(session, "ltp", choices = rv$colNames_ltp,
        selected = rv$ltp)
    updateSelectizeInput(session, "predsCP", choices = rv$colNames_CP_preds,
        selected = rv$preds_CP)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )
  }

  if (eventName == "fta"){
    updateSelectizeInput(session, "fta", choices = rv$colNames_fta,
        selected = rv$fta)
    updateSelectizeInput(session, "ltp", choices = rv$colNames_ltp,
        selected = rv$ltp)
    updateSelectizeInput(session, "predsCP", choices = rv$colNames_CP_preds,
        selected = rv$preds_CP)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )
  }

  if (eventName == "predsCP"){
    updateSelectizeInput(session, "fta", choices = rv$colNames_fta,
        selected = rv$fta)
    updateSelectizeInput(session, "ltp", choices = rv$colNames_ltp,
        selected = rv$ltp)
    updateSelectizeInput(session, "predsCP", choices = rv$colNames_CP_preds,
        selected = rv$preds_CP)
    updateSelectizeInput(session, "class", choices = rv$colNames_size,
      selected = rv$sizeCol
    )
  }

  if (eventName == "run_SE"){
    updateTabsetPanel(session, "analyses_SE", "Model Comparison")
    updateSelectizeInput(session, "outSEp", choices = rv$modNames_SEp)
    updateSelectizeInput(session, "outSEk", choices = rv$modNames_SEk)
    updateSelectizeInput(session, "outSEclass", choices = rv$sizeclasses)
    updateSelectizeInput(session, "DWPCol", choices = rv$colNames_DWP,
      selected = rv$DWPCol)
    if (length(rv$colNames_DWP) == 1){
      updateSelectizeInput(session, "DWPCol", selected = rv$colNames_DWP)
    }
    reset("outgclass")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
  }

  if (eventName == "run_SE_clear"){
    toReset <- c("outSEp", "outSEk", "outsizeclassSE", #DWPCol,       # reset DWPCol after run_SE_clear?
                 "split_SS", "split_CO", "modelChoices_SE1", "outgclass")
    lapply(toReset, reset)
    updateSelectizeInput(session, "modelChoices_SE1", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outSEp", choices = "")
    updateSelectizeInput(session, "outSEk", choices = "")
    updateSelectizeInput(session, "outSEclass", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
  }

  if (eventName == "outSEclass"){
    updateSelectizeInput(session, "outSEp", choices = rv$modNames_SEp)
    updateSelectizeInput(session, "outSEk", choices = rv$modNames_SEk)
  }

  if (eventName == "run_CP"){
    toReset <- c("outgclass")#, "gSearchInterval", "gSearchMax")
    lapply(toReset, reset)
    updateTabsetPanel(session, "analyses_CP", "Model Comparison")
    updateSelectizeInput(session, "outCPl", choices = rv$modNames_CPl)
    updateSelectizeInput(session, "outCPs", choices = rv$modNames_CPs)
    updateSelectizeInput(session, "outCPdist", choices = rv$modNames_CPdist)
    updateSelectizeInput(session, "outCPclass", choices = rv$sizeclasses)
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
  }

  if (eventName == "run_CP_clear"){
    toReset <- c("outCPl", "outCPs", "outCPdist", "outsizeclassCP",
                 "split_SS", "split_CO", "modelChoices_CP1", "outgclass")
    lapply(toReset, reset)
    updateSelectizeInput(session, "modelChoices_CP1", choices = "")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
    updateSelectizeInput(session, "outCPl", choices = "")
    updateSelectizeInput(session, "outCPs", choices = "")
    updateSelectizeInput(session, "outCPdist", choices = "")
    updateSelectizeInput(session, "outCPclass", choices = "")
    updateSelectizeInput(session, "outgclass", choices = "")
  }

  if (eventName == "outCPclass"){
    updateSelectizeInput(session, "outCPl", choices = rv$modNames_CPl)
    updateSelectizeInput(session, "outCPs", choices = rv$modNames_CPs)
    updateSelectizeInput(session, "outCPdist", choices = rv$modNames_CPdist)
  }

  if (eventName == "run_g"){
    updateSelectizeInput(session, "outgclass", choices = rv$sizeclasses_g)
    updateTabsetPanel(session, "analyses_g", "Summary")
  }

  if (eventName == "run_g_clear"){
    reset("outgclass")
    updateSelectizeInput(session, "outgclass", choices = "")
  }

  if (eventName == "run_M"){
    updateSelectizeInput(session, "xID", choices = rv$colNames_xID, selected = rv$xIDcol)
    updateNumericInput(session, "frac", value = rv$frac)
    updateSelectizeInput(session, "split_SS", choices = rv$splittable_SS)
    updateSelectizeInput(session, "split_CO", choices = rv$colNames_CO)
  }

  if (eventName == "run_M_clear"){
    reset("split_SS")
    reset("split_CO")
    updateSelectizeInput(session, "split_SS", choices = "")
    updateSelectizeInput(session, "split_CO", choices = "")
  }


  if (eventName == "split_M_clear"){
    reset("split_SS")
    reset("split_CO")
    updateSelectizeInput(session, "split_SS", choices = rv$splittable_SS)
    updateSelectizeInput(session, "split_CO", choices = rv$colNames_CO)
  }
}
