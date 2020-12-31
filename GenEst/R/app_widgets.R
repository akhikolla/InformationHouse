#' @title Create and manage widgets for data input, function execution, data output
#'
#' @description This is a generalized function for creating a data input 
#'   widget used in the GenEst GUI, based on the data type (\code{dataType}).
#'   Included within the widget is a conditional panel that allows removal of
#'   the specific data file (and clearing of all downstream models) once
#'   it has been loaded.
#'
#' @param dataType Toggle control for the model type of the widget. One of 
#'   "SE", "CP", "SS", "DWP", or "CO".  
#'
#' @param set Name of data set. One of "RP", "RPbat", "cleared", "powerTower",
#'   "PV", "trough", "mock"
#'
#' @param inType Toggle control for the input type of the widget. One of
#'   "nsim", "CL", "class", "obsSE", "predsSE", "kFixed", "ltp", "fta",
#'   "predsCP", "dist", "xID", "frac", "DWPCol", "COdate", "gSearchInterval", or
#'   "gSearchMax".
#'
#' @param Name Name (id) of the widget created.
#'
#' @param Fun Function name (as character) used to create the widget.
#'
#' @param Label Label presented to the user in the UI for the widget.
#'
#' @param Args List of any additional arguments to be passed to the widget
#'   creating function.
#'
#' @param modType Toggle control for the model type of the widget. One of
#'   "SE", "CP", "M", or "g".
#'
#' @param Condition Condition under which the widget is present to the user.
#'
#' @param mods Model Set Size object (from the reactive values list).
#'
#' @param sci Name of carcass class element
#'
#' @param sizeclasses Vector of carcass class names (from the reactive values
#'   list).
#'
#' @return HTML for the data input widget.
#' @name app_widgets
dataInputWidget <- function(dataType){

  if (!dataType %in% c("SE", "CP", "SS", "DWP", "CO")){
    stop(paste0("input dataType (", dataType, ") not supported"))
  }

  okft <- c("text/csv", "text/comma-separated-values", ".csv")
  llab <- "Load file"
  clab <- "Clear file"
  cstyle <- cButtonStyle()

  Label <- switch(dataType, "SE" = "Searcher Efficiency (SE)",
                            "CP" = "Carcass Persistence (CP)",
                            "SS" = "Search Schedule (SS)",
                            "DWP" = "Density Weighted Proportion (DWP)",
                            "CO" = "Carcass Observation (CO)")
  ButtonName <- switch(dataType, "SE" = "file_SE",
                                 "CP" = "file_CP",
                                 "SS" = "file_SS",
                                 "DWP" = "file_DWP",
                                 "CO" = "file_CO")
  ClearButtonName <- switch(dataType, "SE" = "file_SE_clear",
                                      "CP" = "file_CP_clear",
                                      "SS" = "file_SS_clear",
                                      "DWP" = "file_DWP_clear",
                                      "CO" = "file_CO_clear")
  ConditionPrefix <- switch(dataType, "SE" = "output.filename_",
                                      "CP" = "output.filename_",
                                      "SS" = "output.filename_",
                                      "DWP" = "output.filename_",
                                      "CO" = "output.filename_")
  PanelCondition <- paste0(ConditionPrefix, dataType, " != null")

  list(
    h5(b(Label)),
    fileInput(ButtonName, label = NULL, accept = okft, buttonLabel = llab), 
    conditionalPanel(condition = PanelCondition,
      actionButton(ClearButtonName, clab, style = cstyle) 
    )
  )
}

#' @rdname app_widgets
#'
dataDownloadWidget <- function(set){

  if (!set %in% c("cleared", "RP", "RPbat", "powerTower", "PV", "trough",
                  "mock")){
    stop(paste0("input set (", set, ") not supported"))
  }

  setNames <- c("cleared" = "Wind---Cleared plots, bats + birds",
                "RP" = "Wind---Road and pad searches, bats + birds",
                "RPbat" = "Wind---Road and pad searches, bats",
                "powerTower" = "Solar---Power tower",
                "PV" = "Solar---Photovoltaic (PV)",
                "trough" = "Solar---Trough",
                "mock" = "Mock data")

  setName <- setNames[set]
  fluidRow(
    column(6, h4(setName)),
    column(2, actionButton(paste0("load_", set), "Load Data")),
  )
}


#' @rdname app_widgets
#'
modelInputWidget <- function(inType){

  if (!inType %in% c("nsim", "CL", "class", "obsSE", "predsSE",
                     "kFixedInput", "ltp", "fta", "predsCP", "dist",
                     "xID", "frac", "DWPCol", "COdate",
                     "gSearchInterval", "gSearchMax")){
    stop(paste0("input inType (", inType, ") not supported"))
  }

  Name <- inType

  Label <- switch(inType, 
             "nsim" = "Number of Iterations:",
             "CL" = "Confidence Level:", 
             "class" = "Carcass Class Column (optional):",
             "obsSE" = "Observations:", 
             "predsSE" = "Predictor Variables (optional):",
             "kFixedInput" = NULL, 
             "ltp" = "Last Time Present:", 
             "fta" = "First Time Absent:", 
             "predsCP" = "Predictor Variables (optional):",
             "dist" = "Distributions to Include",
             "xID" = "Carcass ID Column (CO)",
             "frac" = "Fraction of Facility Surveyed:",
             "DWPCol" = "Density Weighted Proportion:",
             "COdate" = "Date Found:",
             "gSearchInterval" = "Search Interval (days):",
             "gSearchMax" = "Total Span of Monitoring (days):")

  widgetFun <- switch(inType, 
                 "nsim" = "numericInput",
                 "CL" = "numericInput", 
                 "class" = "selectizeInput",
                 "obsSE" = "selectizeInput", 
                 "predsSE" = "selectizeInput", 
                 "kFixedInput" = "htmlOutput", 
                 "ltp" = "selectizeInput", 
                 "fta" = "selectizeInput", 
                 "predsCP" = "selectizeInput", 
                 "dist" = "checkboxGroupInput",
                 "xID" = "selectizeInput",
                 "frac" = "numericInput",
                 "DWPCol" = "selectizeInput",
                 "COdate" = "selectizeInput",
                 "gSearchInterval" = "numericInput",
                 "gSearchMax" = "numericInput")

  Args <- switch(inType, 
            "nsim" = list(value = 1000, min = 1, max = 10000, step = 1),
            "CL" = list(value = 0.90, min = 0, max = 0.999, step = 0.001),
            "class" = list(c("No data input yet"), multiple = TRUE,
                               options = list(maxItems = 1)), 
            "obsSE" = list(c("No SE data input yet"), multiple = TRUE), 
            "predsSE" = list(c("No SE data input yet"), multiple = TRUE),
            "kFixedInput" = list(NULL),
            "ltp" = list(c("No CP data input yet"), multiple = TRUE,
                      options = list(maxItems = 1)),
            "fta" = list(c("No CP data input yet"), multiple = TRUE,
                      options = list(maxItems = 1)),
            "predsCP" = list(c("No CP data input yet"), multiple = TRUE),
            "dist" = list(choices = CPdistOptions(),
                        selected = unlist(CPdistOptions()), inline = TRUE),
            "xID" =  list(c("No carcass data input yet"), multiple = TRUE,
                            options = list(maxItems = 1)),
            "frac" = list(value = 1.0, min = 0.01, max = 1.0, step = 0.01),
            "DWPCol" = list(c("No DWP data input yet"), multiple = TRUE,
                            options = list(maxItems = 1)),
            "COdate" = list(c("No carcass data input yet"), multiple = TRUE,
                                  options = list(maxItems = 1)),
            "gSearchInterval" = list(value = NULL, min = 1, max = 400, step = 1),
            "gSearchMax" =  list(value = NULL, min = 1, max = 1000, step = 1))

  Condition <- switch(inType, 
                 "nsim" = NULL, 
                 "CL" = NULL, 
                 "class" = NULL,
                 "obsSE" = NULL, 
                 "predsSE" = NULL, 
                 "kFixedInput" = NULL, 
                 "ltp" = NULL, 
                 "fta" = NULL, 
                 "predsCP" = NULL, 
                 "dist" = NULL,
                 "xID" = NULL,
                 "frac" = NULL,
                 "DWPCol" = "output.DWPNeed == 'yes'",
                 "COdate" = NULL,
                 "gSearchInterval" = NULL,
                 "gSearchMax" = NULL)

  widgetMaker(Condition, Name, widgetFun, Label, Args)
}


#' @rdname app_widgets
#'
widgetMaker <- function(Condition, Name, Fun, Label, Args){
  allArgs <- Args
  if (Fun == "htmlOutput"){
    allArgs[["id"]] <- Name
  } else{
    allArgs[["inputId"]] <- Name
  }
  allArgs[["label"]] <- Label
  if(is.null(Condition) || Condition == ""){
    do.call(what = Fun, args = allArgs)
  } else {
    conditionalPanel(condition = Condition, 
      do.call(what = Fun, args = allArgs)
    )
  }
}

#' @rdname app_widgets
#'
modelRunWidget <- function(modType){

  if (!modType %in% c("SE", "CP", "M", "g")){
    stop(paste0("input modType (", modType, ") not supported"))
  }

  rName <- switch(modType,
             "SE" = "run_SE",
             "CP" = "run_CP",
             "M" = "run_M",
             "g" = "run_g")
  rLabel <- switch(modType, 
              "SE" = "Run Model",
              "CP" = "Run Model",
              "M" = "Estimate",
              "g" = "Estimate")
  rCondition <- switch(modType, 
                  "SE" = "input.obsSE != null",
                  "CP" = "input.ltp != null & input.fta != null",
                  "M" = "input.modelChoices_SE1 != null &
                         input.modelChoices_CP1 != null &
                         output.sizeclasses_SE == output.sizeclasses_CP & 
                         output.filename_SS != null &  output.kNeed != 'yes' &
                         input.DWPCol != null &
                         input.COdate != null",
                  "g" = "input.modelChoices_SE1 != null &
                         input.modelChoices_CP1 != null &
                         output.kNeed != 'yes' &
                         output.sizeclasses_SE == output.sizeclasses_CP")


  cName <- switch(modType, 
             "SE" = "run_SE_clear",
             "CP" = "run_CP_clear",
             "M" = "run_M_clear",
             "g" = "run_g_clear")
  cLabel <- switch(modType, 
              "SE" = "Clear Model",
              "CP" = "Clear Model",
              "M" = "Clear Estimate",
              "g" = "Clear Estimate")
  cCondition <- switch(modType, 
                  "SE" = "output.SEModDone == 'OK'",
                  "CP" = "output.CPModDone == 'OK'",
                  "M" = "output.MModDone == 'OK'",
                  "g" = "output.gModDone == 'OK'")

  rButton <- conditionalPanel(condition = rCondition,
               br(), 
               actionButton(rName, rLabel)          
             )
  cButton <- conditionalPanel(condition = cCondition,
               actionButton(cName, cLabel, style = cButtonStyle()),
               br(), br()          
             )

  preText <- preTextMaker(modType)
  lpre <- length(preText)
  preText[[lpre + 1]] <- rButton
  preText[[lpre + 2]] <- cButton
  preText
}


#' @rdname app_widgets
#'
preTextMaker <- function(modType){

  if (!modType %in% c("SE", "CP", "M", "g")){
    stop(paste0("input modType (", modType, ") not supported"))
  }

  Condition <- switch(modType, 
     "SE" = "input.obsSE == null",
     "CP" = "input.ltp == null | input.fta == null",
     "M" = c("input.xID == null",
              "input.modelChoices_SE1 == null |
             input.modelChoices_CP1 == null |
             output.sizeclasses_SE != output.sizeclasses_CP",
             "output.filename_SS == null",
             "input.modelChoices_SE1 != null &
             input.modelChoices_CP1 != null &
             output.sizeclasses_SE == output.sizeclasses_CP & 
             (input.DWPCol == null | input.COdate == null)",
             "output.kNeed == 'yes' &
             input.modelChoices_SE1 != null"),
     "g" = c("input.modelChoices_SE1 == null |
             input.modelChoices_CP1 == null |
             output.sizeclasses_SE != output.sizeclasses_CP",
             "output.kNeed == 'yes' & 
             input.modelChoices_SE1 != null")
  )

  Text <- switch(modType, 
             "SE" = "Select observation columns to run model",
             "CP" = "Select observation columns to run model",
             "M" = c("Select carcass ID column to run model",
                      "Select SE and CP models fit to matching carcass classes to
                      run model",
                     "Input Search Schedule data to run model",
                     "Select input columns to run model",
                     "A value for k is required to estimate mortality.
                      Return to Search Efficiency tab and fix k."),
             "g" = c("Select SE and CP models fit to matching size
                     classes to run model",
                     "A value for k is required to estimate detection
                     probability. Return to Search Efficiency tab and fix k.")
           )

  out <- vector("list", length(Condition))
  for (i in 1:length(Condition)){
    out[[i]] <- conditionalPanel(condition = Condition[i], br(), center(em(Text[i])))
  }
  out
}
#' @rdname app_widgets
#'
modelOutputWidget <- function(modType){

  if (!modType %in% c("SE", "CP", "M", "g")){
    stop(paste0("input modType (", modType, ") not supported"))
  }

  Condition <- switch(modType,
                 "SE" = "output.SEModDone == 'OK'",
                 "CP" = "output.CPModDone == 'OK'",
                 "M" = "output.MModDone == 'OK'",
                 "g" = "output.gModDone == 'OK' & 
                       output.sizeclass_gyn == 'YES'")

  Header <- switch(modType,
                 "SE" = "Table & Figure Selection:",
                 "CP" = "Table & Figure Selection:",
                 "M" = "Splitting Mortality:",
                 "g" = "Table & Figure Selection:")

  sCondition <- switch(modType,
                  "SE" = c("output.sizeclass_SEyn == 'YES'", "", ""),
                  "CP" = c("output.sizeclass_CPyn == 'YES'", rep("", 3)),
                  "M" = c("", ""),
                  "g" = "")

  sName <- switch(modType,
             "SE" = c("outSEclass", "outSEp", "outSEk"),
             "CP" = c("outCPclass", "outCPdist", "outCPl", "outCPs"),
             "M" = c("split_SS", "split_CO"),
             "g" = "outgclass")

  sLabel <- switch(modType,
              "SE" = c("Carcass Class:", "p Model:", "k Model:"),
              "CP" = c("Carcass Class:", "Distribution:", "Location:",
                       "Scale:"),
              "M" = c("Search Schedule (SS) Variable:",
                      "Carcass Observation (CO) Variable:"),
              "g" = "Carcass Class:")

  sArgs <- switch(modType,
              "SE" = list(list(choices = " ", multiple = FALSE),
                          list(choices = " ", multiple = FALSE),
                          list(choices = " ", multiple = FALSE)),
              "CP" = list(list(choices = " ", multiple = FALSE),
                          list(choices = " ", multiple = FALSE),
                          list(choices = " ", multiple = FALSE),
                          list(choices = " ", multiple = FALSE)),
              "M" = list(list(choices = " ", multiple = TRUE, 
                           options = list(maxItems = 1)),
                         list(choices = " ", multiple = TRUE, 
                           options = list(maxItems = 2))),
              "g" = list(list(choices = " ", multiple = FALSE)))

  
  subs <- vector("list", length = length(sCondition))
  for(i in 1:length(sCondition)){
    subs[[i]] <- widgetMaker(sCondition[i], sName[i], "selectizeInput", 
                   sLabel[i], sArgs[[i]]
                 )
  }

  aText <- NULL
  splitButtons <- NULL
  if (modType == "M"){
    aText <- list(em("Max. two total splits, max. one schedule-based split"),
               br(), br()
             )
    splitButtons <- splitButtonWidget()
  }

  conditionalPanel(condition = Condition,
    b(u(big(Header))), br(), br(), aText, subs, splitButtons
  )
}

#' @rdname app_widgets
#'
splitButtonWidget <- function(){
  list(
    fluidRow(
      column(width = 6, actionButton("split_M", "Split Estimate")),
      column(width = 6,
        conditionalPanel(
          condition = "output.MSplitDone == 'OK' & output.nMSplits > 1",
          actionButton("transpose_split", "Transpose")
        )
      )
    ),
    conditionalPanel(condition = "output.MSplitDone == 'OK'", 
        actionButton("split_M_clear", "Clear Split", style = cButtonStyle())
    )
  )
}


#' @rdname app_widgets
#'
modelSelectionWidget <- function(mods, modType){

  if (!any(attr(mods, "class") %in% c("cpmSetSize", "pkmSetSize"))){
    stop("mods must be a cpmSetSize or pkmSetSize object")
  }
  if (!modType %in% c("SE", "CP")){
    stop(paste0("input modType (", modType, ") not supported"))
  }
  nsizeclasses <- length(mods)
  menuHeader <- modelSelectionWidgetHeader(mods)
  modelMenu <- menuHeader
  if (length(mods) > 0){
    for(sci in names(mods)){
      modelMenuRow <- modelSelectionWidgetRow(mods, modType, sci)
      modelMenu <- paste(modelMenu, modelMenuRow)  
    }
  }
  return(renderUI({HTML(modelMenu)}))
}

#' @rdname app_widgets
#'
modelSelectionWidgetHeader <- function(mods){

  if (!any(attr(mods, "class") %in% c("cpmSetSize", "pkmSetSize"))){
    stop("mods must be a cpmSetSize or pkmSetSize object")
  }
  nsizeclasses <- length(mods)
  if (nsizeclasses == 1){
    ntext <- "model"
  } else if (nsizeclasses > 1){
    ntext <- "models"
  } else{
    stop("nsizeclasses input is improper")
  }
  prefix <- "Select "
  suffix <- " for mortality and detection probability estimation"
  menuHeader <- em(paste0(prefix, ntext, suffix))

  menuBreak <- NULL
  if (nsizeclasses > 1){
    menuBreak <- br("")
  } 
  paste(menuHeader, menuBreak)
}

#' @rdname app_widgets
#'
modelSelectionWidgetRow <- function(mods, modType, sci){

  if (!any(attr(mods, "class") %in% c("cpmSetSize", "pkmSetSize"))){
    stop("mods must be a cpmSetSize or pkmSetSize object")
  }
  if (!modType %in% c("SE", "CP")){
    stop(paste0("input modType (", modType, ") not supported"))
  }
  sizeclasses <- names(mods)
  nsizeclasses <- length(mods)
  if (! sci %in% sizeclasses) stop(sci, "not in carcass classes")
  if (modType == "SE"){
    AICcTab <- aicc(mods[[sci]], quiet = TRUE)
  }
  if (modType == "CP"){
    AICcTab <- aicc(mods[[sci]], quiet = TRUE)
  }
  modOrder <- as.numeric(row.names(AICcTab))
  modNames <- names(mods[[sci]])[modOrder]
  modNames <- gsub("; NULL", "", modNames)
  modNames <- gsub("dist: ", "", modNames)
  modNames <- gsub("~ 1", "~ constant", modNames)

  modNames_nchar <- nchar(modNames)
  modNames_maxchar <- max(modNames_nchar)
  modNames_nspaces <- modNames_maxchar - modNames_nchar + 10
  modSpaces <- sapply(modNames_nspaces, 
                     function(x){paste(rep(" ", x), collapse = "")}
                   )
  modDeltaAICcs <- AICcTab[ , "\u0394AICc"]
  modLabels <- paste0(modNames, " (\u0394AICc: ", modDeltaAICcs, ")")
  names(modNames) <- modLabels
  labels_nchar <- nchar(modLabels)
  labels_maxchar <- max(labels_nchar)
  widthval <- max(c(400, labels_maxchar * 7 + 20))
  widthtxt <- paste0(widthval, "px")
  mtuText <- paste0("modelChoices_", modType, which(names(mods) == sci))
  scText <- paste0(sci, ":")
  if (nsizeclasses == 1){
    scText <- ""
  }

  selectizeInput(mtuText, scText, modNames, multiple = TRUE, width = widthtxt, 
    options = list(maxItems = 1)
   )
}

#' @rdname app_widgets
#'
kFixedWidget <- function(sizeclasses){
  widgetHeader <- kFixedWidgetHeader(sizeclasses)

  kFixedMenu <- widgetHeader
  for(sci in sizeclasses){
    kFixedRow <- kFixedWidgetRow(sizeclasses, sci)
    kFixedMenu <- paste(kFixedMenu, kFixedRow)
  }
  renderUI({HTML(kFixedMenu)})
}

#' @rdname app_widgets
#'
kFixedWidgetHeader <- function(sizeclasses){
  nsizeclasses <- length(sizeclasses)
  fluidRow(column(width = 8, align = "center", b("Fixed k (optional)")))
}


#' @rdname app_widgets
#'
kFixedWidgetRow <- function(sizeclasses, sci){
  if (! (sci %in% sizeclasses)) stop(sci, " not in sizeclasses")
  mvText <- paste0("kFixed_val_", sci)
  scText <- paste0(sci, ":")
  rowName <- paste0("string_", sci)

  rowNameStyle <- style(type = "text/css", 
                    paste0("#", rowName, " { margin-top: 10px;}"))
  numStyle <- style(type = "text/css",
                   paste0("#", mvText, " { margin-top: -15px;}"))

  spacerCol <- column(width = 1, div(""))
  rowNameCol <- column(width = 1,  
                  div(id = rowName, b(scText)), align = "right", 
                  rowNameStyle
                )
  if (length(sizeclasses) == 1){
    spacerCol <- NULL
    rowNameCol <- NULL
  }
  numCol <- column(width = 5, align = "center",
              numericInput(mvText, "", value = "", min = 0, max = 1, step = 0.001),
              numStyle
            )

  fluidRow(spacerCol, rowNameCol, numCol)
}


