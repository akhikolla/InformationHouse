#' @title app utilities
#'
#' @description utility functions for simple app rv management
#'
#' @param rv Reactive values list for the GenEst GUI, created by
#'   \code{\link{initialReactiveValues}}, which calls
#'   \code{\link[shiny]{reactiveValues}}
#'
#' @param toReVal Names of elements in \code{rv} to reset to their factory
#'   setting (as defined by \code{\link{initialReactiveValues}}).
#'
#' @param cols column names to select
#'
#' @param data data table
#'
#' @param modNames names of the model to be split off
#'
#' @param pos position in the name to split off
#'
#' @param sizeclasses names of the carcass classes
#'
#' @param parts the component parts of the model's name
#'
#' @param type "SE" or "CP"
#'
#' @param tab logical for if it's the table output for CP
#'
#' @param sizeCol carcass class column name
#'
#' @param choice carcass class chosen
#'
#' @param colNames_size updated vector of size column names in all needed
#'  tables
#'
#' @name app_utilities
#'
NULL
#' @rdname app_utilities
initialReactiveValues <- function(){
  reactiveValues(
    data_SE = NULL, data_CP = NULL, data_SS = NULL, data_DWP = NULL,
    data_CO = NULL,
    filename_SE = NULL, filename_CP = NULL, filename_SS = NULL,
    filename_DWP = NULL, filename_CO = NULL,
    colNames_SE = NULL, colNames_SE_preds = NULL, colNames_SE_preds0 = NULL,
    toRemove_SE_obs = NULL, toRemove_SE_preds = NULL,
    colNames_SE_obs = NULL, colNames_SE_obs0 = NULL,
    colNames_CP = NULL, colNames_CP_preds = NULL, colNames_CP_preds0 = NULL,
    toRemove_CP_preds = NULL,
    colNames_ltp = NULL, colNames_ltp0 = NULL, toRemove_ltp = NULL,
    colNames_fta = NULL, colNames_fta0 = NULL, toRemove_fta = NULL,

    colNames_SS = NULL, splittable_SS = NULL,
    colNames_DWP = NULL,
    colNames_CO = NULL, colNames_COdates = NULL,
    colNames_size = NULL, colNames_size0 = NULL,

    nsim = 1000, CL = 0.90,

    sizeCol = NULL, toRemove_sizeCol = NULL,
    sizeclasses = NULL, sizeclass = NULL, sizeclass_SE = NULL,
    sizeclass_CP = NULL, sizeclass_g = NULL, sizeclass_M = NULL,
    nsizeclasses = 0,

    obsCols_SE = NULL, preds_SE = NULL, predictors_SE = NULL,
    formula_p = NULL, formula_k = NULL, kFixed = NULL,
    mods_SE = NULL, mods_SE_og = NULL, sizeclasses_SE = NULL,
    outSEpk = NULL, AICcTab_SE = NULL, modOrder_SE = NULL, modNames_SE = NULL,
    modNames_SEp = NULL, modNames_SEk = NULL, modSet_SE = NULL,
    best_SE = NULL, modTab_SE = NULL, modTabPretty_SE = NULL,
    modTabDL_SE = NULL, figH_SE = 800, figW_SE = 800,
    sizeclasses_k = NULL, nsizeclasses_k = NULL, cols_SE = .cols_SE,

    ltp = NULL, fta = NULL, preds_CP = NULL, dist = NULL,
    predictors_CP = NULL, formula_l = NULL, formula_s = NULL,
    mods_CP = NULL, mods_CP_og = NULL, CPdls = NULL, outCPdlsfig = NULL,
    outCPdlstab = NULL, sizeclasses_CP = NULL, AICcTab_CP = NULL,
    modOrder_CP = NULL, modNames_CP = NULL, modNames_CPdist = NULL,
    modNames_CP = NULL, modNames_CPs = NULL, modSet_CP = NULL,
    best_CP = NULL, modTab_CP = NULL, figH_CP = 700, figW_CP = 800,

    M = NULL, Msplit = NULL, unitCol = NULL, colNames_xID = NULL, xID = NULL,
    frac = 1, sizeCol_M = NULL, DWPCol = NULL, COdate = NULL,
    SEmodToUse = NULL, CPmodToUse = NULL,
    split_CO = NULL, split_SS = NULL, nsplit_CO = 0, nsplit_SS = 0,
    figH_M = 600, figW_M = 800,

    SS = NULL, avgSI = NULL, SStemp = NULL, gSearchInterval = NULL,
    gSearchMax = NULL, sizeclasses_g = NULL, nsizeclasses_g = NULL,
    gGeneric = NULL, SEmodToUse_g = NULL, CPmodToUse_g = NULL,
    figH_g = 700, figW_g = 800,

    kCheck = NA, kCheck_g = NA, csvformat = ""
  )
}

#' @rdname app_utilities
#'
reVal <- function(rv, toReVal){
  if("xID" %in% toReVal){
    rv$colnames_xID  <- NULL
  }
  if("nsplit_CO" %in% toReVal){
    rv$nsplit_CO <- 0
  }
  if("nsplit_SS" %in% toReVal){
    rv$nsplit_SS <- 0
  }

  if("SS" %in% toReVal){
    rv$SS <- NULL #seq(0, 364, 7)
  }
  if("gSearchInterval " %in% toReVal){
    rv$gSearchInterval  <- NULL#7
  }
  if("gSearchMax" %in% toReVal){
    rv$gSearchMax <- NULL#364
  }
  if("figH_SE" %in% toReVal){
    rv$figH_SE <- 800
  }
  if("figW_SE" %in% toReVal){
    rv$figW_SE <- 800
  }
  if("figH_CP" %in% toReVal){
    rv$figH_CP <- 800
  }
  if("figW_CP" %in% toReVal){
    rv$figW_CP <- 800
  }
  if("figH_M" %in% toReVal){
    rv$figH_M <- 600
  }
  if("figW_M" %in% toReVal){
    rv$figW_M <- 800
  }
  if("figH_M" %in% toReVal){
    rv$figH_g <- 700
  }
  if("figW_M" %in% toReVal){
    rv$figW_g <- 800
  }
  rv
}

#' @title Read in csv files in either format
#'
#' @description Handle reading in of a csv that is either comma-decimal or
#'   semicolon-comma separation style
#'
#' @param path file path
#'
#' @return read in data table
#'
#' @export
#'
readCSV <- function(path){
  ef <- function(x){"_BAD_READ_"}
  out1 <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
            error = ef, warning = ef)

  out2 <- tryCatch(read.csv2(path, stringsAsFactors = FALSE),
            error = ef, warning = ef)
  if (is.null(attr(out1, "class")) & is.null(attr(out2, "class"))){
   stop("File not found or not formatted as a .csv")
  }
  if ("data.frame" %in% attr(out1, "class")){
    if (is.null(attr(out2, "class"))){
      return(out1)
    }
    if ("data.frame" %in% attr(out2, "class")){
      if (ncol(out2) == 1){
        return(out1)
      }
    }
  }
  if ("data.frame" %in% attr(out2, "class")){
    if (is.null(attr(out1, "class"))){
      return(out2)
    }
    if ("data.frame" %in% attr(out1, "class")){
      if (ncol(out1) == 1){
        return(out2)
      }
    }
  }
  return(out1)
}

#' @title Prepare predictors based on inputs
#'
#' @description Prepare predictor inputs from the app for use in the model
#'   function
#'
#' @param preds predictors, as input to the app
#'
#' @return prepared predictors (or 1 if no predictors)
#'
#' @export
#'
prepPredictors <- function(preds = NULL){

 out <- paste(preds, collapse = "*")
 if (is.null(preds)){
   out <- 1
 }
 return(out)
}

#' @rdname app_utilities
#'
setkNeed <- function(rv){
  textout <- "no"
  if(length(rv$obsCols_SE) == 1 & any(is.na(rv$kFixed))){
    textout <- "yes"
  }
  return(renderText(textout))
}

#' @title Select the date columns from a data table
#'
#' @description Simple function to facilitate selection of date columns from
#'   a data table
#'
#' @param data data table potentially containing columns that could be
#'   coerced (via \code{checkDate()}) into a properly formatted date
#'
#' @return column names of columns that can be coerced to dates
#'
#' @export
#'
dateCols <- function(data){

  ncols <- ncol(data)
  dateTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    temp <- tryCatch(
              checkDate(data[ , coli]),
              error = function(x){FALSE}
            )
    dateTF[coli] <- is.Date(temp)
  }
  out <- colnames(data)[dateTF]
  return(out)
}

#' @title Select the potential carcass class columns from a data table
#'
#' @description Simple function to facilitate selection of columns that could
#'   be carcass class values from a data table.
#'
#' @param data data table
#'
#' @return column names of columns that can be carcass class values
#'
#' @export
#'
sizeCols <- function(data){

  if (is.null(data)){
    return(NULL)
  }
  ncols <- ncol(data)
  scTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    tmp <- data[ , coli]
    if (length(unique(tmp)) < nrow(data) & length(unique(tmp)) > 1){
      scTF[coli] <- TRUE
    } else{
      scTF[coli] <- FALSE
    }
  }
  out <- colnames(data)[scTF]
  return(out)
}


#' @title Select the DWP-ok columns from a data table
#'
#' @description Simple function to facilitate selection of columns that could
#'   be DWP values from a data table
#'
#' @param data data table
#'
#' @return column names of columns that can be DWP values
#'
#' @export
#'
DWPCols <- function(data){
  ncols <- ncol(data)
  dwpTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    tmp <- data[ , coli]
    if (!is.factor(tmp) && is.numeric(tmp) &&( all(tmp > 0) & all(tmp <= 1))){
      dwpTF[coli] <- TRUE
    } else{
      dwpTF[coli] <- FALSE
    }
  }
  out <- colnames(data)[dwpTF]
  return(out)
}

#' @title Select the predictor-ok columns from a data table
#'
#' @description Simple function to facilitate selection of columns that could
#'   be predictors for SE or CP models from a data table
#'
#' @param data data table
#'
#' @return column names of columns that can be predictors
#'
#' @export
#'
predsCols <- function(data){
  ncols <- ncol(data)
  predTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    tmp <- data[ , coli]
    cont <- FALSE
    if (is.numeric(tmp) && any(na.omit(tmp %% 1 != 0))){
      cont <- TRUE
    }
    if (length(unique(tmp)) == nrow(data)){
      reps <- FALSE
    } else{
      reps <- TRUE
    }
    if (grepl("[-.]", colnames(data)[coli])){
      okName <- FALSE
    } else{
      okName <- TRUE
    }
    if (!any(is.na(tmp)) & !cont & reps & okName){
      predTF[coli] <- TRUE
    } else{
      predTF[coli] <- FALSE
    }
  }
  out <- colnames(data)[predTF]
  return(out)
}

#' @title Select the columns from a data table that could be SE observations
#'
#' @description Simple function to facilitate selection of columns that could
#'   be observations for an SE model
#'
#' @param data data table
#'
#' @return column names of columns that can be observations
#'
#' @export
#'
obsCols_SE <- function(data){
  ncols <- ncol(data)
  obsTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    tmp <- na.omit(data[ , coli])
    if (any(tmp == 0 | tmp == 1) & all(tmp == 0 | tmp == 1)){
      obsTF[coli] <- TRUE
    } else{
      obsTF[coli] <- FALSE
    }
  }
  out <- colnames(data)[obsTF]
  return(out)
}

#' @title Select the columns from a data table that could be CP Last Time
#'   Present observations
#'
#' @description Simple function to facilitate selection of columns that could
#'   be Last Time Present observations for a CP model
#'
#' @param data data table
#'
#' @return column names of columns that can be observations
#'
#' @export
#'
obsCols_ltp <- function(data){
  ncols <- ncol(data)
  obsTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    tmp <- data[ , coli]
    if (is.numeric(tmp) && is.finite(tmp) && all(na.omit(tmp) >= 0)){
      obsTF[coli] <- TRUE
    } else{
      obsTF[coli] <- FALSE
    }
  }
  out <- colnames(data)[obsTF]
  return(out)
}


#' @title Select the columns from a data table that could be CP First Time
#'   Absent observations
#'
#' @description Simple function to facilitate selection of columns that could
#'   be First Time Absent observations for a CP model
#'
#' @param data data table
#'
#' @return column names of columns that can be observations
#'
#' @export
#'
obsCols_fta <- function(data){
  ncols <- ncol(data)
  obsTF <- rep(NA, ncols)
  for (coli in 1:ncols){
    tmp <- data[ , coli]
    if (is.numeric(tmp) && all(na.omit(tmp) > 0)){
      obsTF[coli] <- TRUE
    } else{
      obsTF[coli] <- FALSE
    }
  }
  out <- colnames(data)[obsTF]
  return(out)
}


#' @title Remove selected columns from column names
#'
#' @description Simple function to facilitate removal of columns selected
#'
#' @param colNames column names from which some could be removed
#'
#' @param selCols selected columns to be removed
#'
#' @return column names without selected columns
#'
#' @export
#'
removeCols <- function(colNames, selCols){
  which_sel <- which(colNames %in% selCols)
  if (length(which_sel) > 0){
    out <- colNames[-which_sel]
  } else{
    out <- colNames
  }
  return(out)
}

#' @rdname app_utilities
#'
updateColNames_size <- function(rv){

  SECPCO <- NULL
  SE <- sizeCols(rv$data_SE)
  CP <- sizeCols(rv$data_CP)
  CO <- names(rv$data_CO)

  SECP <- which(SE %in% CP)
  SECO <- which(SE %in% CO)
  CPSE <- which(CP %in% SE)
  CPCO <- which(CP %in% CO)
  COSE <- which(CO %in% SE)
  COCP <- which(CO %in% CP)
  alltogether <- c(SECP, SECO, CPSE, CPCO, COSE, COCP)

  if (length(alltogether) == 0){
    if (is.null(SE) + is.null(CP) + is.null(CO) == 2){
      SECPCO <- unique(c(SE, CP, CO))
    }
  } else{
    if (is.null(SE) + is.null(CP) + is.null(CO) == 1){
      SECPCOa <- c(SE[SECP], SE[SECO], CP[CPSE], CP[CPCO], CO[COSE], CO[COCP])
      SECPCO <- unique(SECPCOa)
    } else{
      SECP <- SE[SE %in% CP]
      SECPCO <- CO[CO %in% SECP]
    }
  }

  return(SECPCO)
}

#' @rdname app_utilities
#'
selectData <- function(data, cols){
  if (is.null(data)){
    return(NULL)
  }

  colNames <- colnames(data)
  selectedTab <- data[ , which(colNames %in% cols)]
  selectedDF <- data.frame(selectedTab)
  if (length(cols) == 1){
    colnames(selectedDF) <- cols
  }
  return(selectedDF)
}

#' @rdname app_utilities
#'
modNameSplit <- function(modNames, pos){
  modNames_split <- modNames
  nmod <- length(modNames)
  if (nmod > 0){
    for (modi in 1:nmod){
      modNames_split[modi] <- strsplit(modNames[modi], "; ")[[1]][pos]
    }
  }
  modNames_split <- gsub("NULL", "s ~ 1", modNames_split)
  modNames_split <- gsub("~ 1", "~ constant", modNames_split)
  modNames_split <- gsub("dist:", "", modNames_split)
  return(modNames_split)
}

#' @title Count the minimum number of carcasses in the cells
#'
#' @description Count the minimum number of carcasses in all of the cells
#'   within a \code{_SetSize} model complex
#'
#' @param mods model output from the \code{_SetSize} version of a function
#'
#' @return the minimum number of carcasses in the cells
#'
#' @export
#'
countCarcs <- function(mods){
  nsizeclasses <- length(mods)
  nmods <- sum(unlist(lapply(mods, length)))
  if (nsizeclasses > 0 & nmods > 0){
    ncarc <- rep(NA, nmods)
    counter <- 0
    for (sci in names(mods)){
      for (modi in 1:length(mods[[sci]])){
        counter <- counter + 1
        if (!grepl("Failed model fit", mods[[sci]][[modi]][1])){
          ncarc[counter] <- min(table(mods[[sci]][[modi]]$carcCell))
        }
      }
    }
    ncarc <- min(na.omit(ncarc))
  } else {
    ncarc <- Inf
  }
  return(ncarc)
}

#' @rdname app_utilities
#'
prepSizeclassText <- function(sizeclasses){
  return(renderText(paste(sizeclasses, collapse = " ")))
}

#' @rdname app_utilities
#'
modNamePaste <- function(parts, type = "SE", tab = FALSE){
  if (tab & parts[1] == " exponential"){
    out <- paste(c(parts[1:2], "NULL"), collapse = "; ")
  } else{
    out <- paste(parts, collapse = "; ")
  }
  if (type == "CP"){
    out <- paste("dist:", out, sep = "")
  }
  out <- gsub("~ constant", "~ 1", out)
  return(out)
}

#' @title Produce the options for the distributions in the CP model
#'
#' @description Simply make the named list for the disributions in the CP
#'   model
#'
#' @return list with named elements of the distributions
#'
#' @export
#'
CPdistOptions <- function(){
  list("exponential" = "exponential", "weibull" = "weibull",
    "lognormal" = "lognormal", "loglogistic" = "loglogistic"
  )
}

#' @rdname app_utilities
#'
plotNA <- function(type = "model"){
  if (type == "model"){
    badText <- "Selected model was not fit successfully."
  }
  if (type == "split"){
    badText <- "Second split too fine for plotting. Consider transposing."
  }
  plot(1, 1, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "",
    ylab = "", ylim = c(0, 1), xlim = c(0,1))
  text(0.01, 0.9, badText, adj = 0)
}

#' @rdname app_utilities
#'
updateSizeclasses <- function(data, sizeCol){
  if (is.null(sizeCol)){
    return("all")
  }
  return(as.character(sort(unique(data[ , sizeCol]))))
}

#' @rdname app_utilities
#'
pickSizeclass <- function(sizeclasses, choice){

  sizeclass <- NULL
  if (!(choice %in% sizeclasses)){
    choice <- sizeclasses[1]
  }
  sizeclass <- sizeclasses[which(sizeclasses == choice)]
  return(sizeclass)
}

#' @rdname app_utilities
#'
updatesizeCol <- function(sizeCol, colNames_size){
  if (!is.null(sizeCol)){
    if (!(sizeCol %in% colNames_size)){
      NULL
    } else{
      sizeCol
    }
  } else {
    NULL
  }
}
