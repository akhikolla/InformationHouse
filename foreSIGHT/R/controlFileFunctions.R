# parameter variation definitions
# link namelist input with shortened version in modelTag
parVariationDef <- c(ann = "annual",
                     har = "harmonic",
                     har26 = "harmonic",
                     seas = "seasonal")

# There are multiple harmonic models - define the harmonic tags for each variable 
# this is required to create the modelTag from namelist 
# Anjana: this can be eliminated once all the modelTags for "Temp", "PET", and "Radn" are all updated to be "har" instead of "har26"
harDefinition <- c(P = "har",
                   Temp = "har26",
                   PET = "har26",
                   Radn = "har26")

# TO DO: add an argument 'compatibleAtts' (default FALSE)
# if set to 'TRUE' all attributes of the variable will be shown in the printed data.frame 
# - with T/F to indicate if the attribute is compatible to use with the model

#' Prints the available stochastic model options
#'
#' \code{viewModels} prints the stochastic model options available for the different hydroclimatic variables in foreSIGHT.
#'                   These options may be used to create an controlFile for input to function \code{generateScenarios}.
#' @param variable String; the variable name. Type \code{viewModels()} without arguments to view the valid variable names.
#' @param compatibleAtts TRUE/FALSE (default is FALSE). Whether the attributes compatible with each model should be printed.
#' @seealso \code{writeControlFile}, \code{generateScenarios}
#' @examples
#' # To view the valid variable names use the function without arguments
#' viewModels()
#' 
#' # Examples to view the model options available for different variables
#' viewModels("P")
#' viewModels("Temp")
#' viewModels("Radn")
#' viewModels("PET")
#' @export
viewModels <- function(variable = NULL, compatibleAtts = FALSE) {
  
  # vector of stochastic model tags - exclude scaling
  stochModels <- modelTaglist[-(modelTaglist=="Simple-ann")]
  
  # get stochastic model Info
  varNames <- sapply(strsplit(stochModels, "-"), `[[`, 1)
  parVariation <- sapply(strsplit(stochModels, "-"), `[[`, 2)
  modelType <- sapply(strsplit(stochModels, "-"), `[[`, 3)
  
  # check if the model has additional info
  tagLength <- sapply(strsplit(stochModels, "-"), length)
  for (m in 1:length(modelType)){
    if (tagLength[m] > 3) {
      for (i in 4:tagLength[m]) {
        modelType[m] <- paste0(modelType[m], "-", strsplit(stochModels[m], "-")[[1]][i])
      }
    }
  }
  modelTimeStep <- modelTimeStep[stochModels]
  defaultModel <- rep(FALSE, length(stochModels))
  defaultModel[which(stochModels %in% defaultModelTags)] <- TRUE
  
  if(is.null(variable) || !(variable %in% varNames)) {
    print("Please select a valid variable. The valid variable names are:")
    print(unique(varNames))
  } else {
    # get index
    ind <- which(varNames == variable)
    modelTimeStep <- modelTimeStep[ind]   # subsetted ind here so that the column would be named modelTimeStep in the data.frame output by the function
    defaultModel <- defaultModel[ind]
    
    # model information dataframe
    varModels <- data.frame(modelType[ind], parVariationDef[parVariation[ind]], modelTimeStep, defaultModel)
    colnames(varModels)[1:2] <- mdlFields
    rownames(varModels) <- NULL
    
    if (compatibleAtts) {
      compAtts <- model_attribute_comb[[variable]]
      rownames(compAtts) <- compAtts[ ,1]; compAtts <- t(compAtts[ ,-1])
      if (length(ind) == 1) {
        varModels <-data.frame(c(varModels, compAtts[stochModels[ind], ]))
      } else {
        varModels <- cbind(varModels, compAtts[stochModels[ind], ])
        rownames(varModels) <- NULL
      }
    }
    
    return(varModels)
  }
}

# function to get valid values of modelType and modelParameterVariation; this is an internal function
getModelSpec <- function(colName, var = NULL) {
  modelSpec <- NULL
  
  # If variable is not specified - return data for all variables
  if(is.null(var)){
    var <- fSVars
  }
  for (v in var) {
    modelSpec <- c(modelSpec, as.character(viewModels(v)[[colName]]))
  }
  return(unique(modelSpec))
}

# Anjana: fix the expected length of penaltyWeights in the code - this has to be equal to the length of penaltyAttributes


#################################################################
# NOTE: modelParameterBounds is not fully functional as of now
#     =======================
# To set this, the user would need to know the parameters of the stochastic model and the existing bounds
# So, this is only a developer option as of now - no checks in place - not part of foreSIGHT 1.0 release
# There is an indirect way to know the parameters of the default model options - by using writeControlFile() (with basic = FALSE) to write a sample advanced controlFile.
#################################################################
# Additional option - if the simulation method is "scaling" - why not just say scaling instead of having a json file input
# with simulationMethod = "scaling" - that way we can eliminate the "simulationMethod" field in the namelist input


# Model definifition fields of the namelist
#======================================================================================
# These are the columns used to define a model - printed using output of viewModels()
# The possible values of these fields should always be a string/number (i.e., not a vector)
mdlFields <- c("modelType", "modelParameterVariation") #NOTE: Both modelType & modelParameterVariation have to be specified


mdlBounds <- c("modelParameterBounds")


# Optimization fields - can be specified independent of - mdlFields & mdlBounds
#=======================================================================================
optFields <- c("optimisationArguments", "penaltyAttributes", "penaltyWeights")

# Function to return the masterNamelist options used to check controlFile
getNamelistMaster <- function(){
  
  outList <- list()
  for (f in mdlFields) {
    outList[[f]] <- getModelSpec(f)
  }
  for (f in mdlBounds) {
    outList[[f]] <- list()
  }
  for (f in optFields) {
    if (f == "optimisationArguments") {
      outList[[f]] <- names(optimArgsdefault)[!(names(optimArgsdefault)=="lambda.mult")]
    } else {
      outList[[f]] <- list()
    }
  }
  return(outList)
}

# Function to read controlFile. Will be used in generateScenario
readNamelist <- function(jsonfile) {
  if(!file.exists(jsonfile)) {
    stop("controlFile json file does not exist")
  }
  nml <- jsonlite::fromJSON(txt = jsonfile)
  checkControlFile(nml)
  return(nml)
}

# Function to get the user model choice of "v" from namelist, input is a list - read in from the namelist json file
getModelChoice <- function(nml, v) {
  
  # initialize modelChoice
  modelChoice <- viewModels(v)[FALSE, 1:2]
  modelChoice[1, ] <- NA
  
  # get modelChoice
  for (f in mdlFields) {
    fvalue <- nml[[f]][[v]]
    if (length(fvalue) > 1) {
      stop(paste0("controlFile: ", f, " of variable ", v, " is invalid. Type viewModels(\"", v, "\") to view the valid options"))
    }
    modelChoice[[f]] <- fvalue
  }
  # This a data.frame with columns = mdlFields
  return(modelChoice)
}

# Function to convert namelist modelChoice to modelTag
# Function will return the default modelTag (of corresponding var) if mdlField is not specified in the namelist
# Can be used after reading Obs to create a vector of modelTags
# modelTag = "v-parTag-modelType"
getModelTag <- function(nml, v) {
  modelChoice <- getModelChoice(nml, v)
  
  if (ncol(modelChoice) == 0) {
    # use default
    modelTag <- defaultModelTags[[v]]
  } else {
    # create from nml specifications
    parTagLong <- modelChoice[["modelParameterVariation"]]
    
    # Get the harmonic part of the modelTag based on the variable
    if(parTagLong == "harmonic") {
      parTag <- harDefinition[[v]] 
    } else {
      parTag <- names(parVariationDef[parVariationDef == parTagLong])
    }
    
    modelType <- modelChoice[["modelType"]]
    modelTag <- paste(c(v, parTag, modelType), collapse = "-")
  }
  return(modelTag)
}

# To create the original 'modelInfoMod' list based on the nml input
# 'modelInfoMod' depends only on the vars specified as fields of mdlBounds
getModelInfoMod <- function(nml) {
  
  field <- "modelParameterBounds"
  if(!is.null(nml[[field]])) {
    vars <- names(nml[[field]])
    
    modelInfoMod <- NULL
    for (v in vars){
      modelTag <- getModelTag(nml, v)
      modelInfoMod[[modelTag]] <- modifyParBounds(nml, v)
    }
    return(modelInfoMod)
  } else {
    #cat(paste0("\ncontrolFile does not contain ", field, ". Using foreSIGHT defaults.\n"))
    return(list())
  }
  
}

getOptimArgs <- function(nml) {
  
  field1 <- "optimisationArguments"
  field2 <- "penaltyWeights"
  
  optimArgs <- list()
  
  if(!is.null(nml[[field1]])) {
    optimArgs <- nml[[field1]]
  } else {
    #cat(paste0("\ncontrolFile does not contain ", field1, ". Using foreSIGHT defaults.\n"))
  }
  
  if(!is.null(nml[[field2]])) {
    optimArgs[["lambda.mult"]] <- nml[[field2]]
  }
  
  return(optimArgs)
}

getAttPenalty <- function(nml) {
  
  field <- "penaltyAttributes"
  return(nml[[field]])
  
}



# Given a namelist & variable, return integrated (modelInfo + new bounds) minimum and maximum bounds in a list
modifyParBounds <- function(nml, v) {
  
  field <- "modelParameterBounds"
  boundsIn <- nml[[field]][[v]]
  parsIn <- names(boundsIn)
  
  modelTag <- getModelTag(nml, v)
  parNam <- get.model.info(modelTag)[["parNam"]]
  minBound <- get.model.info(modelTag)[["minBound"]]
  maxBound <- get.model.info(modelTag)[["maxBound"]]
  
  for (p in parsIn) {
    ind <- which(parNam == p)
    minBound[ind] <- boundsIn[[p]][1]
    maxBound[ind] <- boundsIn[[p]][2]
  }
  
  return(list(minBound = minBound,
              maxBound = maxBound))
}


# Function to write the namelist from R list to json file
# Intended use: to write sample 

#' Writes a sample \code{controlFile.json} file
#'
#' \code{writeControlFile()} writes a sample \code{controlFile.json} file. The \code{controlFile.json} file is used to specify alternate model and optimisation options and used as an input to the function \code{generateScenarios}.
#'                          The user may use the sample file created by this function as a guide to create an "\code{controlFile.json}" file for their application.
#' @param jsonfile string; to specify the name of the json file to be written. The default name of the sample file is "sample_controlFile.json".
#'                 The file will be written to the working directory of the user.
#' @param basic logical (\code{TRUE/FALSE}); used to specify whether a "basic" or "advanced" sample file is to be written. The default is \code{TRUE}.
#'               A "basic" controlFile does not contain modelParameterBounds, and is sufficient for most applications.
#' @param nml list; the namelist to be written to the json file, as an R list. This argument may be used to create a JSON file using an controlFile from an existing simulation.
#'            If this argument is set to NULL, the function writes the default model/optimisation options defined in the package to the json file.
#' @return A json file. The file may be used as an example to create an "\code{controlFile.json}" file for input to \code{generateScenarios}. 
#'         An "\code{controlFile.json}" file may contain any subset of the fields listed below. The user may delete the unused fields from the file.
#'         The exception cases where it is mandatory to specify two fields together in controlFile are detailed as part of the list below.
#' \itemize{
#' \item {\strong{\code{modelType}}} {: a list by variable. Each element of the list is a string specifying the type of stochastic model. if \code{modelType} is specified for a variable in controlFile, 
#'                        \code{modelParameterVariation} should also be specified. This is because these two fields together define the stochastic model.
#'                        Use \code{viewModels()} to view the valid options for \code{modelType} by variable.}
#' \item {\strong{\code{modelParameterVariation}}} {: a list by variable. Each element of the list is a string specifying the type of the parameter variation (annual, seasonal, harmonic etc.) of the stochastic model. 
#'                                      if \code{modelParameterVariation} is specified for a variable in controlFile, \code{modelType} should also be specified. 
#'                                      This is because these two fields together define the stochastic model. 
#'                                      Use \code{viewModels()} to view valid options for \code{modelParameterVariation} by variable.}
#' \item {\strong{\code{modelParameterBounds}}} {a nested list by variable. Each element is a list containing the bounds of the parameters of the chosen stochastic model.
#'                                   This field exists to provide an option to overwrite the default bounds of the parameters of the stochastic model.
#'                                   Careful consideration is recommended prior to setting \code{modelParameterBounds} in the controlFile to overwrite the defaults provided in the package.}
#' \item {\strong{\code{optimisationArguments}}} {: a list. Contains the optimisation options used by function \code{ga} from the \code{ga} package. Brief definitions are given below.
#'       \itemize{
#'       \item \code{pcrossover} a value of probability of crossover. Defaults to 0.8.
#'       \item \code{pmutation} a value of probability of mutation. Defaults to 0.1.
#'       \item \code{maxiter} a value of the maximum number of generations. Defaults to 50.
#'       \item \code{maxFitness} a value of the stopping criteria. Defaults to -0.001.
#'       \item \code{popSize} a value of the population size. Defaults to 500.
#'       \item \code{run} a value of an alternative stopping criteria, consecutive runs without improvement in fitness. Defaults to 20.
#'       \item \code{seed} a value of the random seed. Defaults to NULL.
#'       \item \code{parallel} specifies if parallel computing should be used. Defaults to False. Can be set to the number of desired cores, or \code{TRUE}, where it will detect the number of available cores and run.
#'       \item \code{keepBest} specifies if the optimisation should keep the best solution in each generation. Defaults to TRUE.
#'       \item \code{suggestions} suggestions for starting values of parameters for optimisation.
#'       }
#'       }
#' \item {\strong{\code{penaltyAttributes}}} {: a character vector of climate attributes to place specific focus on during targeting via the use of a penalty function during the optimisation process.}
#'                                The \code{penaltyAttributes} should belong to \code{attPerturb} or \code{attHold} that are specified in the exposure space used as input to \code{generateScenarios}.
#'                                If \code{penaltyAttributes} are specified in the controlFile, \code{penaltyWeights} should also be specified.
#' \item {\strong{\code{penaltyWeights}}} {: a numeric vector; the length of the vector should be equal to the length of \code{penaltyAttributes}.
#'                             \code{penaltyWeights} are the multipliers of the corresponding \code{penaltyAttributes} used during the optimisation.}
#' }
#' @details The function may be used without any input arguments to write a "basic" sample controlFile.
#' @seealso \code{generateScenarios}, \code{viewModels}, \code{viewDefaultOptimArgs}
#' @examples
#' \dontrun{
#' # To write a sample controlFile
#' writeControlFile()
#' 
#' # To write an advanced sample controlFile
#' writeControlFile(jsonfile = "sample_controlFile_advanced.json", basic = FALSE)
#' }
#' @export
writeControlFile <- function(jsonfile = "sample_controlFile.json", basic = TRUE, nml = NULL) {
  
  # get defaults
  if (is.null(nml)) {
    modelTag <- defaultModelTags
    names(modelTag) <- NULL
    modelInfo <- get.multi.model.info(modelTag)
    optimArgs <- optimArgsdefault
    attPenalty <- NULL
    nml <- toNamelist(modelTag, modelInfo, optimArgs, attPenalty)
    checkControlFile(nml)
  } else {
    checkControlFile(nml)
  }
  
  if(basic) {
    nml[[mdlBounds]] <- NULL
    nml[["optimisationArguments"]] <- NULL  # removed optimArgs from basic controlFile; add if needed
  }
  json_nml <- jsonlite::toJSON(nml, pretty = TRUE, auto_unbox = TRUE)
  write(json_nml, file = jsonfile)
  
  return(invisible())
  
}

# Function to create a namelist from a number of foreSIGHT options
# A typical use case is after inverse simulation:
#--------------------------------------------------------------------------------
# 1. modelTag - input namelist information + tags for additional variables
# 2. optimArgs - contains input optimArgs merged with default optimArgs
# 3. penaltyAtt - input namelist
# 4. modelInfo - input namelist bounds + exisiting default bounds for unspecified parameters
# The output is of type list - it can be passed to a writeNamelist function to write a json file
toNamelist <- function(modelTag, modelInfoMod = NULL, optimArgs = NULL, attPenalty = NULL) {
  
  if (!is.null(attPenalty)) penaltyWeights <- optimArgs[["lambda.mult"]]
  optimArgs[["lambda.mult"]] <- NULL
  
  # modelParameterVariation
  #----------------------------------------------------
  varNames <- sapply(strsplit(modelTag, "-"), `[[`, 1)
  parTag <- sapply(strsplit(modelTag, "-"), `[[`, 2)
  parVar <- rep(NA, length(parTag))
  for (i in 1:length(parTag)) {
    parVar[i] <- parVariationDef[[parTag[i]]]
  }
  
  # modelType
  #----------------------------------------------------
  mType <- sapply(strsplit(modelTag, "-"), `[[`, 3)
  # check if the model has additional info
  tagLength <- sapply(strsplit(modelTag, "-"), length)
  for (m in 1:length(mType)) {
    if (tagLength[m] > 3) {
      for (i in 4:tagLength[m]) {
        mType[m] <- paste0(mType[m], "-", strsplit(modelTag[m], "-")[[1]][i])
      }
    }
  }
  
  # List with varnames
  #-----------------------------------------------------
  modelParametervariation <- list()
  modelType <- list()
  for (v in 1:length(varNames)) {
    modelParametervariation[[varNames[v]]] <- parVar[v]
    modelType[[varNames[v]]] <- mType[v]
  }
  
  # modelParameterBounds
  #-----------------------------------------------------
  if (!is.null(modelInfoMod) & !is.null(names(modelInfoMod))) {
    modelInfoDef <- get.multi.model.info(modelTag)
    modelInfo <- modifyList(modelInfoDef, modelInfoMod)
    modelParameterBounds <- getModelParBounds(modelTag, modelInfo)
  } else {
    modelParameterBounds <- NULL
  }
  
  # Final Namelist
  #-----------------------------------------------------
  nml <- list()
  nml[["modelType"]] <- modelType
  nml[["modelParameterVariation"]] <- modelParametervariation
  nml[["modelParameterBounds"]] <- modelParameterBounds
  if (!is.null(optimArgs) & !is.null(names(optimArgs))) {
    nml[["optimisationArguments"]] <- optimArgs
    for (n in names(optimArgs)) {
      if (is.null(optimArgs[[n]])) {
        nml[["optimisationArguments"]][[n]] <- NULL
      }
    }
  }
  if(!is.null(attPenalty)) {
    nml[["penaltyAttributes"]] <- attPenalty
    nml[["penaltyWeights"]] <- penaltyWeights
  }
  return(nml)
  
}

# Internal function to get the parameters and their bounds using modelInfo and modelTag
# Used to get information to write namelist out
getModelParBounds <- function(modelTag, modelInfo = NULL) {
  
  varNames <- sapply(strsplit(modelTag, "-"), `[[`, 1)
  
  if (is.null(modelInfo)) {
    # Get the default modelInfo
    modelInfo <- get.multi.model.info(modelTag)
  }
  
  # Anjana: Is modelInfo named using varNames or modelTags? - need to recheck
  modelParameterBounds <- list()
  for (i in 1:length(modelTag)) {
    parNam <- modelInfo[[modelTag[i]]][["parNam"]]
    
    parBounds <- matrix(rep(NA, length(parNam)*2), ncol = length(parNam))
    parBounds[1, ] <- modelInfo[[modelTag[i]]][["minBound"]]
    parBounds[2, ] <- modelInfo[[modelTag[i]]][["maxBound"]]
    
    for (p in 1:length(parNam)) {
      modelParameterBounds[[varNames[i]]][[parNam[p]]] <- parBounds[ ,p]
    }
  }
  return(modelParameterBounds)
  
}

# Function to view the parameter bounds of specific models - intended for user
# Prints the bounds of the parameters
# 1. controlFile is NULL : prints bounds of default model for the variable
# 2. else : prints the bounds of the model specified in the namelist (json file)
# Anjana: use the function 'getModelParBounds' to get the specific function
# viewModelParameters <- function(variable, controlFile = NULL) {
#   
#   #############
#   # TO DO
#   #############
#   
# }

#########################################################
# Functions to check controlFile input by the user
# Will be used inside generateScenario
#########################################################

# Checks on namelist required
# Do the combination of model choices exists as a model Tag?
# Are the names of the model parameters valid for that model Tag? (Can add detailed checks later)
# Is the minimum value of bound lower than the maximum? (Can add detailed checks later)
# Are the optimization argument fields valid (i.e., belongs to optimArgsdefault) (can add detailed checks later)

checkControlFile <- function(nml) {
  
  # Master Namelist
  nmlMaster <- getNamelistMaster()
  
  # Check : names of the fields in input namelist
  namesIn <- names(nml)
  for (n in namesIn) {
    if (!(n %in% names(nmlMaster))) {
      stop(paste0("controlFile field \"", n, "\" unrecognized"))
    }
  }
  
  checkMdlFields(nml)
  checkOptFields(nml)
  checkMdlBounds(nml)
  
}


# check mdlFields
checkMdlFields <- function(nml) {
  
  # Master Namelist
  nmlMaster <- getNamelistMaster()
  
  # Check :  model fields are valid
  for (field in mdlFields) {
    master <- nmlMaster[[field]]
    for (i in unlist(nml[[field]])) {
      if(length(i) > 1 | !is.character(i)) stop(paste0("controlFile ", field, " input", i, " should be a string. Type viewModels() to view the valid options."))
      if(!(i %in% master)) stop(paste0("controlFile ", field, " input", i, "unrecognized"))
    }
  }
  
  # Check : model fields contain the same variables
  vars1 <- names(nml[[mdlFields[1]]])
  for (v in vars1) {
    if (!(v %in% fSVars)) stop(paste0("controlFile variable", v, " specified in ", mdlFields[1], " is unrecognized. ", "Type viewModels() to view the valid variables"))
  }
  
  for (i in 2:length(mdlFields)) {
    varsi <- names(nml[[mdlFields[i]]])
    for (v in varsi) {
      if (!(v %in% fSVars)) stop(paste0("controlFile variable", v, " specified in ", mdlFields[i], " is unrecognized. ", "Type viewModels() to view the valid variables"))
    }
    if (!(all(vars1 %in% varsi) & all(varsi %in% vars1))) stop(paste0("Specify ", paste(mdlFields, collapse="/")," for every variable in controlFile"))
  }
  
  # Check : combination of model choices are okay
  for (v in vars1) {
    modelChoice <- getModelChoice(nml, v)
    # compare with available options
    if (nrow(merge(modelChoice, viewModels(v)[ ,1:2])) != 1) {
      stop(paste0("controlFile: combination of ", paste(mdlFields, collapse = "/"), " specified for variable ", v, " is not valid. ",
                  "Type viewModels(\"", v, "\") to view the valid options"))
    }
  }
  return(invisible())
}


# check mdlBounds
# This function needs to get modelTags from namelist if mdlFields are specified
# Otherwise it needs to get modelTags from the defaults for each variable
checkMdlBounds <- function(nml) {
  
  field <- mdlBounds
  
  vars <- names(nml[[field]])
  for (v in vars) {
    
    if (!(v %in% fSVars)) stop(paste0("controlFile variable ", v, " specified in ", field, " is unrecognized. Type viewModels() to view the valid variables."))
    
    parIn <- names(nml[[field]][[v]])
    modelTag <- getModelTag(nml, v)
    parNam <- get.model.info(modelTag)[["parNam"]]
    
    for (p in parIn) {
      if (!(p %in% parNam)) stop("controlFile parameter ", p,  " specified in ", field, " for variable ", v, " is unrecognized.")
      bounds <- nml[[field]][[v]][[p]]
      if (length(bounds) != 2) stop("controlFile parameter ", p, " specified in ", field, " for variable ", v, " should have minimum and maximum bounds.")
      if(bounds[1] > bounds[2]) stop("controlFile parameter ", p, " specified in ", field, " for variable ", v, ": minimum bound is greater than maximum.")
    }
  }
  return(invisible())
}


# check optFields
checkOptFields <- function(nml) {
  
  # Check names of optimArgs
  field <- "optimisationArguments"
  
  # Master Namelist
  nmlMaster <- getNamelistMaster()
  
  namesIn <- names(nml)
  if(field %in% namesIn) {
    master <- nmlMaster[[field]]
    for (i in names(nml[[field]])) {
      if(!(i %in% master)) stop(paste0("controlFile ", field, " input", i, "unrecognized"))
    }
  }
  
  # Check length of penaltyAtt
  if(! (length(nml[["penaltyWeights"]]) == length(nml[["penaltyWeights"]]))) {
    stop("controlFile: specify penaltyWeights for all penaltyAttributes")
  }
  return(invisible())
}





