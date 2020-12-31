# Note: Have good defaults so that most users won't have to input a controlFile

#' Prints the default optimisation arguments
#'
#' \code{viewDefautOptimArgs()} prints the default values of optimisation arguments (\code{optimisationArguments}) used by \code{generateScenarios}
#' @details The function does not take any input arguments. This a helper function that prints the default values of the optimisation arguments.
#' The user may specify alternate values of these arguments in fields named according to the corresponding argument name nested under 
#' \code{optimisationArguments} in a JSON file to use as the \code{controlFile} input to the \code{generateScenarios} function.
#' @seealso \code{writeControlFile}
#' @examples
#' # To view the default optimisation arguments
#' viewDefaultOptimArgs()
#' @export
viewDefaultOptimArgs <- function() {
  optimArgs_toPrint <- optimArgsdefault
  optimArgs_toPrint[["lambda.mult"]] <- NULL
  optimArgs_toPrint[["suggestions"]] <- NULL
  print(optimArgs_toPrint)
}

optimArgsdefault=list(pcrossover= 0.8,   # list of a parameters used by the ga optimiser (if used)
                      pmutation=0.1,
                      maxiter=50,
                      maxFitness=-0.001,
                      popSize = 500,
                      run=20,
                      seed = NULL,
                      parallel = FALSE,
                      keepBest=TRUE,
                      lambda.mult=0.0,
                      suggestions=NULL)

varShortToLong <- c("P" = "Precipitation",
                   "Temp" = "Temperature",
                   "PET" = "Evapotranspiration",
                   "Radn" = "Radiation")

varUnits <- c("P" = "mm",
              "Temp" = "\u00B0C",
              "PET" = "mm",
              "Radn" = "MJ/m2")

#' Prints the names of and units of valid variables
#'
#' \code{viewVariables()} prints the names of valid variables and their units in the package. The user should input these variable
#' in the same units.
#' @details The function does not take any input arguments.
#' @seealso \code{generateScenarios}
#' @examples
#' # To view the valid variables
#' viewVariables()
#' @export
viewVariables <- function() {
  # vector of stochastic model tags - exclude scaling
  stochModels <- modelTaglist[-(modelTaglist=="Simple-ann")]
  # get variable name
  shortName <- unique(sapply(strsplit(stochModels, "-"), `[[`, 1))
  longName <- varShortToLong[shortName]
  units <- varUnits[shortName]
  outData <- cbind(shortName, longName, units); rownames(outData) <- NULL
  return(outData)
 }



# function to get a vector of varUnits given varNames
getVarUnits <- function(varNames) {
  varUnits <- NA
  for (v in 1:length(varNames)) {
    if (varNames[v] %in% c("Temp")) {
      varUnits[v] <- "\u00B0C"  #expression(paste0(~degree, "C"))
    } else {
      varUnits[v] <- "fraction"
    }
  }
  return(varUnits)
}

# Update the Temp, PET and Radn models to have the har-wgen version, like precipitation
# Temp-har-wgen, PET-har-wgen, Radn-har-wgen, 
modelTaglist=c("Simple-ann",
               "P-ann-wgen",
               "P-seas-wgen",
               "P-har-wgen",
               "Temp-har26-wgen",
               "Temp-har26-wgen-wd",
               "Temp-har26-wgen-wdsd",
               "PET-har26-wgen",
               #"PET-har12-wgen",  # Can I remove this?
               "PET-har26-wgen-wd",
               "Radn-har26-wgen",
               "P-ann-latent",
               "P-har-latent")

# currently, the model time steps are all daily - and the following vector is not used anywhere, 
# but it could be if we add weather generators that use other time steps
# the following data could even be in a format other than a vector for ease of implementation
# assuming that stochastic models of other scale if added would have a separate model tag that includes "-" time step info
modelTimeStep=c("daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily",
                "daily"
                )
names(modelTimeStep) <- modelTaglist

defaultModelTags <- c(P = "P-har-wgen",
                      Temp = "Temp-har26-wgen",
                      PET = "PET-har26-wgen",
                      Radn = "Radn-har26-wgen")

# existing foreSIGHT variables
fSVars <- unique(sapply(strsplit(modelTaglist[-(modelTaglist=="Simple-ann")], "-"), `[[`, 1))

#' Prints the names of valid attributes
#'
#' \code{viewAttributes()} prints the names of valid attributes that may be used to create an exposure space.
#' @details The function does not take any input arguments. The valid attributes that may be specified as \code{attPerturb} or \code{attHold} to create an exposure space using the function \code{createExpSpace} are printed.
#' @seealso \code{createExpSpace}, \code{viewAttributeDef}
#' @examples
#' # To view the valid attributes.
#' viewAttributes()
#' 
#' # To view the definition of any valid attribute
#' viewAttributeDef("P_ann_tot_m")
#' @export
viewAttributes <- function() {
  print(names(attribute.funcs))
}

#' Prints the definition of an attribute
#'
#' \code{viewAttributeDef} prints the short definition of a valid attribute
#' @param attribute A string; the name of the attribute.
#' @seealso \code{viewAttributes}, \code{createExpSpace}
#' @examples
#' # To view the definition of any valid attribute
#' viewAttributeDef("P_ann_tot_m")
#' # To view the valid attributes
#' viewAttributes()
#' @export
viewAttributeDef <- function(attribute) {
  if(attribute %in% names(attribute.funcs)) {
    print(tagBlender(attribute))
  } else {
    print("Please choose a valid attribute: The valid attributes are:")
    viewAttributes()
  }
}

#' Prints the names and bounds of the parameters of the stochastic models
#'
#' \code{viewModelParameters} prints the names of the parameters of the stochastic model and its default minimum and maximum bounds. 
#' The stochastic model is specified using the function arguments.
#' @param variable A string; the name of the variable. Type \code{viewModels()} to view valid variable names
#' @param modelType A string; the model type. Use \code{viewModels} to view the valid values.
#' @param modelParameterVariation A string; the parameter variation. Use \code{viewModels} to view the valid values.
#' @details The available stochastic models can be viewed using the function \code{viewModels()}.
#' This function prints the default ranges of the parameters of the stochastic model specified the
#' stochastic model of interest. 
#' @seealso \code{viewModels}, \code{writeControlFile}
#' @examples
#' viewModelParameters("P", "wgen", "annual")
#' viewModelParameters("P", "wgen", "harmonic")
#' @export
viewModelParameters <- function(variable, modelType, modelParameterVariation) {
  nml <- list()
  nml[["modelType"]] <- list()
  nml[["modelParameterVariation"]] <- list()
  nml[["modelType"]][[variable]] <- modelType
  nml[["modelParameterVariation"]][[variable]] <- modelParameterVariation
  
  modelTag <- getModelTag(nml, variable)
  modelInfo <- get.model.info(modelTag)
  modelPars <- data.frame(parameter = modelInfo[["parNam"]], 
                          min_bound = modelInfo[["minBound"]], 
                          max_bound = modelInfo[["maxBound"]])
  if (nrow(modelPars) < 1) print("Are the input arguments a valid combination of variable|modelType|modelParameterVariation?")
  # colnames(modelPars) <- c("parameter", "minimum bound", "maximum bound")
  print(modelPars)
}


# Anjana: Consider storing modelInfo in sysdata.rda - as part of the model_attribute_comb data.frame
#get.model.info() - based on model tag gets general model information (e.g. nperiods in a year, no. harmonic cycles fitted)
#Get info for individual models
get.model.info<-function(modelTag=NULL #string used to specify model for stochastic generation
){
  
  modelInfo=list()
  #SET UP MODEL RELATED PARAMETERS
  switch(modelTag,
         "Simple-ann"  = {modelInfo$simVar="All"
         modelInfo$simPriority=1 
         modelInfo$nperiod=1
         },
         
         "P-seas-wgen" = {modelInfo$simVar="P"
         modelInfo$simPriority=1
         modelInfo$nperiod=4       # 4 periods in a year
         modelInfo$fixedPars=NA    # No fixed pars
         modelInfo$ncycle=NA       # No harmonic fit
         modelInfo$npars=modelInfo$nperiod*4 #par vector is of length 16
         modelInfo$parNam=c("pdd_1","pdd_2","pdd_3","pdd_4",
                            "pwd_1","pwd_2","pwd_3","pwd_4",
                            "alpha_1","alpha_2","alpha_3","alpha_4",
                            "beta_1","beta_2","beta_3","beta_4")
         modelInfo$minBound=c(0.389, 0.334, 0.375, 0.277, 0.078, 0.079, 0.084, 0.036, 
                              0.295,	0.303, 0.309,	0.257, 0.043,	0.046, 0.048, 0.034) #Aus 3stdev hard bounds
         modelInfo$maxBound=c(0.997, 0.989, 0.994, 0.998, 0.85, 0.714, 0.714, 0.808, 
                              0.998, 0.998, 0.998, 0.998, 15.716, 30.08, 27.877, 21.193)
         #bounds here?????????????
         #npar.optim???? - then split into max, min bounds
         },
         "P-ann-wgen" = {modelInfo$simVar="P"
         modelInfo$simPriority=1
         modelInfo$nperiod=1
         modelInfo$fixedPars=NA
         modelInfo$ncycle=NA
         modelInfo$npars=modelInfo$nperiod*4       #par vector is of length 4
         modelInfo$parNam=c("pdd","pwd","alpha","beta")
         modelInfo$minBound=c(0.427, 0.088, 0.313, 0.043) #Aus 3stdev hard bounds
         modelInfo$maxBound=c(0.998, 0.824, 0.998, 25.46)
         },
         "P-har-wgen" = {modelInfo$simVar="P"
         modelInfo$simPriority=1
         modelInfo$nperiod=365
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=12  #par vector is of length 12 
         #modelInfo$parNam=c("pdd", "pwd", "alpha", "beta")
         modelInfo$parNam=c("pdd_m","pdd_amp","pdd_ang",
                            "pwd_m","pwd_amp","pwd_ang",
                            "alpha_m","alpha_amp","alpha_ang",
                            "beta_m","beta_amp","beta_ang")
         
         # modelInfo$minBound=c(0.476, 0.006, 0.730,
         #                      0.093, 0.004, 0.543,
         #                      0.33, 0.002, 4.108,
         #                      0.085, 0.028, 1.348) #Aus 3stdev hard bounds
         # modelInfo$maxBound=c(0.950, 0.257, 0.733,
         #                      0.728, 0.319, 0.545,
         #                      0.950, 0.200, 4.110,
         #                      15.00, 6.50, 1.350)
         modelInfo$minBound=c(0.476, 0.006, 0,
                              0.093, 0.004, 0,
                              0.33, 0.002, 0,
                              0.085, 0.028, 0) # Culley 2019 I have widened the 3stdev bounds above, to allow them to be manually controlled with user input for faster testing (no re-compiling needed)
         modelInfo$maxBound=c(0.950, 0.557, 6.28,
                              0.728, 0.519, 6.28,
                              0.950, 0.600, 6.28,
                              15.00, 10, 6.28)
         },
         "P-ann-latent" = {modelInfo$simVar="P"
         modelInfo$simPriority=1
         modelInfo$nperiod=1
         modelInfo$fixedPars=NA
         modelInfo$ncycle=NA
         modelInfo$npars=modelInfo$nperiod*4
         modelInfo$parNam=c("alpha", "sigma", "mu", "lambda")
         modelInfo$minBound=c(0, 0.001, -15, 1)
         modelInfo$maxBound=c(0.999, 10, 0, 2)
         },
         "P-har-latent" = {modelInfo$simVar="P"
         modelInfo$simPriority=1
         modelInfo$nperiod=365
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=12  #par vector is of length 12 since each par has a mean, amplitude & phase angle
         #modelInfo$parNam=c("alpha", "sigma", "mu", "lambda")
         modelInfo$parNam=c("alpha_m","alpha_amp","alpha_ang",
                            "sigma_m","sigma_amp","sigma_ang",
                            "mu_m","mu_amp","mu_ang",
                            "lambda_m","lambda_amp","lambda_ang")
         modelInfo$minBound=c(0, 0, 0,
                              0.001, 0, 0,
                              -15, 0, 0,
                              1, 0, 0)
         modelInfo$maxBound=c(0.999, 0, 0,
                              10, 5, 6.28,
                              0, 8, 6.28,
                              2, 0, 0)
         },
         "Temp-har26-wgen-wd" = {modelInfo$simVar="Temp"
         modelInfo$simPriority=2
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=TRUE  #conditioned on wet/dry status
         modelInfo$wdCycle="All"
         modelInfo$nperiod=26
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=4*(1+modelInfo$ncycle*2)+1             #par vector is of length  13
         modelInfo$parNam=c("cor0",
                            "W-mCycle-m","W-mCycle-amp","W-mCycle-ang",
                            "W-sCycle-m","W-sCycle-amp","W-sCycle-ang",
                            "D-mCycle-m","D-mCycle-amp","D-mCycle-ang",
                            "D-sCycle-m","D-sCycle-amp","D-sCycle-ang")
         
         modelInfo$minBound=c(0.45,7.0,1.0,-0.05,0.9,0.1,-1.6,7.0,1.0,-0.05,0.9,0.1,-1.6) #Placeholder bounds
         modelInfo$maxBound=c(0.90,28.0,9.0,0.81,4.9,1.4,3.15,28.0,9.0,0.81,4.9,1.4,3.15)
         },
         "Temp-har26-wgen" = {modelInfo$simVar="Temp"
         modelInfo$simPriority=2
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=FALSE  #conditioned on wet/dry status
         modelInfo$wdCycle=FALSE
         modelInfo$nperiod=26
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=2*(1+modelInfo$ncycle*2)+1             #par vector is of length  7
         modelInfo$parNam=c("cor0",
                            "WD-mCycle-m","WD-mCycle-amp","WD-mCycle-ang",
                            "WD-sCycle-m","WD-sCycle-amp","WD-sCycle-ang")
         modelInfo$minBound=c(0.45,7.0,1.0,-0.05,0.9,0.1,-1.6) #Placeholder bounds
         modelInfo$maxBound=c(0.9,28.0,9.0,0.81,4.9,1.4,3.15)
         },
         "Temp-har26-wgen-wdsd" = {modelInfo$simVar="Temp"
         modelInfo$simPriority=2
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=TRUE  #conditioned on wet/dry status
         modelInfo$wdCycle="sCycle"
         modelInfo$nperiod=26
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=3*(1+modelInfo$ncycle*2)+1             #par vector is of length  10
         modelInfo$parNam=c("cor0",
                            "WD-mCycle-m","WD-mCycle-amp","WD-mCycle-ang",
                            "W-sCycle-m","W-sCycle-amp","W-sCycle-ang",
                            "D-sCycle-m","D-sCycle-amp","D-sCycle-ang")
         modelInfo$minBound=c(0.45,7.0,1.0,-0.05,0.9,0.1,-1.6,0.9,0.1,-1.6) #aus bounds
         modelInfo$maxBound=c(0.90,28.0,9.0,0.81,4.9,1.4,3.15,4.9,1.4,3.15)
         },
         "PET-har12-wgen" = {modelInfo$simVar="PET"
         modelInfo$simPriority=2
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=FALSE  #conditioned on wet/dry status
         modelInfo$wdCycle=FALSE
         modelInfo$nperiod=12
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=2*(1+modelInfo$ncycle*2)+1             #par vector is of length  7
         modelInfo$parNam=c("cor0",
                            "WD-mCycle-m","WD-mCycle-amp","WD-mCycle-ang",
                            "WD-sCycle-m","WD-sCycle-amp","WD-sCycle-ang")
         modelInfo$minBound=c(0.0, 
                              0,0.01,0.2,
                              0.01,0.01,0.2)  #NB: Placeholder bounds
         modelInfo$maxBound=c(0.9,
                              6,5,0.3,
                              3,3,0.3)
         },
         "PET-har26-wgen" = {modelInfo$simVar="PET"
         modelInfo$simPriority=2
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=FALSE  #conditioned on wet/dry status
         modelInfo$wdCycle=FALSE
         modelInfo$nperiod=26
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=2*(1+modelInfo$ncycle*2)+1             #par vector is of length  7
         modelInfo$parNam=c("cor0",
                            "WD-mCycle-m","WD-mCycle-amp","WD-mCycle-ang",
                            "WD-sCycle-m","WD-sCycle-amp","WD-sCycle-ang")
         modelInfo$minBound=c(0.0  ,0  ,0.01 ,0.2 ,1 ,0.4 ,0.2)  #NB: Placeholder bounds
         modelInfo$maxBound=c(0.9  ,6  ,5    ,0.3 ,3    ,3    ,0.3)
         },
         "PET-har26-wgen-wd" = {modelInfo$simVar="PET"
         modelInfo$simPriority=2
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=TRUE  #conditioned on wet/dry status
         modelInfo$wdCycle="All"
         modelInfo$nperiod=26
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=4*(1+modelInfo$ncycle*2)+1             #par vector is of length  13
         modelInfo$parNam=c("cor0",
                            "W-mCycle-m","W-mCycle-amp","W-mCycle-ang",
                            "W-sCycle-m","W-sCycle-amp","W-sCycle-ang",
                            "D-mCycle-m","D-mCycle-amp","D-mCycle-ang",
                            "D-sCycle-m","D-sCycle-amp","D-sCycle-ang")
         
         modelInfo$minBound=c(0.001,
                              0.01,0.01,0.95,
                              0.01,0.01,0.9,
                              0.01,0.01,0.95,
                              0.01,0.01,0.9) #Placeholder bounds
         modelInfo$maxBound=c(0.95,
                              30.0,10.0,1.1,
                              10.0,10.0,1.05,
                              30.0,9.0,1.1,
                              10.0,10.0,1.05)
         },
         "Radn-har26-wgen" = {modelInfo$simVar="Radn"
         modelInfo$simPriority=3
         modelInfo$nAssocSeries=0
         modelInfo$WDcondition=FALSE  #conditioned on wet/dry status
         modelInfo$wdCycle=FALSE
         modelInfo$nperiod=26
         modelInfo$fixedPars=NA
         modelInfo$ncycle=1 
         modelInfo$npars=2*(1+modelInfo$ncycle*2)+1             #par vector is of length  7
         modelInfo$parNam=c("cor0",
                            "WD-mCycle-m","WD-mCycle-amp","WD-mCycle-ang",
                            "WD-sCycle-m","WD-sCycle-amp","WD-sCycle-ang")
         modelInfo$minBound=c(0.45,7.0,1.0,-0.05,0.9,0.1,-1.6) #Placeholder bounds
         modelInfo$maxBound=c(0.9,29.0,10.0,0.81,4.9,1.4,3.15)
         },
         #--- MORE VERSIONS COMING ---
         
         # "P-2har26-wgen-FS" = {modelInfo$simVar="P"
         # modelInfo$nperiod=26
         #                       modelInfo$fixedPars="phase.angle"
         #                       modelInfo$ncycle=2
         #                       modelInfo$npars=4*(1+modelInfo$ncycle*1)  #par vector is of length 12
         # },
         # versions where occurence w/d is kept the same as current
         
         -999
  )
  return(modelInfo)
  
}

get.attribute.info <- function(modelTag = NULL){

  if (modelTag == "Simple-ann") {
    validAttributes <- c("P_ann_tot_m", "Temp_ann_avg_m", "PET_ann_avg_m", "Radn_ann_avg_m")
  } else {

    # Identify variable
    varType <- get.varType(modelTag, sep = "-")

    # Use the appropriate data.frame to read the validAttributes
    valid_indices <- which(model_attribute_comb[[varType]][[modelTag]] == TRUE)
    validAttributes <- model_attribute_comb[[varType]]$validAttributes[valid_indices]

  }
  return(validAttributes)
}



