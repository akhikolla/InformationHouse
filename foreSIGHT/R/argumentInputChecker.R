################################
##  CHECKS FOR ARGUMENT TAGS  ##
################################

# CONTAINS
  # get.varType() - grabs first part of string, sep specifiable  
  # argument_check()- check duplicate attribute/model tags + check requests for two model types of one variable + checking other master control arguments
  
#---------------------------------------------------------------------------
#FUNCTIONS

#Function to split string and extract first component
get.varType<-function(attrib=NULL, # attribute name
                      sep="_"){
  varType=strsplit(x = attrib,split=sep)[[1]][1]
  return(varType)
}

# New checks for the createExpSpace function
check_attributes <- function(attHold = NULL,
                             attPerturb = NULL,
                             attTargetsFile = NULL,
                             attPerturbSamp = NULL,
                             attPerturbBy = NULL,
                             attPerturbMin = NULL,
                             attPerturbMax = NULL,
                             attribute.funcs = NULL
) {
  
  attributelist <- names(attribute.funcs)

  # 1. Checks on the names of attributes specified
  #------------------------------------------------------------------
  attSel <- c(attPerturb, attHold)
  if(is.null(attPerturb)){
    stop("No attributes nominated for perturbation")
  }
  
  if(is.null(attHold)){
    message("Note: There are no attributes held at historical levels")
  }
  
  if (anyDuplicated(attSel)!=0) {
    stop("There are multiple entries of the same attribute")
  }
  
  if(!is.null(attHold)){
    for (i in 1:length(attHold)){
      if(sum(attHold[i] %in% attributelist)==0){
        stop(paste0("attHold [",i,"] unrecognised"))
      }
    }
  }
  
  for (i in 1:length(attPerturb)){
    if(sum(attPerturb[i] %in% attributelist)==0){
      stop(paste0("attPerturb [",i,"] unrecognised"))
    }
  }
  
  # 2. Checks on arguments used for creating the sample space
  #------------------------------------------------------------------
  if (!is.null(attTargetsFile)) {
    if (!is.character(attTargetsFile)) { stop("attTargetsFile should be the path of the csv file with targets")}
  }
  
  if(is.character(attTargetsFile)) {
    # READING FROM FILE
    targetMat <- read.table(file = attTargetsFile, sep = ",", header = TRUE)
    att_frmFile <- names(targetMat)

    for (i in 1:length(att_frmFile)) {
      if(sum(att_frmFile[i] %in% attSel)==0){
        stop("There is a mismatch in attributes specified in attPerturb & attHold and attTargetsFile")
      }
    }
    
    if (length(att_frmFile) != length(attSel)) {
      stop("Ensure that targets for attPerturb & attHold are specified in attTargetsFile")
    }

  } else {
    
    if (!is.null(attPerturbSamp)) {
      if (length(attPerturb) != length(attPerturbSamp)) {
        stop("attPerturbSamp should be specified for each attribute in attPerturb")
      }
      if (!is.numeric(attPerturbSamp)) stop("Enter numeric values for attPerturbSamp")
    }
    
    if (!is.null(attPerturbSamp)) {
      if (length(attPerturb) != length(attPerturbSamp)) {
        stop("attPerturbSamp should be specified for each attribute in attPerturb")
      }
      if (!is.numeric(attPerturbSamp)) stop("Enter numeric values for attPerturbSamp")
      if (any(attPerturbSamp < 0) | !all((attPerturbSamp %% 1) == 0)) stop("Enter positive integers for attPerturbSamp") 
    }
    
    if (!is.null(attPerturbBy)) {
      if (length(attPerturb) != length(attPerturbBy)) {
        stop("attPerturbby should be specified for each attribute in attPerturb")
      }
      if (!is.numeric(attPerturbBy)) stop("Enter numeric values for attPerturbBy")
      if (any(attPerturbBy < 0)) stop("Enter positive values for attPerturbBy") 
    }
    
    if (length(attPerturb) != length(attPerturbMin)) {
      stop("attPerturbMin should be specified for each attribute in attPerturb")
    }
    
    if (length(attPerturb) != length(attPerturbMax)) {
      stop("attPerturbMax should be specified for each attribute in attPerturb")
    }
    
    if (!is.numeric(attPerturbMin)) stop("Enter numeric values for attPerturbMin")
    if (!is.numeric(attPerturbMax)) stop("Enter numeric values for attPerturbMax")
    
    if (!all(attPerturbMin <= attPerturbMax)) {
      stop("attPerturbMin should be less than or equal to attPerturbMax")
    }
    
  }
  
  return(invisible(NULL))
  
}

# Function to check supplied arguments
# Anjana: Revisit argument checks to create a small common function to make the checking if-conditions more compact
check_duplicates_mismatch<-function(obs=NULL,
                         attSel=NULL,        
                         attPrim=NULL,
                         attHold=NULL,
                         attPerturb=NULL,
                         modelTag=NULL,
                         optimArgs=NULL,
                         file
                         ){
 
  # variables in the input data
  names <- names(obs)
  names<-names[names!="year"];names<-names[names!="month"];names<-names[names!="day"]
  
  # Anjana - commented after createExpSpace
  # # Perturbed attributes should exist
  # if(is.null(attPerturb)){
  #   logfile("Error: No attributes nominated for perturbation",file)
  #   logfile("Program terminated",file)
  #   stop("No attributes nominated for perturbation")
  # }
  
  # Simple scaling : no attHeld, no attPrim, single perturbed attribute per variable
  if (modelTag[1]=="Simple-ann") {
    
    if(length(attHold)!=0) {
      logfile("Error: Invalid - Simple scaling cannot hold attributes constant",file)
      logfile("Program terminated",file)
      stop("Simple scaling cannot hold attributes constant")
    }
    
    if(length(attPrim)!=0) {
      logfile("Error: Simple scaling uses no primary attributes",file)
      logfile("Program terminated",file)
      stop("Simple scaling uses no primary attributes")
    }
    
    if(length(attPerturb)!=length(names)){
      logfile("Error: There is a mismatch between number of variables and number of attributes. These should be the same for simple scaling, which only has multiplicative or additive changes",file)
      logfile("Program terminated",file)
      stop("There is a mismatch between number of variables and number of attributes. These should be the same for simple scaling, which only has multiplicative or additive changes")
    }
  
    
  # Checks for stochastic models  
  } else {
    
    if(is.null(attHold)){
      warn("No attributes held at historical levels",file)
    }
    
    # Anjana - commented after createExpSpace
    # # CHECK FOR DUPLICATE TAGS
    # if (anyDuplicated(attSel)!=0) {
    #   logfile("Error:There are multiple entries of the same attribute",file)
    #   logfile("Program terminated",file)
    #   stop("There are multiple entries of the same attribute")
    # }
    
    if (anyDuplicated(attPrim)!=0) {
      logfile("Error: There are multiple entries of the same primary attribute",file)
      logfile("Program terminated",file)
      stop("There are multiple entries of the same primary attribute")
    }
    
    # Check that modelTag and attribute names are recognized
    for(i in 1:length(modelTag)){
      if(sum(modelTag[i] %in% modelTaglist)==0){
        logfile("Error: modelTag unrecognised",file)
        logfile("Program terminated",file)
        stop(paste0("modelTag ",i," unrecognised"))
      }
    }
    
    # Anjana - commented after createExpSpace
    # if(!is.null(attHold)){
    #   for (i in 1:length(attHold)){
    #     if(sum(attHold[i] %in% attributelist)==0){
    #       logfile("Error: attHold unrecognised",file)
    #       logfile("Program terminated",file)
    #       stop(paste0("attHold ",i," unrecognised"))
    #     }
    #   }
    # }
    
    #CHECKS FOR TWO REQUESTED MODEL TYPES
    modelVars<-sapply(modelTag,get.varType,USE.NAMES=FALSE,sep="-")
    
    if (anyDuplicated(modelVars)!=0) {
      logfile("Error: There are multiple entries of a model type for one variable",file)
      logfile("Program terminated",file)
      stop("There are multiple entries of a model type for one variable")
    }
    
    # Checks for columns of data without model tags.
    if (length(which((names %in% modelVars)==FALSE))>0) {
      message("reference contains more variables than the specified attributes or models. Stochastic series will only be produced for the specified settings.")
      #warn("There is a mismatch between provided model types and supplied variables. Stochastic series will only be produced for supplied model tags",file)
      array<-c("year","month","day",modelVars)
      obs=obs[array]
    }
    
  }
  
  # Anjana - commented after createExpSpace 
  # # Checks common to Simple scaling and stochastic models
  # for (i in 1:length(attPerturb)){
  #   if(sum(attPerturb[i] %in% attributelist)==0){
  #     logfile("Error: attPerturb unrecognised",file)
  #     logfile("Program terminated",file)
  #     stop(paste0("attPerturb ",i," unrecognised"))
  #   }
  # }
  
  if (anyDuplicated(modelTag)!=0) {
    logfile("Error: There are multiple entries of the same model tag",file)
    logfile("Program terminated",file)
    stop("There are multiple entries of the same model tag")
  }
  
  # Other checks for Stochastic models
  # CHECKS FOR LAMBDA VALUES
  
  # attPrim should exist in attSel
  if(!is.null(attPrim)){
    for (i in 1:length(attPrim)){
      if(sum(attPrim[i] %in% attSel)==0){
        logfile(paste0("contolFile: penaltyAttribute [",i,"] does not exist in the expSpace"), file)
        logfile("Program terminated",file)
        stop(paste0("contolFile: penaltyAttribute [",i,"] does not exist in the expSpace"))
      }
    }
  }
  
  if((length(attPrim!=0)) & (length(attPrim)!=length(which(optimArgs$lambda.mult>0)))) {
    warn("contolFile: There are specified penaltyAttributes with a lambda value of zero",file)
  }
  
  if(length(attPrim)>length(optimArgs$lambda.mult)){        # NO. OF ATTPRIM IS GREATER THAN LAMBDA VECTOR
    warn("There are more specified penaltyAttributes than lambda values",file)
    logfile("Error: check number of supplied lambda values",file)
    logfile("Program terminated",file)
    stop("Ensure a lambda value is entered for each Primary attribute")
  }else{
    note=paste0("Lambda(",attPrim,"): ",optimArgs$lambda.mult,collapse = ", ")
    progress(note,file)
    logfile(note,file)
  }
  
  # Anjana: Commented as I think these checks are not required after moving createExpSpace to separate function
  # # Add more checks for ExSpArgs - ?
  # 
  # if(is.character(exSpArgs)==TRUE){
  #   # READING FROM FILE (ASSUMED SELF-CHECKED)
  #   # AMEND THIS
  # }else {
  #   
  #   # CHECKS FOR BOUNDS
  #   boundVars<-sapply(names(exSpArgs$bounds),get.varType,USE.NAMES=FALSE,sep="_")
  #   attVars<-sapply(attSel,get.varType,USE.NAMES = FALSE)
  #   boundNames<-names(exSpArgs$bounds)
  #   
  #   # Anjana: There is a problem with these checks since it assumes that the ORDER of attributes (names) in exSpArgs are 
  #   # same as that in attPerturb which need not be the case - revisit
  #   
  #   if (modelTag[1]=="Simple-ann") {
  #     if (!isTRUE(all(boundVars==attVars))) {
  #       logfile("Error: Ensure bounds are entered for each variable in provided data",file)
  #       logfile("Program terminated",file)
  #       stop("Ensure bounds are entered for each variable in provided data")
  #     }
  #   } else {
  #     # switch to attPerturb
  #     if (!isTRUE(all(boundNames==attSel))) { 
  #       logfile("Error: Enter bounds for each attribute in attPerturb",file)
  #       logfile("Program terminated",file)
  #       stop("Enter bounds for each attribute in attPerturb")
  #     }
  #   } 
  #   
  #   }
  
  return(invisible())
  
}

#############################################
##  LOGIC CHECKS FOR ATTRIBUTE/MODEL TAGS  ##
#############################################

#CONTAINS
# argument_logic_check()
#checks for simple scaling attributes
#checks for matching variable, attribute and model lists

check_models_attributes<-function(names=NULL,
                                  attSel=NULL,        
                                  attPrim=NULL,
                                  modelTag=NULL,
                                  file
){
  
  nam<-names[-c(1:3)]
  
  modelVars<-sapply(modelTag,get.varType,USE.NAMES=FALSE,sep="-")
  
  # Can't do rain dependent Temp or PET with no P
  if (sum("P" %in% modelVars)==0) {
    if(sum("Temp-har26-wgen-wd" %in% modelTag)==1) {
      logfile("Error: Cannot simulate stochastic wet/dry dependent temperature without a rainfall model",file)
      stop("Cannot simulate stochastic wet/dry dependent temperature without a rainfall model")
    } else if(sum("PET-har26-wgen-wd" %in% modelTag)==1) {
      logfile("Error: Cannot simulate stochastic wet/dry dependent PET without a rainfall model",file)
      stop("Cannot simulate stochastic wet/dry dependent PET without a rainfall model")
    }
    
  }
  
  if (modelTag[1]=="Simple-ann") {
    
    validAtts <- get.attribute.info(modelTag = "Simple-ann")
    if(sum(attSel %in% validAtts)!=length(attSel)) {
      logfile("Error: Simple scaling cannot perturb selected attributes",file)
      logfile("Program terminated",file)
      stop("Simple scaling cannot perturb selected attributes. Choose a stochastic model")
    }
    
  } else {
    
    validAtts=("temp")
    for(i in 1:length(modelVars)) {
      temp=get.attribute.info(modelTag=modelTag[i])
      validAtts=append(validAtts,temp)
    }
    validAtts=validAtts[-1]
    
    if(sum(attSel %in% validAtts)!=length(attSel)) {
      logfile("Error: Model combinations cannot perturb selected attributes",file)
      logfile("Program terminated",file)
      stop("Model combinations cannot perturb or hold selected attributes. Change attPerturb or attHold selection.")
    }
    
    progress("You have selected the following penalty attributes:",file)
    # cat("     ")
    # cat(attPrim,sep=", ")
    # cat("\n")
    # cat("\n")
    logfile(attPrim,file)
    progress("These attributes will be perturbed with model types:",file)
    # cat("     ")
    # cat(modelTag,sep=", ")
    # cat("\n")
    # cat("\n")
    logfile(modelTag,file)
    progress("The scenarios will include the following attributes in the objective function:",file)
    # cat("     ")
    # cat(attSel,sep=", ")
    # cat("\n")
    # cat("\n")
    logfile(attSel,file)
  
  }
  
  return(invisible())
  
  #model assessor
  
  
}
