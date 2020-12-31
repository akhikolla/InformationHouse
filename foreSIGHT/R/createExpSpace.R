
# Should we use a namelist file to specify the attributes and samples?
# Note: Both attPerturb and attHold have to specified in the target file

#' Creates exposure space of hydroclimatic targets for generation of scenarios using 'generateScenarios'
#'
#' \code{createExpSpace} returns a list containing the targets (\code{targetMat}) and the metadata (input arguments) used to create the exposure space. 
#'
#' @param attPerturb A char vector; the names of the attributes to be perturbed. 
#'                   This vector can contain attributes of different hydroclimatic variables.
#' @param attPerturbSamp An integer vector; the number of samples for each attribute \code{attPerturb}. 
#'                       The length of this vector should be equal to the length of \code{attPerturb}.
#' @param attPerturbMin A numeric vector; the minimum bounds for sampling of \code{attPerturb}. 
#'                      The length of this vector should be equal to the length of \code{attPerturb}. 
#'                      For variables like precipitation, evapotranspiration, radiation, etc. \code{attPerturbMin} should specified as a fraction of the original 
#'                      (eg: 0.9 = 90\% of the original attribute). For temperature, \code{attPerturbMin} should be specified in K (eg: 0.9 = 0.9 K).
#' @param attPerturbMax A numeric vector; the maximum bounds for sampling of \code{attPerturb}. 
#'                      The length of this vector should be equal to the length of \code{attPerturb}. 
#'                      For variables like precipitation, evapotranspiration, radiation, etc. \code{attPerturbMax} should specified as a fraction of the original 
#'                      (eg: 0.9 = 90\% of the original attribute). For temperature, \code{attPerturbMax} should be specified in K (eg: 0.9 = 0.9 K).
#'                      Note that to create a single sample of the attribute, \code{attPerturbSamp} could be specified as 1 with \code{attPerturbMin} and \code{attPerturbMax} specified as equal. 
#' @param attPerturbType A string to specify the type of sampling, defaults to regular spacing. Valid sampling types are:
#' \itemize{
#' \item "regGrid" a regular grid sampling all the attributes specified in \code{attPerturb} simultaneously
#' \item "OAT" one-at-a-time sampling of the attributes specified in \code{attPerturb}
#' }
#' @param attPerturbBy A numeric vector; increment of values to create samples between \code{attPerturbMin} and \code{attPerturbMax}.
#' If \code{attPerturbBy} is specified, attPerturbSamp should be set as \code{NULL}.
#' @param attHold A char vector; the names of the attributes to be held at historical levels. 
#'                This vector can contain attributes of different hydroclimatic variables.
#' @param attTargetsFile String specifying the full path to a CSV file containing the target exposure space. 
#'                       The column names in the file should correspond to the attributes specified in \code{attPerturb} and \code{attHold}. 
#'                       \code{attTargetsFile} is alternate way to specify exposure space targets that do not form a regular grid. 
#'                       If \code{attTargetsFile} is specified, the inputs arguments \code{attPerturbSamp}, \code{attPerturbMin}, \code{attPerturbMax}, 
#'                       and \code{attPerturbType} should be set to \code{NULL} and will not be used by the function.
#' @return The exposure space as a list containing the following fields:
#' \itemize{
#' \item \code{targetMat} a dataframe or matrix; each column is a perturb/hold attribute, each row is a point in the exposure space.
#' \item \code{attRot} a char vector containing the one-at-a-time ("OAT") attributes associated with \code{targetMat}, \code{attRot} is \code{NULL} for other types of sampling.
#' \item \code{attPerturb}, \code{attHold}, \code{attPerturbSamp}, \code{attPerturbMin}, \code{attPerturbMax}, \code{attPerturbType} in the  function input arguments, if not \code{NULL}.
#' }
#' @details The list of valid attributes that may be specified using \code{attPerturb} or \code{attHold} can be viewed using the function
#' \code{viewAttributes()}. The definition of the attribute can be viewed using the function \code{viewAttributeDef}.
#' @seealso \code{generateScenarios}, \code{viewAttributes}, \code{viewAttributeDef}
#' @examples
#' # To view the valid attributes. The function does not take any input arguments. 
#' viewAttributes()
#' 
#' # To view the definition of any valid attribute
#' viewAttributeDef("P_ann_tot_m")
#' 
#' # To create an exposure space of points on a regular grid
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m")
#' attPerturbType <- "regGrid"
#' attPerturbSamp <- c(3, 1, 1)
#' attPerturbMin <- c(0.9, 1, 1)
#' attPerturbMax <- c(1.1, 1, 1)
#' attHold <- c("P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", 
#' "P_MAM_tot_m", "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
#' expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = attPerturbSamp, 
#' attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, 
#' attPerturbType = attPerturbType, attHold = attHold, attTargetsFile = NULL)
#' 
#' # Using attPerturbBy to specify the increment of perturbation (attPerturbSamp set to NULL)
#' 
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m")
#' attPerturbType <- "regGrid"
#' attPerturbMin <- c(0.9, 1, 1)
#' attPerturbMax <- c(1.1, 1, 1)
#' attPerturbBy <- c(0.1, 0, 0)
#' attHold <- c("P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m", 
#' "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
#' expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = NULL, 
#' attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, attPerturbType = attPerturbType, 
#' attPerturbBy = attPerturbBy, attHold = attHold, attTargetsFile = NULL)
#' 
#' # To create an exposure space of observed attributes without perturbation
#' # Note that attPerturbMin and attPerturbMax values are set to 1 for variables like precipitation, 
#' # and 0 for temperature 
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m", "Temp_DJF_avg_m")
#' attPerturbType <- "regGrid"
#' attPerturbSamp <- c(1, 1, 1, 1)
#' attPerturbMin <- c(1, 1, 1, 0)
#' attPerturbMax <- c(1, 1, 1, 0)
#' expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = attPerturbSamp, 
#' attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, attPerturbType = attPerturbType, 
#' attHold = NULL, attTargetsFile = NULL)
#' @export

createExpSpace <- function(attPerturb,
                           attPerturbSamp,
                           attPerturbMin,
                           attPerturbMax,
                           attPerturbType = "regGrid",
                           attPerturbBy = NULL,
                           attHold = NULL,
                           attTargetsFile = NULL # If this file is specified, use this, else create based on sample space
) {
  
  # print("CHECKING INPUT ARGUMENTS")
  
  if (!is.null(attPerturbBy)) {
    if (!is.null(attPerturbSamp)) {
      stop("Since attPerturbBy is specified, attPerturbSamp should be set to NULL")
    }
  } else {
    if (is.null(attTargetsFile)) {
      if (is.null(attPerturbSamp)) {
        stop("Either attPerturbSamp or attPerturbBy should be specified.")
      }
    }
  }
  
  check_attributes(attHold = attHold,
                   attPerturb = attPerturb,
                   attTargetsFile = attTargetsFile,
                   attPerturbSamp = attPerturbSamp,
                   attPerturbBy = attPerturbBy,
                   attPerturbMin = attPerturbMin,
                   attPerturbMax = attPerturbMax,
                   attribute.funcs = attribute.funcs)
  
  # convert "by" to "samp" since downstream code is set up to work with samp
  if (!is.null(attPerturbBy)) {
    attPerturbSamp <- floor((attPerturbMax - attPerturbMin)/attPerturbBy) + 1
    attPerturbSamp[attPerturbBy == 0] <- 1
    attPerturbMax <- attPerturbMin + (attPerturbSamp-1)*attPerturbBy
  }
  
  # create exSpArgs
  if (!is.null(attTargetsFile)) {
    exSpArgs <- attTargetsFile
    } else {
      exSpArgs <- list()
      exSpArgs$samp <- attPerturbSamp
      exSpArgs$type <- attPerturbType
      exSpArgs$bounds <- createBounds(attPerturb, attPerturbMin, attPerturbMax)
      exSpArgs <- addExpArgs_attHold(attPerturb, attHold, exSpArgs)
    }
  
  attSel <- c(attPerturb, attHold)
  # create attInfo
  attInfo=list()
  #attribute name chopper function 
  attInfo$varType=vapply(attSel,FUN = get.attribute.varType,FUN.VALUE=character(1),USE.NAMES = FALSE) #drop use of names as comes ordered anyway
  #ASSIGN TARGET TYPE (IF P USE "FRAC", IF T USE "DIFF")
  attInfo$targetType=vapply(attInfo$varType,FUN=get.target.type,FUN.VALUE=character(1),USE.NAMES=FALSE)
  
  # create temporary log file
  file <- paste0(tempdir(), "/generateExpSpace_log.txt")
  
  # print("CREATING EXPOSURE SPACE")
  # create the space
  spaceInfo=expSpaceSampManager(exSpArgs=exSpArgs,attInfo=attInfo,attSel=attSel,file=file)
  
  # save input arguments into the list
  spaceInfo$attPerturb <- attPerturb
  spaceInfo$attHold <- attHold
  spaceInfo$attTargetsFile <- attTargetsFile
  spaceInfo$attPerturbSamp <- attPerturbSamp
  spaceInfo$attPerturbMin <- attPerturbMin
  spaceInfo$attPerturbMax <- attPerturbMax
  spaceInfo$attPerturbType <- attPerturbType
  spaceInfo$attPerturbBy <- attPerturbBy
  
  return(spaceInfo)
  
}

createBounds <- function(attPerturb, 
                         attPerturbMin, 
                         attPerturbMax
                         ) {
  
  bounds <- list()
  for (i in 1:length(attPerturb)) {
    if (attPerturbMin[i] == attPerturbMax[i]) {
      bounds[[i]] <- attPerturbMin[i]
    } else {
      bounds[[i]] <- c(attPerturbMin[i], attPerturbMax[i])
    }
  }
  names(bounds) <- attPerturb
  return(bounds)
}

addExpArgs_attHold <- function(attPerturb = attPerturb, attHold = attHold, exSpArgs = exSpArgs) {
  
  # make attSel 
  attSel=c(attPerturb,attHold)
  attVars<-sapply(attSel,get.varType,USE.NAMES = FALSE)

  # Code for creating a default set of historical bounds
  # The bounds are set to '0' for attributes of Temperature (absolute values) & '1' for other variables
  n=length(attVars)                       
  boundsdefault=vector(length = n)
  for (i in 1:n) {
    if(attVars[i]=="Temp"){
      boundsdefault[i]=0
    } else {
      boundsdefault[i]=1
    }
  }
  tmp=as.list(boundsdefault)
  names(tmp)=attSel
    
  # this fills in samp/bounds
  # A default list of exSpArgs - specific to the particular case setup
  exSpArgsdefault=list(type="regGrid",
                       samp=rep(1,n),
                       bounds=tmp)
    
  # Combine user specified and defaults
  exSpArgs=modifyList(exSpArgsdefault,exSpArgs)
  exSpArgs$samp=c(exSpArgs$samp,rep(1,length(attHold)))
  
  return(exSpArgs)
}