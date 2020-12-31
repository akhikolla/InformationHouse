#attribute calculator#

#' Calculates the attributes of the hydroclimate time series
#'
#' \code{calculateAttributes} calculates the specified attributes of the input daily hydroclimate time series.
#' @param climateData data.frame; daily climate data, the attributes of which are to be calculated. \code{climateData} should be a data.frame with columns named \emph{year} \emph{month} \emph{day} \emph{*variable_name1*} \emph{*variable_name2*}. 
#'            Use \code{viewModels()} to view the valid variable names. Note that the first three columns of the data.frame contain the year, month, and day of the data. The columns have to be named as specified.
#'            Please refer data provided with the package that may be loaded using \code{data(tankDat)} for an example of the expected format of \code{climateData}.
#' @param attSel a vector; specifying the names of the attributes to be calculated.
#' @param startYr a number (default \code{NULL}); to specify the starting year to subset \code{climateData} if required. 
#' If \code{NULL}, \code{startYr} is starting year in the input \code{climateData}.
#' @param endYr a number (default \code{NULL}); to specify the ending year to subset \code{climateData} if required. 
#' If \code{NULL}, \code{endYr} is last year in the input \code{climateData}.
#' @return The function returns a vector of length equal to the number of selected attributes (\code{attSel}), named using the names of the attributes.
#' @seealso \code{viewAttributes}, \code{viewAttributeDef}
#' @examples
#' # view the names of available attributes
#' viewAttributes()
#' # load example climate data available in the package
#' data("tankDat")
#' attSel <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
#' tank_obs_atts <- calculateAttributes(tank_obs, attSel = attSel)
#' @export

calculateAttributes<-function(climateData,                    # input data in the format of tank_obs (can be reference, obs, or, future projections)
                              attSel,                         # vector of selected attributes
                              startYr = NULL,                 # changed slice & window to startYr and endYr
                              endYr = NULL                    #       - can specify one without the other as well 
                              
){
  IOmode="verbose"
  arrayID=NULL
  simLengthNyrs=NULL
  file<-filename(IOmode=IOmode,arrayID=arrayID)
  
  input<-input_check(climateData,file,simLengthNyrs) #Checks for missing values/full years of data
  obs<-input$data
  
  if(!is.null(startYr)) {
    if (!is.null(endYr)) {
      if(startYr > endYr) stop("startYr should be less than or equal to endYr.")
    }
    if (startYr < obs$year[1]) warning("startYr is less than the starting year of climateData.")
  } else {
    startYr <- obs$year[1]
  }
  
  if(!is.null(endYr)) {
    if (endYr > tail(obs$year,1)) warning("endYr is greater than the last year of climateData.")
  } else {
    endYr <- tail(obs$year,1)
  }
  obs<-obs[which(obs$year>=startYr&obs$year<=endYr),]
  
  # if(!is.null(slice)){                       #For non historical records     ***** old code using 'slice' and 'window'
  #   start=slice-window
  #   text<-paste("Note: Window is set ",window," years before slice",sep="")
  #   #    progress(text,file)
  #   obs<-obs[which(obs$year>=start&obs$year<=slice),]
  # }
  
  #Get necessary variables for historical atts
  attInfo=attribute.info.check(attSel=attSel)              # vector of selected attributes (strings)
  simVar<-attInfo$varType
  simVar<-unique(simVar)
  nvar<-length(simVar)
  
  #  banner("INDEXING DATES",file)
  #  progress("Indexing dates...",file)                       # USE NEW APPENDED/CHECKED DATA
  yy=obs$year;mm=obs$month;dd=obs$day  
  # STORE DATE VECTORS
  
  datInd=list()
  datInd[["obs"]]=get.date.ind(dd=dd,mm=mm,yy=yy,nperiod=12,southHemi=TRUE)
  attInd=get.att.ind(attInfo=attInfo,simVar=simVar)
  
  #  progress("Dates indexed OK",file)
  
  #  banner("OBSERVED BASELINE ATTRIBUTE CALCULATION",file)
  #  progress("Calculating attributes...",file)
  
  attObs=list()                            #make this into its own function (also inserted into model sequencer)
  for(i in 1:nvar){
    attObs[[i]]=attribute.calculator(attSel=attSel[attInd[[i]]],
                                     data=obs[[simVar[i]]],
                                     datInd=datInd[["obs"]],
                                     attribute.funcs=attribute.funcs) 
  }
  
  #  progress("Attributes calculated OK",file)   #NEED SOME ACTUAL CHECKING HERE BEFORE PRONOUNCING OK
  attObs=unlist(attObs)
  names(attObs)=attSel
  #  progress(attObs)
  
  return(attObs)
}