#Function to check supplied arguments
argument_check_calibrator<-function(names=NULL,
                                    obs=NULL,
                                    modelTag=NULL
){
  
  names<-names[names!="year"];names<-names[names!="month"];names<-names[names!="day"]
  
  #CHECKS FOR MODELTAGS
  if (modelTag[1]=="Simple-ann") { stop("Simple scaling does not require calibration - invalid request")}
  if (anyDuplicated(modelTag)!=0) {stop("There are multiple entries of the same model tag")}

    for(i in 1:length(modelTag)){
      if(sum(modelTag[i] %in% modelTaglist)==0){
        stop(paste0("modelTag ",i," unrecognised"))
      }
    }
  
  #CHECKS FOR TWO REQUESTED MODEL TYPES
  modelVars<-sapply(modelTag,get.varType,USE.NAMES=FALSE,sep="-")
  if (anyDuplicated(modelVars)!=0) {stop("There are multiple entries of a model type for one variable")}
  
  ### Checks for columns of data without model tags.
  if (length(which((names %in% modelVars)==FALSE))>0) {
    cat("There is a mismatch between provided model types and supplied variables. Calibration will only be executed for supplied model tags")
    array<-c("year","month","day",modelVars)
    obs=obs[array]
  }
    
 return(obs)
}


input_check_calibrator<- function(obs
                                   ){  #This input has its dates in first three columns
  
  #TRUNCATE TO START AND END AT WHOLE YEAR
  first=obs$year[1]
  if(obs$month[1]!=1 && obs$day[1]!=1){
    obs<-obs[obs$year!=first,]
  }
  length=nrow(obs)
  last=obs$year[length]
  if(obs$month[length]!=12 && obs$day[length]!=31){
    obs<-obs[obs$year!=last,]
  }
  length=nrow(obs)
  
  #CHECK FOR COMMON ERROR FLAGS
  obs[obs <= -99]<-NA
  
  #HOW MANY VARIABLES
  n=ncol(obs)-3
  if(n<1){stop("Not enough columns provided in input data")}
  
  #ERROR VALUES/MISSING VALUES
  M <- sapply(obs, function(x) sum(is.na(x)))
  for (Var in names(M)){
    if (M[Var]>0) {
      entry <- which(is.na(obs[,Var]))
      dates <- cbind(obs[entry,"year"],obs[entry,"month"],obs[entry,"day"])
      cat(paste0("Missing entries in ",Var),sep="\n")
      cat(apply(dates,1,paste,collapse="-"),sep="\n")
      cat("\n")
    }
  }
  
  #EXTRACT FIRST DAY AND LAST DAY
  dateS=paste(obs$year[1],str(obs$month[1],2,"0"),str(obs$day[1],2,"0"),sep="-")
  dateF=paste(obs$year[length],str(obs$month[length],2,"0"),str(obs$day[length],2,"0"),sep="-")
  temp=as.numeric(length(date_gen <- seq(as.Date(dateS),as.Date(dateF),by="day")))
  
  if(sum(M)>0){stop("There are missing data entries in the variables")} 
  if (temp<3650){cat("You have provided less than 10 years of data",sep="\n")} 
  if (length!=temp) {stop("There are missing dates from the provided data. Ensure leap days are included")}
  
  return(obs)
} #END



