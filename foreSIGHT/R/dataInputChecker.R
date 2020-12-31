#########################################
##   CHECKS FOR DATA INPUT TIMESERIES  ##
#########################################

# CONTAINS
  #input checking for data frame

#setwd("D:\\1_Code\\ScenarioNeutralGenerator\\zz_draftPackage\\data\\")
#obs<-read.csv("datJan2001-2016.csv")

input_check <- function(obs,
                        file,
                        simLengthNyrs=NULL
                        ){            #This input has its dates in first three columns

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
  
  #ERROR VALUES/MISSING VALUES
  M <- sapply(obs, function(x) sum(is.na(x)))
  for (Var in names(M)) {
    if (M[Var]>0) {
      entry <- which(is.na(obs[,Var]))
      dates <- cbind(obs[entry,"year"],obs[entry,"month"],obs[entry,"day"])
      warn(p("Missing entries in ",Var),file)
      cat(apply(dates,1,paste,collapse="-"),sep="\n")
      cat("\n")
      logfile(apply(dates,1,paste,collapse="-"),file)
      
    }
  }
  
  #EXTRACT FIRST DAY AND LAST DAY
  dateS=paste(obs$year[1],str(obs$month[1],2,"0"),str(obs$day[1],2,"0"),sep="-")
  dateF=paste(obs$year[length],str(obs$month[length],2,"0"),str(obs$day[length],2,"0"),sep="-")
  temp=as.numeric(length(date_gen <- seq(as.Date(dateS),as.Date(dateF),by="day")))
  
  
  #CHECK IF SIMLENGTHNYRS> SUPPLIED DATA
  if(!is.null(simLengthNyrs)){
    last=obs$year[1]+simLengthNyrs-1
    obsLast=obs$year[length]
    if(last<obsLast){ 
      logfile("Error: Simulation length requested is shorter than observed data",file)
      logfile("Program terminated",file)
      stop("Simulation length requested is shorter than observed data")
      }
  }
  
  
  if(sum(M)>0){
    logfile("Error: There are missing data entries in the variables",file)
    logfile("Program terminated",file)
    stop("There are missing data entries in the variables")
  
  } 
  if (temp<3650){
    warn("You have provided less than 10 years of data",file)
  } 
  if (length!=temp) {
    logfile("Error: There are missing dates from the provided data. Ensure leap days are included",file)
    logfile("Program terminated",file)
    stop("There are missing dates from the provided data. Ensure leap days are included")
  }
  
  return(list(data=obs))
} #END
