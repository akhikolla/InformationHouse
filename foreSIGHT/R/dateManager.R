#######################################
##   DATE MANAGER FUNCTION LIBRARY   ##
#######################################

#CONTAINS
  # dateExtender()
  # extendDates()
  # makeDates() - produces dates data.frame (year, month, day columns)
  # mod.get.date.ind() - grab a variety of date indices from multiple modelTags
  # mod.get.date.ind.extnd() - grab a variety of date indices from multiple modelTags with date extension
  # get.date.ind() - master function to grab a variety of date indices
  # get.month.ind()
  # get.year.ind()
  # get.seas.ind() - (note seasons currently correspond to southern hemisphere)
  # julian.day()
  # split.ts()  - divide year into even(ish) chunks for harmonic fit
  # get.period.ind() - groups indices with the same period assignment

##############################################################################
dateExtender<-function(obs=NULL,
                       simLengthNyrs=NULL,
                       file=NULL,
                       modelTag=NULL
                       # Anjana: Commented - removed multiple ways of extending dates for simplicity
                       #datStart=NULL,
                       #datFinish=NULL
                       ){
  #EXTEND DATES IF NEEDED
  if(!is.null(simLengthNyrs)){
    if(modelTag[[1]] != "Simple-ann"){
      dateExtnd=extendDates(simLengthNyrs=simLengthNyrs,dd=obs$day,mm=obs$month,yy=obs$year)
      progress("Extending dates",file)
    }else{
      dateExtnd=obs[,c("year","month","day")]                                              # make the same as observed
      progress("Length of time series cannot be increased using simple scaling",file)
    }
  # }else if(!is.null(datStart)){
  #   if(modelTag[[1]] != "Simple-ann"){
  #     dateExtnd=makeDates(datStart=datStart,datFinish=datFinish)
  #     progress("Extending dates",file)
  #   }else{
  #     dateExtnd=obs[,c("year","month","day")]                                              # make the same as observed
  #     progress("Length of time series cannot be increased using simple scaling",file)
  #   }
  } else {
    dateExtnd=obs[,c("year","month","day")]                                               # make the same as observed
  }
  return(dateExtnd)
}

extendDates<-function(simLengthNyrs=NULL,
                      dd=NULL,
                      mm=NULL,
                      yy=NULL
){
  ndays=length(mm)
  dateS=paste(yy[1],str(mm[1],2,"0"),str(dd[1],2,"0"),sep="-")
  dateF=paste((yy[1]+simLengthNyrs-1),str(mm[ndays],2,"0"),str(dd[ndays],2,"0"),sep="-")
  date_gen <- seq(as.Date(dateS),as.Date(dateF),by="day")
  
  day <- as.numeric(format(date_gen,"%d"))
  month<- as.numeric(format(date_gen,"%m"))
  year<- as.numeric(format(date_gen,"%Y"))
  dates=data.frame(year,month,day)
  return(dates)
}
#TEST
# tester=extendDates(simLengthNyrs=100,dd=obs$day,mm=obs$month,yy=obs$year)

makeDates<-function(datStart=NULL,
                    datFinish=NULL){
  date_gen=seq(as.Date(datStart),as.Date(datFinish),by="day")
  day <- as.numeric(format(date_gen,"%d"))
  month<- as.numeric(format(date_gen,"%m"))
  year<- as.numeric(format(date_gen,"%Y"))
  dates=data.frame(year,month,day)
  return(dates)
}

#get date info across multiple models
mod.get.date.ind<-function(obs=NULL,
                           modelTag=NULL,
                           modelInfo=NULL,
                           southHemi=TRUE){
  nMod=length(modelTag)   #how many models
  datInd=list()
  datInd[["obs"]]=get.date.ind(dd=obs$day,mm=obs$month,yy=obs$year,nperiod=12,southHemi=southHemi)              #make obs based datInd
  for(i in 1:nMod){
    datInd[[modelTag[i]]]=get.date.ind(dd=obs$day,mm=obs$month,yy=obs$year,nperiod=modelInfo[[modelTag[i]]]$nperiod,southHemi=southHemi)    # FROM dateManager.R
    # datInd[[modelTag[i]]]$i.mod=datInd[[modelTag[i]]]$i.pp  #add on i.mod
  }
  return(datInd)
}

#get date info across multiple models - dates extended
mod.get.date.ind.extnd<-function(obs=NULL,
                                 dateExtnd=NULL,
                                 modelTag=NULL,
                                 modelInfo=NULL,
                                 southHemi=TRUE,
                                 simLengthNyrs=NULL,
                                 file=NULL
                                 ){
  datInd=list()
  datInd[["obs"]]=get.date.ind(dd=obs$day,mm=obs$month,yy=obs$year,nperiod=12,southHemi=southHemi)              #make obs based datInd
  for(i in 1:length(modelTag)){
    datInd[[modelTag[i]]]=get.date.ind(dd=dateExtnd$day,mm=dateExtnd$month,yy=dateExtnd$year,nperiod=modelInfo[[modelTag[i]]]$nperiod,southHemi=TRUE)          # FROM dateManager.R
  }
  return(datInd)
}

#Get dat indices
get.date.ind<-function(dd=NULL,
                       mm=NULL,
                       yy=NULL,
                       nperiod=NULL,     #information surrounding fitted model
                       southHemi=TRUE
                       ){
  
  ndays=length(dd)                  # get number of days on record
  nyr=yy[ndays]-yy[1]+1             # get number of years on record
  i.mm=get.month.ind(mm=mm)         #get indices for months
  i.yy=get.year.ind(yy=yy,nyr=nyr,n=ndays)  #get indices for years
  if(southHemi==TRUE){
    i.ss=get.seas.ind(i.mm=i.mm)    #get indices for seasons
  }else{
    print("warning check seasons")  #warning not southern hemisphere
  }
  
  dateS=paste(yy[1],mm[1],dd[1],sep="-")
  dateF=paste(yy[ndays],mm[ndays],dd[ndays],sep="-")
  jj=julian.day(dateS=dateS,dateF=dateF)
  
  i.pp=list()
  if((nperiod==1)|(nperiod==4)|(nperiod==12)){
    if(nperiod==1){i.pp[[1]]=seq(1,ndays)}  # annual model case
    if(nperiod==4){i.pp=i.ss}             # seasonal model case
    if(nperiod==12){i.pp=i.mm}            # monthly model case (not currently in use)
  }else{
    harInd= split.ts(nperiod=nperiod,jj=jj)  # alternative period split - calculate indices
    i.pp=get.period.ind(har.period=harInd,nperiod=nperiod)
  }

  datInd=list(ndays=ndays,
              nyr=nyr,
              i.mm=i.mm,
              i.yy=i.yy,
              i.ss=i.ss,
              i.pp=i.pp)
  return(datInd)
}

get.month.ind<-function(mm=NULL  # ts vector of months
){
  i.mm=NULL
  for(m in 1:12) i.mm[[m]]=which(mm==m)         # CREATE MONTHLY INDICES
  return(i.mm)
}

get.year.ind<-function(yy=NULL,   # ts vector of years
                       nyr=NULL,  # nyears oon record
                       n=NULL #no. of days on record
                       ){
  years=seq(yy[1],yy[n])                                   # GET VECTOR OF YEARS
  i.yy=NULL
  for(Y in 1:nyr) i.yy[[Y]]=which(yy==years[Y])            # CREATE MONTHLY INDICES
  return(i.yy)
}

get.seas.ind<-function(i.mm=NULL # list of days sorted by month
){
  
  #NOTE SOUTHERN HEMISPHERE HERE
  # define months belonging to each season
  seas <- t(matrix(c(9,10,11,#SPR -SON
                     12,1,2,#SUM -DJF
                     3,4,5,#AUT -MAM
                     6,7,8),#WIN -JJA
                   nrow=3,ncol=4))
  
  i.ss=NULL
  for(s in 1:4) i.ss[[s]]=c(i.mm[[seas[s,1]]],i.mm[[seas[s,2]]],i.mm[[seas[s,3]]]) #CREATE SEASONAL INDICES i.ss[[1]]-i.ss[[4]] (not contiguous, needs a sort)
  for(s in 1:4){
    tmp=i.ss[[s]]
    tmp=sort(tmp)
    i.ss[[s]]=tmp  #put in daily order
  }
  rm(tmp)
  return(i.ss)
  
}

julian.day<-function(dateS=NULL,  #start date e.g. "1995-01-01"
                     dateF=NULL
){
  date_gen <- seq(as.Date(dateS),as.Date(dateF),by="day")
  ndays <- length(date_gen)
  jj <- as.numeric(format(date_gen,"%j"))
  return(jj)
}
# jj=julian.day(dateS="1995-01-01",dateF="2004-12-31")

split.ts<-function(nperiod=26,   #no. of periods to divide year over
                   jj=NULL       #vector of julian day values
){
  nd.per=floor(366/nperiod)  #no. days in period
  nd.year=nperiod*nd.per     #est no. days in year
  short=366-nd.year               #no. days short from 366
  
  indl <- NULL
  for (i in 1:nperiod) {indl <- c(indl,rep(i,nd.per))}
  indl <- c(indl,rep(nperiod,short)) #add missing days on end
  
  harInd <- rep(NA,length(jj))
  for (i in 1:366) {
    tmpInd=which(jj==i)
    harInd[tmpInd] <- indl[i]       #determine which day belongs to which period
  }
  return(harInd)
}
#harInd= split.ts(nperiod=26,jj=jj)

get.period.ind<-function(har.period=NULL,  # ts vector of period assigned
                         nperiod=NULL      # number of periods used
){
  i.hh=NULL
  for(h in 1:nperiod) i.hh[[h]]=which(har.period==h)         # CREATE MONTHLY INDICES
  return(i.hh)
}
#i.hh=get.period.ind(har.period=harInd,nperiod=26)

