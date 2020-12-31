#DEMO TANK MODEL FOR APP  

# #test data
# fnam="tankhist.csv"
# dat=read.csv(fnam,header=T)
# A=Sys.time()
# out=tankPerformance(data=dat,
#                 roofArea=50,
#                 nPeople=1,
#                 tankVol=3000,
#                 write.file=TRUE,
#                 fnam="tank_demo_sam.csv")
# B=Sys.time()
# print(B-A)
# 
# systemArgs<-list(roofArea=50,   
#                  nPeople=1,
#                  tankVol=3000,
#                  firstFlush=1,
#                  write.file=T,
#                  fnam="inception.csv",
#                  metric="reliability")


#' Wrapper function for a rain water tank system model 
#'
#' \code{tankWrapper} is a wrapper function for a rainwater tank system model in foreSIGHT. This function is used in examples in function help files and vignettes.
#' This function may also be used as an example to create wrapper functions for other system models with scenarios generated using foreSIGHT in \code{R} or other programming languages. 
#' @param data data.frame; contains observed daily precipitation and temperature to be used to run the rain water tank system model in a data.frame with columns named \emph{year} \emph{month} \emph{day} \emph{P} \emph{Temp}. 
#'            Note that the first three columns of the data.frame contain the year, month, and day of observation. The columns have to be named as specified.
#'            Please refer data provided with the package that may be loaded using \code{data(tankDat)} for an example of the expected format of \code{data}.
#' @param systemArgs a list; contains the input arguments to the rain water tank system model. The valid fields in the list are:
#' \itemize{
#' \item {\code{roofArea}} {: numeric; the roof area in sq.m}
#' \item {\code{nPeople}} {: integer; number of people using water}
#' \item {\code{tankVol}} {: numeric; volume of the tank in L}
#' \item {\code{firstFlush}} {: numeric; first flush depth over roof in mm}
#' \item {\code{write.file}} {: logical; indicates whether output is to be written to file}
#' \item {\code{fnam}} {: string; name of the output file}
#' }
#' @param metrics string vector; the metrics of performance of the system model to be reported. The valid strings may be viewed using the function \code{viewTankMetrics()}
#' @return The function returns a list containing the calculated values of the performance metrics specified in \code{metrics} after running the system model.
#' @seealso \code{runSystemModel}, \code{viewTankMetrics}
#' @examples
#' # view available performance metrics
#' viewTankMetrics()
#' # load example climate data to run the system model
#' data(tankDat)
#' systemArgs <- list(roofArea = 205, nPeople = 1, tankVol = 2400, 
#' firstFlush = 2.0, write.file = FALSE)
#' tankWrapper(tank_obs, systemArgs, 
#' metrics = c("average daily deficit (L)", "reliability (fraction)"))
#' @export

tankWrapper<-function(data,
                      systemArgs,
                      metrics
) {
  
  performance<-tankPerformance(data=data,
                               roofArea=systemArgs$roofArea,   
                               nPeople=systemArgs$nPeople,
                               tankVol=systemArgs$tankVol,
                               firstFlush=systemArgs$firstFlush,
                               write.file=systemArgs$write.file,
                               fnam=systemArgs$fnam)
  
  performanceSubset <- performance[metrics]
  return(performanceSubset)
  
}


#' Prints the names of the performance metrics of the rain water tank system model 
#'
#' \code{viewTankMetrics} prints the names of the performance metrics available in the example rain water tank system model. T
#' @details This is a helper function that does not take any input arguments. The user may specify one or more of the metric names 
#' as the \code{metric} argument of \code{tankWrapper} to select the performance metrics from the tank system model.
#' to select the performance metrics.
#' @seealso \code{tankWrapper}
#' @examples
#' viewTankMetrics()
#' @export
viewTankMetrics <- function() {
  print(tankMetrics)
}


#---------------------------------------------------------
tank_model<-function(roofArea=50,   #Roof area in m2
                     nPeople=1,     #No of people(?)
                     tankVol=3000,  #tank volume in L
                     firstFlush=1,  #first flush diverter size in mm
                     rainTS=NULL,   #daily rain time series
                     tempTS=NULL,   #daily temperature time series
                     date=NULL     #date - data.frame c("year","month","day")
                     
){
  nday=length(rainTS)  #how many days to simulate
  
  #Set up vectors
  tankVolSim=rep(NA,nday)
  tankInflow=rep(NA,nday)
  tankSpill=rep(NA,nday)
  supply=rep(NA,nday)
  
  #Get initial tank volume
  initTankVol=tankVol/2.0 #start tank half full
  
  #DEMAND...
  indoorDemand=10*nPeople         #constant indoor demand (toilet flushing)
  outdoorDemandSum=739    
  seasPattern=c(0.229,0.188,0.142,0.064,0.030,0,0.001,0.007,0.014,0.049,0.107,0.171) #proportions from Goyder household water use report
  
  outdoorDemand=rep(NA,nday)
  
  for(i in 1:12){
    ind=which(date$month == i)
    outdoorDemand[ind]=seasPattern[i]*outdoorDemandSum
  }
  
  #upper temp threshold - lower temp threshold
  lowTempThresh=12; upTempThresh=28
  indUp=which(tempTS>=upTempThresh)              # made up temperature multipliers (Masters report foudn above 28 degs more watering)
  outdoorDemand[indUp]=outdoorDemand[indUp]*1.25
  indLow=which(tempTS<=lowTempThresh)
  outdoorDemand[indLow]=outdoorDemand[indLow]*0.6
  
  #combined demand
  demand=rep(indoorDemand,nday)+outdoorDemand
  
  #FIRST FLUSH  -removed at the start of each storm 
  divertAmount=firstFlush*roofArea   #1mm x roof area (Ls)
  
  
  #HOW MUCH FLOW?
  inflow=rainTS*roofArea       #mm x m2 (Ls) (100% conversion to runoff)
  roofFlow=inflow              #save as roof flow
  
  #WET DRY DAY PATTERN
  stormDur=rep(NA,nday)
  if(inflow[1]>0){stormDur[1]=1}else{stormDur[1]=0} #assume day prior to start is dry
  for(i in 2:nday){
    if(inflow[i]>0){ # wet day today
      stormDur[i]=stormDur[i-1]+1       #wet -wet or dry-wet
    }else{           # dry day today
      stormDur[i]=0  # dd or wd pattern
    }
  }
  
  #REMOVE FIRST FLUSH FROM EACH STORM
  for(i in 1:nday){
    if(stormDur[i]==1){       #first wet day in storm  (divert Amount re-set to max as emptied in between storms)
      temp=inflow[i]-divertAmount
      if(temp<0){
        divertRemainder=abs(temp)
        inflow[i]=0           #no left over flow
      }else{
        divertRemainder=0
        inflow[i]=temp        #new flow (minus divert)
      }
    }else{
      if(stormDur[i]>1){       #not first day of storm
        temp=inflow[i]-divertRemainder
        if(temp<0){
          divertRemainder=abs(temp)
          inflow[i]=0           #no left over flow
        }else{
          divertRemainder=0
          inflow[i]=temp        #new flow (minus divert)
        }
      }else{
        divertRemainder=0
      }
    }
  }
  
  #CALCULATE STARTING POINT
  if((initTankVol+inflow[1]-demand[1])<0){           # Tank has run dry
    tankVolSim[1]=0                                  # MINIMUM IS ZERO
    tankSpill[1]=0
    supply[1]=initTankVol+inflow[1]                  # Gave what you could
    tankInflow[1]=inflow[1]
  }else{
    if((initTankVol+inflow[1]-demand[1])>tankVol){ # Tank is overtopped
      tankVolSim[1]=tankVol
      tankSpill[1]=(initTankVol+inflow[1]-demand[1])-tankVol
      tankInflow[1]=tankVol-(initTankVol-demand[1])
      supply[1]=demand[1]
    }else{                                           # Tank is partially filled
      tankVolSim[1]=initTankVol+inflow[1]-demand[1]
      tankSpill[1]=0
      supply[1]=demand[1]
      tankInflow[1]=inflow[1]
    }
  }
  # loop over days
  for(i in 2:nday){
    if((tankVolSim[i-1]+inflow[i]-demand[i])<0){            # Tank has run dry
      tankVolSim[i]=0
      tankSpill[i]=0
      supply[i]=tankVolSim[i-1]+inflow[i]                   # Gave what you could
      tankInflow[i]=inflow[i]                               # all inflow travels through tank
    }else{
      if((tankVolSim[i-1]+inflow[i]-demand[i])>tankVol){    # Tank is overtopped
        tankVolSim[i]=tankVol
        tankSpill[i]=(tankVolSim[i-1]+inflow[i]-demand[i])-tankVol
        supply[i]=demand[i]
        tankInflow[i]=tankVol-(tankVolSim[i-1]-demand[i])     # took what could fit (once demand also taken)
      }else{                                                # Tank is partially filled
        tankVolSim[i]=tankVolSim[i-1]+inflow[i]-demand[i]
        tankSpill[i]=0
        supply[i]=demand[i]
        tankInflow[i]=inflow[i]                             # all inflow travels through tank
        
      }
    }
  }
  
  return(list(rainTS=rainTS,roofFlow=roofFlow,inflow=inflow,tankInflow=tankInflow,tankSpill=tankSpill,tankVolSim=tankVolSim,demand=demand,supply=supply))
}

# placed outside functions so that it can viewed using a helper function
tankMetrics <- c("volumetric reliability (fraction)",
                 "reliability (fraction)",
                 "system efficiency (%)",
                 "storage efficiency (%)",
                 "average tank storage (L)",
                 "average daily deficit (L)")
  

tankPerformance<-function(data=NULL,
                          roofArea=50,   
                          nPeople=1,
                          tankVol=3000,
                          firstFlush=1,
                          write.file=TRUE,
                          fnam="tankperformance.csv"
                          
){
  out=tank_model(roofArea=roofArea,
                 nPeople=nPeople,
                 tankVol=tankVol,
                 firstFlush=firstFlush,
                 rainTS=data$P,
                 tempTS=data$Temp,
                 date=data[,c("year","month","day")])
  
  # to test need to write to csv
  if(write.file==TRUE){
    write.csv(out, file=fnam, quote=F, row.names=F)  #write tank info to file
  }
  
  #-----tank model eval - sort of based on university of warwick metrics--------
  #reliability - fraction of days (total demand) met nSatisfied/nday
  reliability = length(which(out$demand==out$supply))/length(out$supply)
  
  #Totals (first year removed)
  Qin_tot=sum(out$tankInflow)
  Quse_tot=sum(out$supply)
  Qspill_tot=sum(out$tankSpill)
  Qroof_tot=sum(out$roofFlow)
  
  #system efficiency = (1- Quse/Qin)*100
  sysEff=Quse_tot/Qroof_tot*100
  #storage efficiency spilled v captured
  storEff=(Qspill_tot/Qin_tot)*100
  #volumetric reliability
  volRel=sum(out$supply)/sum(out$demand)
  #Av. water stored by tank on a daily basis (av. tank level)
  avTankStor=sum(out$tankVolSim)/length(out$tankVolSim)
  #avDeficit (max(demand-supply,0)/ndays)
  temp=(out$demand-out$supply);temp[which(temp<0)]=0
  avDeficit=sum(temp)/length(temp)
  
  outList <- list(volRel, reliability, sysEff, storEff, avTankStor, avDeficit)
  names(outList) <- tankMetrics
  return(outList)
  # return(list(tankMetrics[1]=volRel,tankMetrics[2]=reliability,tankMetrics[3]=sysEff,
  #             tankMetrics[4]=storEff,tankMetrics[5]=avTankStor,tankMetrics[6]=avDeficit))
}

