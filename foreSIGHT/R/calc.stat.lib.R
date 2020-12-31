#STATS COMPILER LIBARY
#library(moments)   #get skewness calculator

#FORTRAN FUNCS
R10calc <- function(x) {
  temp=get.nwet(data=x,threshold=10)
  return(temp)
  }
# R10calc <- function(x) {
#   n=length(x)
#   out <- .Fortran("R10calc",
#                   x=as.double(x),y=0,n=as.integer(n),package="foreSIGHT")
#   
#   return(out$y)
# }


# CDDcalc <- function(x) {
#   x[x!=0]=1
#   n=length(x)
#   out <- .Fortran("CDDcalc",
#                   x=as.double(x),y=0,n=as.integer(n),package="foreSIGHT")
#   return(out$y)
# }
# 
# 
# CDWcalc <- function(x) {
#   x[x!=0]=1
#   n=length(x)
#   out <- .Fortran("CDWcalc",
#                   x=as.double(x),y=0,n=as.integer(n),package="foreSIGHT")
#   return(out$y)
# }



# F0calc <- function(x) {
#   n=length(x)
#   out <- .Fortran("F0calc",
#                   x=as.double(x),y=0,n=as.integer(n),package="foreSIGHT")
#   return(out$y)
# }

F0calc <- function(x) {
    temp=get.below(data=x,threshold=0)
    return(temp)
}
  
# GSLcalc <- function(x) {
#   n=length(x)
#   out <- .Fortran("GSLcalc",
#                   x=as.double(x),y=0,n=as.integer(n),package="foreSIGHT")
#   
#   return(out$y)
# }
# 
# CSLcalc <- function(x) {
#   n=length(x)
#   out <- .Fortran("CSLcalc",
#                   x=as.double(x),y=0,n=as.integer(n),package="foreSIGHT")
#   
#   return(out$y)
# }

#CSL CALCULATION
CSLcalc<-function(x){
  n=length(x)
  m=n*1
  half=floor((m/2.0)-1.0)
  len=n-5.0
  sum=sum2=0
  
  for(i in 6:half){
    if((length(which(x[(i-5):i]<17)))==6){
      sum=half-i
      next
    }
  }
  
  for(i in (half+1):len){
    if((length(which(x[(i):(i+5)]>17)))==6){
      sum2=i-(half+1)
      next
    }
  }
  y=sum+sum2
  return(y)
}


# dat=c(rep(20,10),rep(0,10),rep(20,10))
# CSLcalc(x=dat,n=length(dat))

GSLcalc<-function(x){
  n=length(x)
  m=n*1
  half=floor((m/2.0)-1.0)
  len=n-5.0
  sum=sum2=0
  
  for(i in 6:half){
    if((length(which(x[(i-5):i]>5)))==6){
      sum=half-i
      next
    }
  }
  
  for(i in (half+1):len){
    if((length(which(x[(i):(i+5)]<5)))==6){
      sum2=i-(half+1)
      next
    }
  }
  y=sum+sum2
  return(y)
}

# dat=c(rep(20,10),rep(0,10),rep(20,10))
# GSLcalc(x=dat,n=length(dat))



# BASIC FORMAT & PLOTTING FUNCTIONS
d3=function(x){format(x,digits=3)} # FORMAT FUNCTION, SHORTEN TO 3 DIGITS


#Pad string in front
str<-function(x,n,pad=" "){
  temp<-as.character(x)
  nlen<-nchar(temp)
  if(nlen<n){
    for(i in 1:(n-nlen)){temp<-paste(pad,temp,sep="")}
  }
  temp
}
#Pad string at end
str.end<-function(x,n,pad=" "){
  temp<-as.character(x)
  nlen<-nchar(temp)
  if(nlen<n){
    for(i in 1:(n-nlen)){temp<-paste(temp,pad,sep="")}
  }
  temp
}

# FUNCTION TO ADD BOXPLOT WITH PROBLIMS
boxplot.func=function(z,at.pt=NULL,decile.low="10%",decile.high="90%",col=NULL,medcol=NULL,boxwex=0.7){
  boxplot.info <- boxplot(z, plot=FALSE,na.action=T);deciles <- quantile(z, probs=seq(0,1,0.05),na.rm=T)
  boxplot.info$stats[1] <- deciles[decile.low]; boxplot.info$stats[5] <- deciles[decile.high]
  bxp(boxplot.info,at=at.pt,add=T,col=col,na.action=T,range=0,boxwex=boxwex,outline=F,medcol=medcol,yaxt='n',boxfill=col)
}

# FUNCTION TO ALLOW TRANSPARENCY OF COLOUR
add.alpha <- function(COLORS, ALPHA){
  if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
  RGB <- col2rgb(COLORS, alpha=TRUE)
  RGB[4,] <- round(RGB[4,]*ALPHA)
  NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
  return(NEW.COLORS)
}

#outersect
outersect=function(x,y){
  sort(c(setdiff(x,y),setdiff(y,x)))
}

###########################################################################################################################
#CONTROLLER FUNCT - MATCHES LISTED ATT'S WITH CALCULATOR


#INPUTS - TS, INDEXES,LIST OF REQUESTED STATS
#GENERIC EXTRACTOR FUNCTION
extractor=function(func=NULL,data=NULL,indx=NULL,...){ # returns a number
  extractor.out=func(data[indx],...)
  return(extractor.out)
}

extractor.reps=function(func=NULL,data=NULL,indx=NULL,nReps=NULL,...){  # returns a vector
  temp=rep(0,nReps)
  for(rep in 1:nReps){
    dummy=data[,rep]
    temp[rep]=extractor(func=func,data=dummy,indx=indx,...)
  }
  return(temp)
}

#EXTRACTOR FOR MULTIPLE PERIODS (TEMPORARY FUNCTION here)
extractor.summaryMean<-function(func=NULL,
                                data=NULL,
                                indx=NULL,
                                nperiod=NULL,...){
  sim.series=rep(NA,nperiod)
  for(p in 1:nperiod){
    sim.series[p]=extractor(func=func,data=data,indx=indx[[p]],...)
  }
  m.series=mean(x=sim.series,na.rm=TRUE)
  return(m.series)
}


#EXTRACTOR FOR MULTIPLE PERIODS (TEMPORARY FUNCTION here)
extractor.summarySD<-function(func=NULL,
                                data=NULL,
                                indx=NULL,
                                nperiod=NULL,...){
  sim.series=rep(NA,nperiod)
  for(p in 1:nperiod){
    sim.series[p]=extractor(func=func,data=data,indx=indx[[p]],...)
  }
  m.series=sd(x=sim.series,na.rm=TRUE)
  return(m.series)
}

extractor.summaryMin<-function(func=NULL,
                              data=NULL,
                              indx=NULL,
                              nperiod=NULL,...){
  sim.series=rep(NA,nperiod)
  for(p in 1:nperiod){
    sim.series[p]=extractor(func=func,data=data,indx=indx[[p]],...)
  }
  m.series=min(x=sim.series,na.rm=TRUE)
  return(m.series)
}

extractor.summaryMax<-function(func=NULL,
                               data=NULL,
                               indx=NULL,
                               nperiod=NULL,...){
  sim.series=rep(NA,nperiod)
  for(p in 1:nperiod){
    sim.series[p]=extractor(func=func,data=data,indx=indx[[p]],...)
  }
  m.series=max(x=sim.series,na.rm=TRUE)
  return(m.series)
}

#MULTIPLE EXTRACTOR
extractor.multPeriod<-function(func=NULL,
                              data=NULL,
                              indx=NULL,
                              nperiod=NULL,...){
  
  tmp=rep(NA,nperiod)
  for(p in 1:nperiod){
    tmp[p]=extractor(func=func,data=data,indx=indx[[p]],...)
  }
  return(tmp)
}


#EXTRACTOR COEFFICIENT OF VARIATION
extractor.cv<-function(func=NULL,
                       data=NULL,
                       indx=NULL,
                       nperiod=NULL,
                       ...
  
){
  tmp=extractor.multPeriod(func=func,data=data,indx=indx,nperiod=nperiod,...)
  cv=sd(tmp,na.rm=TRUE)/mean(tmp,na.rm=TRUE)
  
  return(cv)
}

#test
#test seasons
# extractor.multPeriod(func=get.avg.tot,data=tmp$P,indx=datInd$i.ss,nperiod=4,nblocks=datInd$nyr)
# extractor.cv(func=get.avg.tot,data=tmp$P,indx=datInd$i.ss,nperiod=4,nblocks=datInd$nyr)

#test months
# extractor.multPeriod(func=get.avg.tot,data=tmp$P,indx=datInd$i.mm,nperiod=12,nblocks=datInd$nyr)
# extractor.cv(func=get.avg.tot,data=tmp$P,indx=datInd$i.mm,nperiod=12,nblocks=datInd$nyr)

#FUNCTIONS WILL PULL OUT INFORMATION FROM INPUT VECTOR
get.perc.above.thresh=function(data=NULL,
                               threshold=NULL){
  temp=length(which(data>threshold))
  if(identical(temp,integer(0))){temp=0}
  temp=temp/length(data)*100 #get percent of record above threshold   
  return(temp)
}

#FUNCTION TO DETERMINE NUMBER OF INSTANCES ABOVE A THRESHOLD - nwet
get.nwet=function(data=NULL,threshold=NULL){
  temp=length(which(data>threshold))
  if(identical(temp,integer(0))){temp=0}
  return(temp)
}

#FUNCTION TO DETERMINE NUMBER OF INSTANCES Below A THRESHOLD - nwet
get.below=function(data=NULL,threshold=NULL){
  temp=length(which(data<threshold))
  if(identical(temp,integer(0))){temp=0}
  return(temp)
}

#FUNCTION TO EXTRACT ALL AMOUNTS ABOVE A THRESHOLD
get.wet.amounts=function(data=NULL,threshold=NULL){
  temp=data[which(data>threshold)]
  return(temp)
}

#FUNCTION TO GET AVERAGE ABOVE A THRESHOLD
get.wet.average=function(data=NULL,threshold=NULL){
  ind=which(data>threshold)
  if(identical(length(ind),integer(0))){
    temp=0                          #if no wet days
  }else{
    temp=mean(data[ind],na.rm=T)
  }

  return(temp)
}

#FUNCTION TO GET STANDARD DEVIATION ABOVE A THRESHOLD
get.wet.sd=function(data=NULL,threshold=NULL){
  ind=which(data>threshold)
  if(identical(length(ind),integer(0))){
    temp=0                          #if no wet days
  }else{
    temp=sd(data[ind],na.rm=T)
  }
  return(temp)
}

#FUNCTIONS TO GET TOTALS ABOVE A THRESHOLD
get.wet.tot=function(data=NULL,threshold=NULL){
  ind=which(data>threshold)
  if(identical(length(ind),integer(0))){
    temp=0                          #if no wet days
  }else{
    temp=sum(data[ind],na.rm=T)
  }
  return(temp)
}

#FUNCTIONS TO GET TOTALS
get.tot=function(data=NULL){
  temp=sum(data,na.rm=T)
  return(temp)
}

#GET AVERAGE TOT
get.avg.tot=function(data=NULL,nblocks=NULL){
  temp=get.tot(data)/nblocks
  return(temp)
}

#FUNCTION TO GET MAXIMA ABOVE A THRESHOLD
get.wet.max=function(data=NULL,threshold=NULL){
  ind=which(data>threshold)
  if(identical(length(ind),integer(0))){
    temp=0                          #if no wet days
  }else{
    temp=max(data[ind],na.rm=T)
  }
  return(temp)
}

get.median.wet=function(data=NULL,threshold=NULL){
  ind=which(data>threshold)
  if(identical(length(ind),integer(0))){
    temp=0                          #if no wet days
  }else{
    temp=median(data[ind],na.rm=T)
  }
  return(temp)
}

get.medians=function(data=NULL){  
  temp=median(data,na.rm=T)
  return(temp)
}

get.quantile=function(data=NULL,  #vector
                      quant=NULL  #quantile (between 0.001-0.999)
                      ){  
  temp=quantile(x=data,probs=quant,na.rm=TRUE,names = FALSE)
  return(temp)
}

get.quantile.rng=function(data=NULL  #vector
                          ){  
  temp=quantile(x=data,probs=0.95,na.rm=TRUE,names = FALSE)[1]-quantile(x=data,probs=0.05,na.rm=TRUE,names = FALSE)[1]
  return(temp)
}


get.quantile.wet=function(data=NULL,  #vector
                          quant=NULL,  #quantile (between 0.001-0.999)
                          threshold=NULL  #wet day threshold
){  
  data[which(data<=threshold)]=NA
  temp=quantile(x=data,probs=quant,na.rm=TRUE,names = FALSE)
  return(temp)
}

#Culley 2020 function for amplitude
get.amplitude=function(data=NULL,
                          datInd=NULL){
  #First, get monthly totals
  nperiod=12
  monthlyTotal=rep(NA,nperiod)
  for(h in 1:nperiod){
    monthlyTotal[h]=sum(data[datInd$i.mm[[h]]],na.rm=TRUE)/datInd$nyr
  }
  #Fit harmonic to totals
  harmonicParams<-fit.harmonic.opts(nperiod=nperiod,v.stat=monthlyTotal)
  
  amplitude=harmonicParams$amp
  return(amplitude)
  
}

#skewness on wet days
get.wet.skewness=function(data=NULL,threshold=NULL){
  temp=data[which(data>threshold)]
  temp=moments::skewness(x=temp)
  return(temp)
}


#dry spell calculator - adapted from D.Guo
get.cdd<-function(data=NULL,       # w-dry status or rain vector
                  i.yy=NULL,       # year indices
                  nyr=NULL         # no. years
                  ){
  CDD <- matrix(NA,nyr,1)

  for(i in 1:nyr){
    chunk=data[i.yy[[i]]]  #chop out wd series
    Dss <- cumul_zeros(chunk) # a function to count all the lengths of continuous 0's (included in later script)
    Dss <- Dss[(which(Dss==0)-1)]
    CDD[i] <- mean(Dss[Dss!=0]) # calculate the average length of dry-spells for each year
  }
  avgCDD <- mean(CDD)
  return(avgCDD)
}

# Danlu func
# this function is to count the number of continuous 0's within a period (for calculating average length of dry spells CDD) 
cumul_zeros <- function(x)  {
  x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the 
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}

get.spell.lengths<-function(data=NULL,  # vector of rain
                            thresh=NULL,  # wetness threshold, all values below or equal to deemed dry
                            type="wet"    # get wet or dry spell length
){
  above=rep(0,length(data)) 
  ind=which(data>thresh)
  above[ind]=1 # record entries above threshold as 1
  tmp=rle(above)
  switch(type,
         "wet" ={ind.wet=which(tmp$values==1)
         spell.len=tmp$lengths[ind.wet]
         },
         "dry" ={ind.dry=which(tmp$values==0)
         spell.len=tmp$lengths[ind.dry]
         },
         -999.00)
  return(spell.len)
}

get.spell.lengths.max<-function(data=NULL,  # vector of rain
                            thresh=NULL,  # wetness threshold, all values below or equal to deemed dry
                            type="wet"    # get wet or dry spell length
){
  temp=max(get.spell.lengths(data=data,thresh=thresh,type=type),na.rm=TRUE)
  temp
}

# series=c(0.5,0.5,0.5,0.01,0.01,0.01,0.8,0.5,0.5,0.5,0,0,0,2,2,0,2)
# get.spell.lengths(data=series,thresh=0.01,type="dry")
# mean(get.spell.lengths(data=series,thresh=0.01,type="wet"),na.rm=TRUE)

p<-function(...){paste(...,sep="")}   # PASTE FUNCTION

get.tag.varType<-function(attrib=NULL, # attribute name
                          sep="-"){
  varType=strsplit(x = attrib,split=sep)[[1]][1]
  return(varType)
}

#categorise func
categ.fun=function(perf.lim=c(5,10), # performance limits (<=5% good, <=10& fair, >10% poor)
                   rel.diff=NULL     #relative difference to classify
  
){
  perf="poor"                        #start off at "poor"
  if(abs(rel.diff)<=perf.lim[1]){
    perf="good"
  }else{
    if(abs(rel.diff)<=perf.lim[2]){
      perf="fair"
    }else{
      perf="poor"
    }
  }
  
  return(perf)
}


#plot classifer chart element - one of many in grid of classifiers

plot.attrib.perf.solo=function(rel.diff,                   #relative difference - scalar
                               perf.lim=c(5,10),           #performance limits - good, fair, poor beyond
                               targetType=NULL,
                               att.name=NULL,              #string that will label plot
                               prim.lab=NULL,              #primary label
                               cex.mult=3,
                               cex.mult.sub=1.1,
                               y.text=0.05,
                               mtext.line=0.35
                               ){
  #COLOR RAMP
  # traffic.col=c("chartreuse3","gold1","red1")
  
  #MAKE VECTOR X
  att.cat=categ.fun(perf.lim,rel.diff)
  if(att.cat=="poor"){ind=3};if(att.cat=="fair"){ind=2};if(att.cat=="good"){ind=1}
  x=rep(0,3)  #make blank x vector
  x[ind]=100  #update to reflect att.cat
  
  
  #PLOT BARPLOT
  barplot(height = cbind(x = x/100),horiz=T,xaxt='n',yaxt='n',
          beside = FALSE,width = c(0.1),col = traffic.col,
          args.legend = list(x = "topleft"))
  
  #ADD TEXT ANNOTATIONS
  if(ind==3){bg.col="black"; front.col="white"}else{bg.col="white";front.col="black"}   #text background colour updater
  if((targetType == "pc")|(targetType == "frac")){
    text(labels=paste(format(rel.diff,digits=2),"%",sep=""),x=0.5,y=y.text,cex=cex.mult,pos=3,col=front.col,bg=bg.col)
  }else{
    text(labels=paste(format(rel.diff,digits=2)," delta",sep=""),x=0.5,y=y.text,cex=cex.mult,pos=3,col=front.col,bg=bg.col)
  }
  
  #ADD SUBTITLE
  mtext(text=att.name,side=1,at=0.5,cex=cex.mult.sub,line=mtext.line)
  
  #ADD TITLE
  if(!is.null(prim.lab)){
    mtext(text=prim.lab,side=3,at=0.5,cex=cex.mult.sub,line=mtext.line)
  }
}

#calculation of percentage change
pc.calc<-function(sim=NULL,     #simulate point 
                  target=NULL   #target point
                  ){
  
  pc.diff=(sim-target)/target*100 #calc percen diff from target
  
  }
#pc.calc(sim,target)

#calculation of percentage change
abs.diff.calc<-function(sim=NULL,     #simulate point 
                       target=NULL    #target point
){
  
  abs.diff=(sim-target)      #calc abs diff from target
  
}
 
########################################
# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}