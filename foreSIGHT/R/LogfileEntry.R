##############################################
##   Wrapper functions for logfile entries  ##
##############################################

#CONTAINS
  #banner function for headings in log file
  #progress function to print to screen and log file
  #warning wrapper to print to screen and log file


#INPUT FILE CHECKING

#?ADD SUB DIRECTORY?



###CREATE LOG FILE
#    abs<-getwd()
#    path=paste(getwd(),"/logfile/",sep="")
#    dir.create(path,showWarnings=FALSE)
filename<-function(IOmode=NULL,arrayID=NULL) {
  if(IOmode!="dev" ){
    time<-format(Sys.time(), format="%Hh%Mm")
    day<-strsplit(as.character(Sys.Date()),split="-")[[1]][3]
    mon<-as.numeric(strsplit(as.character(Sys.Date()),split="-")[[1]][2])
    filename<-paste("logfile-",day,"-",month.abb[mon],"-",time,".txt",sep="")
  }else{
    filename<-paste0("logfile_",arrayID,".txt")
  }
  return(filename)
  #    file=paste(path,filename,sep="")
}



logfile <- function(x,file) {
  
  
  write(x,file=file,append=TRUE)
  y=""
  write(y,file=file,append=TRUE)
  
}

banner <- function(x,file) {
  
  n=nchar(x)
  write(paste(replicate(n+8,"#"),collapse=""),file=file,append=TRUE)
  write(p("##  ",x,"  ##"),file=file,append=TRUE)
  write(paste(replicate(n+8,"#"),collapse=""),file=file,append=TRUE)
  write("",file=file,append=TRUE)
}

progress <- function(x,file) {
  # cat(x,"\n")
  logfile(x,file)
  # cat("\n")
}

warn <- function(x,file) {
  cat(p("Warning: ",x),"\n")
  logfile(p("Warning: ",x),file)
  cat("\n")
}

