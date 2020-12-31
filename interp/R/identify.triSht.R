identify.triSht<-function(x,...)
  {
    if(!inherits(x,"triSht"))
      stop("x must be of class \"tri\"")
    labels<-paste("(",round(x$x,5),",",round(x$y,5),")", sep ="")
    identify(x$x,x$y,labels=labels)
  }
