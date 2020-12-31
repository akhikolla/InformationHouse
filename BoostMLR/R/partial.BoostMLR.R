#-------------------------------------------------------------------------------------------
# Partial plot function

partial.BoostMLR <- function(Object,
                             xvar.name,
                             n.x = 10,
                             n.tm = 10,
                             x.unq = NULL,
                             tm.unq = NULL,
                             Mopt,
                             plot.it = TRUE,
                             path_saveplot = NULL,
                             Verbose = TRUE,
                             ...)
{
  
  if (missing(Object)) {
    stop("Object is missing")
  }
  
  x  <- Object$x
  tm <- Object$tm
  id <- Object$id
  
  if (missing(xvar.name)) {
    stop("xvar.name is missing" )
  }
  
  if(length(xvar.name) > 1){
    stop("Use single covariate in xvar.name")
  }
  
  x_Names <- Object$x_Names
  xvar.name <- intersect(xvar.name, x_Names)
  if (length(xvar.name) == 0) {
    stop("xvar.name do not match original variable names")
  }
  
  user.option <- list(...)
  prob_min <- is.hidden.prob_min(user.option)
  prob_max <- is.hidden.prob_max(user.option)
  
  if( is.null(x.unq) ){
    xvar <- x[,xvar.name,drop = TRUE]
    n.x.unq <- length(unique(xvar))
    if(n.x.unq <= n.x){
      x.unq <- sort(unique(xvar))
    } else 
      {
      x.unq <- sort(unique(xvar))[unique(as.integer(quantile(1: n.x.unq ,probs = seq(prob_min,prob_max,length.out = min(n.x,n.x.unq)   ))))]  
    }
  }
  n.x <- length(x.unq)

  if( is.null(tm.unq) ){
    n.tm.unq <- length(unique(tm))
    if(n.tm.unq <= n.tm){
      tm.unq <- sort(unique(tm))
    } else 
      {
        tm.unq <- sort(unique(tm))[unique(as.integer(quantile(1: length(unique(tm)) ,probs = seq(0,0.9,length.out = min(n.tm,n.tm.unq) ))))]
    }
  }
  n.tm <- length(tm.unq)
  
  if(missing(Mopt)){
    Mopt <- Object$M
  }
  
  L <- Object$Grow_Object$Dimensions$L
  
  p.obj <- lapply(1:n.x,function(i){
    new_x <- x
    new_x[,xvar.name] <- x.unq[i]
    lapply(1:n.tm,function(j){
      tm <- rep(tm.unq[j],nrow(x))
      colMeans(predictBoostMLR(Object = Object,x = new_x,tm = tm,id = id,M = Mopt,importance = FALSE)$mu,na.rm = TRUE)
    })    
  })
  
  pList <- lapply(1:L,function(l){
    pMat <- matrix(NA,nrow = n.x,ncol = n.tm)
    for(i in 1:n.x){
      for(j in 1:n.tm){
        pMat[i,j] <- p.obj[[i]][[j]][l]
      }
    }
    pMat
  })
  
  sList <- lapply(1:L,function(l){
    sMat <- matrix(NA,nrow = n.x,ncol = n.tm)
    for(i in 1:n.x){
      y.lo <- lowess(x = tm.unq,y = pList[[l]][i,,drop = TRUE] )$y
      sMat[i,] <- y.lo
    }
    sMat
  })
  
  if(plot.it){

  if(is.null(path_saveplot)){
      path_saveplot <- tempdir()
  }
  pdf(file = paste(path_saveplot,"/","PartialPlot.pdf",sep=""),width = 14,height = 14)
    for(l in 1:L){
  filled.contour(x = tm.unq,
                 y = x.unq,
                 z = t(sList[[l]]),
                 xlim = range(tm.unq, finite = TRUE) + c(-0, 0),
                 ylim = range(x.unq, finite = TRUE) + c(-0, 0),
                 color.palette =
                   colorRampPalette(c("yellow", "red")),
                 xlab = "Time",
                 ylab = "x",
                 main = "PartialPlot",
                 cex.main = 2,
                 cex.lab = 1.5,
                 plot.axes = {
                 axis(1,cex.axis = 1.5)
                 axis(2,cex.axis = 1.5)})
    }
    dev.off()
      if(Verbose){
     cat("Plot will be saved at:",path_saveplot,sep = "")
  }
  }
  obj <- list(x.unq = x.unq,
              tm.unq = tm.unq,
              pList = pList,
              sList = sList)
  
  invisible(obj)
}
