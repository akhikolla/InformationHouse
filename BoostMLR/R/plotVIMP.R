plotVIMP <- function(vimp_Object,
                     xvar.names = NULL,
                     cex.xlab = NULL,
                     ymaxlim = 0,
                     yminlim = 0,
                     main = "Variable Importance (%)",
                     col = grey(.80),
                     cex.lab = 1.5,
                     ylbl = NULL,
                     legend_placement = NULL,
                     plot.it = TRUE,
                     path_saveplot = NULL,
                     Verbose = TRUE)
{
  
  ymaxtimelim <- 0
  subhead.cexval <- 1
  yaxishead <- NULL
  xaxishead <- NULL
  subhead.labels <- c("Time-Interactions Effects","Main Effects")
  Adj <- -0.5
  cex.main <- 1
  seplim <- NULL
  
  if(missing(vimp_Object) ){
    stop("vimp is not provided")
  }
  
  if(is.list(vimp_Object)){
    L <- length(vimp_Object)
  } else 
    {
    stop("vimp_Object must be a list of length equal to number of response variables")
  }
  
  vimp <- matrix(unlist(lapply(1:L,function(l){
    vimp_Object[[l]][,1]
  })),ncol = L,byrow = FALSE)
  colnames(vimp) <- names(vimp_Object)
  
  if(is.null(xvar.names)){
    xvar.names <- rownames(vimp_Object[[1]])
  }
  MainEffect <- TRUE
  vimp <- vimp*100
  if(plot.it){
    
    if(is.null(path_saveplot)){
      path_saveplot <- tempdir()
    }
    pdf(file = paste(path_saveplot,"/","VIMP.pdf",sep=""),width = 14,height = 14)
    
    if(MainEffect){
      for(l in 1:L){
        if(is.null(ylbl)){
          ylbl <- ""
        }
        ylim <- range(vimp[,l]) + c(yminlim,ymaxlim)
        yaxs <- pretty(ylim)
        yat <- abs(yaxs)
        bp <- barplot(as.matrix(vimp[,l]),beside=T,col=col,ylim=ylim,yaxt="n",ylab = ylbl,main = main,cex.main = cex.main,cex.lab=cex.lab)
        text(c(bp), if(is.null(legend_placement)) pmax(as.matrix(vimp[,l]),0) else legend_placement, rep(xvar.names, 3),srt=90,adj=Adj,cex= if(!is.null(cex.xlab)) cex.xlab else 1 )
        axis(2,yaxs,yat)
      }
    } else
    {
      p <- 1
      n.vimp <- 1
      for(l in 1:L){
        vimp.x <- vimp[1:p,l]
        vimp.time <- vimp[-c(1:p),l]
        ylim <- max(c(vimp.x,vimp.time)) * c(-1, 1) + c(-ymaxtimelim,ymaxlim)
        yaxs <- pretty(ylim)
        yat <- abs(yaxs)
        if(is.null(yaxishead)){
          yaxishead <- c(-ylim[1],ylim[2])
        }
        if(is.null(xaxishead)){
          xaxishead <- c(floor(n.vimp/4),floor(n.vimp/4))
        }
        bp1 <- barplot(pmax(as.matrix(vimp.x),0),beside=T,col=col,ylim=ylim,yaxt="n",ylab = "",cex.lab=cex.lab,
                       main = main,cex.main = cex.main)
        text(c(bp1), pmax(as.matrix(vimp.x),0), rep(xvar.names, 3),srt=90,adj=Adj,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
        text(xaxishead[2],yaxishead[2],labels = subhead.labels[2],cex = subhead.cexval)
        bp2 <- barplot(-pmax(as.matrix(vimp.time),0),beside=T,col=col,add=TRUE,yaxt="n")
        text(c(bp2), -pmax(as.matrix(vimp.time),0), rep(xvar.names, 3),srt=90,adj=Adj,yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
        text(xaxishead[1],-yaxishead[1],labels = subhead.labels[1],cex = subhead.cexval)
        axis(2,yaxs,yat)
      }
    }
    dev.off()
    if(Verbose){
      cat("Plot will be saved at:",path_saveplot,sep = "")
    }
  }
}
