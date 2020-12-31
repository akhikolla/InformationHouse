#' @rdname kcpRS
#' @param x An object of the type produced by \code{kcpRS}
#' @param ... Further plotting arguments.
#' @importFrom graphics plot abline axis legend lines mtext par
#' @importFrom RColorBrewer brewer.pal
#' @export

plot.kcpRS<-function(x,...){

    RS<-as.data.frame(x$RS)
    cps<-x$changePoints
    wsize<-x$wsize
    main_lab<-x$RS_name

    N<-nrow(RS)+wsize-1
    v<-ncol(RS)
    start<-wsize/2
    end<-N-(wsize/2)
    time<-seq(start,end,1)

    #xlim_max


    #ylim max
      max<-max(RS[,1])            #max of the first variable
      nPanels=v+1                 #no. of panels: no. of variables (v) + 1(for the lower values of the variable at the bottom)
      add<-ifelse(max>1,max-1,0)  #additional height for the topmost panel
      ylim_max<-nPanels+add

    #ylim min
      min<-min(RS[,v])              #min of the last variable
      add<-ifelse(min<(-1),min+1,0) #what to add to extend the last panel
      ylim_min<-add                 #additional height (downward direction) for the lowermost panel


    #plot
      par(mar=c(3,5,2,2))
      plot(time,RS[,1],type="n",main=paste0("Running ", main_lab, "s"),ylab=""
           ,xlab="",yaxt="n",ylim=c(ylim_min,ylim_max)
           ,xlim=c(0,N+N/6), cex.axis=.8)

          abline(h=seq(1,v,1),col="gray") #horizontal (zero) lines for each var
          y_ticks=seq(v,1,-1)
          y_labs=colnames(RS)

          axis(2, at=y_ticks, labels=y_labs,las=1,cex.axis=.8,hadj=1,font=1)
          mtext("Time", side=1, line=2,cex = .7)

          #liner colors: to differentiate lines since sometimes they can overlap (if colors are the same, difficult to monitor)
              kulay=brewer.pal(n=9, name = 'Set1')    #set color scheme
              kulay=kulay[-1]                         #erase the first color (RED), since change points will be colored RED

              if (v<8){kulay_final=kulay[1:v]}        #color set only contains 8 unique colors
              if (v>8){nrep=ceiling(v/8)              #replicate colors if there are more than 8 variables
              kulay_final=rep(kulay,nrep)}


      for (i in 1:v){lines(time,RS[,i]+((v)-i+1),col=kulay_final[i],lwd=2)}   #lines
      if (length(cps)>0){abline(v=cps,lty=2,col=2)}                           #change points

      legend(N-(wsize/2)+(N/30),ylim_max,y_labs,lwd=5,col=kulay_final,cex=.6)
}

