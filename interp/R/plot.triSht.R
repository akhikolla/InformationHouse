
plot.triSht<-function(x,add=FALSE,xlim=range(x$x),
                      ylim=range(x$y),do.points=TRUE,
                      do.labels=FALSE, isometric=TRUE,
                      do.circumcircles=FALSE,
                      segment.lty="dashed", circle.lty="dotted", ...)
{
    if(!inherits(x,"triSht"))
        stop("x must be of class \"triSht\"")

    if(isometric){
        xlim=range(x$x)
        ylim=range(x$y)
        xrange <- diff(xlim)
        yrange <- diff(ylim)
        maxrange <- max(xrange,yrange)
        midx <- sum(xlim)/2
        midy <- sum(ylim)/2
        xlim <- midx+(xlim-midx)/xrange*maxrange
        ylim <- midy+(ylim-midy)/yrange*maxrange
    }

    if(!add)
        plot(x$x,x$y,type="n", xlim=xlim, ylim=ylim)
    if(do.points)
        points(x$x,x$y)


    segments(x$x[x$arcs[,"from"]],x$y[x$arcs[,"from"]],
             x$x[x$arcs[,"to"]],x$y[x$arcs[,"to"]], lty=segment.lty, ...)
    if(do.labels){
        midsx <- 1/3*(x$x[x$trlist[,1]]+x$x[x$trlist[,2]]+
                      x$x[x$trlist[,3]])
        midsy <- 1/3*(x$y[x$trlist[,1]]+x$y[x$trlist[,2]]+
                      x$y[x$trlist[,3]])
        text(midsx,midsy,1:x$nt,...)
        text(x$x+0.025*diff(xlim),x$y+0.025*diff(ylim),1:x$n,font=4, ...)
        arcmidsx <- 1/2*(x$x[x$arcs[,"from"]]+x$x[x$arcs[,"to"]])
        arcmidsy <- 1/2*(x$y[x$arcs[,"from"]]+x$y[x$arcs[,"to"]])
        text(arcmidsx+0.025*diff(xlim),arcmidsy+0.025*diff(ylim),
             1:x$narcs,font=3,...)
    }
    if(do.circumcircles)
        circles(x$cclist[,"x"],x$cclist[,"y"],x$cclist[,"r"],
                lty=circle.lty, ...)
}
