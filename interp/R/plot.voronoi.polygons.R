"plot.voronoi.polygons" <- function(x,which,
                                    color=TRUE,
                                    isometric=TRUE,...){
  lx <- length(x)
  if(missing(which))
      which <- 1:lx
  ## exclude border polygons represented as NULL,
  ## intersect ensures working behaviour for eventuelly
  ## given argument "which" (otherwise 1:lx)
  which <-  intersect(which,(1:lx)[!(unlist(lapply(x, is.null)))])
  if(any(is.na(which)))
      stop("border polygons may not be choosen to plot")
  lw <- length(which)
  lmax <- function(x)
    apply(x,2,max)
  lmin <- function(x)
    apply(x,2,min)
  lmean <- function(x)
    apply(x,2,mean)
  xy.max <- apply(sapply(x[which],lmax),1,max)
  xy.min <- apply(sapply(x[which],lmin),1,min)
  xy.mean <- sapply(x[which],lmean)

  xlim=c(xy.min["x"]-
         0.1*(xy.max["x"]-xy.min["x"]),
         xy.max["x"]+
         0.1*(xy.max["x"]-xy.min["x"]))
  ylim=c(xy.min["y"]-
         0.1*(xy.max["y"]-xy.min["y"]),
         xy.max["y"]+
         0.1*(xy.max["y"]-xy.min["y"]))
  if(isometric){
      xrange <- diff(xlim)
      yrange <- diff(ylim)
      maxrange <- max(xrange,yrange)
      midx <- sum(xlim)/2
      midy <- sum(ylim)/2
      xlim <- midx+(xlim-midx)/xrange*maxrange
      ylim <- midy+(ylim-midy)/yrange*maxrange
  }
  plot(x[[which[1]]],type="n",xlim=xlim,
       ylim=ylim,...)
  colors <- heat.colors(lw)
  j <- 0
  for(i in which){
      j <- j+1
      polygon(x[[i]],col=colors[j])
      text(xy.mean[,j]["x"],xy.mean[,j]["y"],i)
  }
  title(paste("plot of",deparse(substitute(x))))
}
