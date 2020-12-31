### plotting segmentation results and segments



### PLOT UTILITIES
#' Switch between plot devices.
#' @param file.name file name without suffix (.png, etc)
#' @param type plot type: png, jpeg, eps, pdf, tiff or svg
#' @param width figure width in inches
#' @param height figure height in inches
#' @param res resolution in ppi (pixels per inch), only for types png,
#' jpeg and tiff
#' @export
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100) {
  file.name <- paste(file.name, type, sep=".")
  if ( type == "png" )
    grDevices::png(file.name, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    grDevices::postscript(file.name, width=width, height=height,paper="special",
                          horizontal = FALSE, onefile = FALSE)
  if ( type == "pdf" )
    grDevices::pdf(file.name, width=width, height=height)
  if ( type == "tiff" )
    grDevices::tiff(file.name, width=width, height=height, units="in", res=res)
  if ( type == "svg" )
    grDevices::svg(file.name, width=width, height=height)
  if ( type == "jpeg" )
    grDevices::jpeg(file.name, width=width, height=height, units="in", res=res)
}

## from lib TeachingDemos; used in plotFeatures
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}


### PLOTTING RESULTS - from genomeBrowser 20161102

## plot features as arrows; kept private since official version
## of this function is maintained in package genomeBrowser
## 20170207 : using ... to pass lwd on to arrow plots: useful!
segment.plotFeatures <- function(data, coors, types, strand,
                         typord=FALSE, cuttypes=FALSE,
                         names=FALSE, legend=FALSE, axis1=FALSE, ylab=NA,
                         args,line=1,ycx=1,tcx=1,tpos=NULL,
                         columns=c(name="name",chr="chr",strand="strand",
                                   start="start",end="end",type="type",
                                   color="color"),
                         ...) {

  ## parse arguments, which will also override other 
  if ( !missing(args) ) {    
    for ( i in seq_along(args) )       
      assign(names(args)[i],args[[i]])
  }

  ## get all available types
  if ( missing(types) ) 
    types <- unique(data[,columns["type"]])

  ## use only features from the requested strand
  if ( !missing(strand) )
    data <- data[data[,columns["strand"]]==strand,]

  ## coordinates
  chr <- as.character(coors[1])
  xlim <- sort(coors[2:3])

  ## cut to xlim
  coor <- cbind(as.numeric(data[,columns["start"]]),
                as.numeric(data[,columns["end"]]))
    
  ## find and filter features overlapping with requested coors
  lovl <- apply(coor,1,max) >= xlim[1] & apply(coor,1,max) <= xlim[2]
  rovl <- apply(coor,1,min) >= xlim[1] & apply(coor,1,min) <= xlim[2]
  insd <- lovl & rovl
  span <- apply(coor,1,min) <= xlim[1] & apply(coor,1,max) >= xlim[2]
  ovl <- ( lovl | rovl ) | span
  ## index features to determine where to plot names
  if ( names )
    data <- cbind(data, lovl=lovl, rovl=rovl, insd=insd, span=span)

    ## filter overlapping features
    ## to avoid factor confusion!
    if ( "chr"%in%names(columns) ) {
        chrs <- as.character(data[,columns["chr"]])
        ovl <- chrs == chr & ovl
    }
    feat <- data[ovl,, drop=FALSE]

    ## get new coordinates
    coor <- cbind(as.numeric(feat[,columns["start"]]),
                  as.numeric(feat[,columns["end"]]))

    ## order start/end by strand (arrow direction!)
    if ( "strand"%in%names(columns) ) {
        str <- as.character(feat[,columns["strand"]])
        ## TODO draw different arrow for strand=.
        coor[str==".",] <-
            t(apply(coor[str==".",,drop=FALSE],1,
                    function(x) sort(x,decreasing=FALSE)))
        coor[str=="+"|str=="+1",] <-
            t(apply(coor[str=="+"|str=="+1",,drop=FALSE],1,
                    function(x) sort(x,decreasing=FALSE)))
        coor[str=="-"|str=="-1",] <-
            t(apply(coor[str=="-"|str=="-1",,drop=FALSE],1,
                    function(x) sort(x,decreasing=TRUE)))
    }

  ## TODO: redirect to plotFeatureBlocks if too many features

  ## only take types present in the current set
  ## NOTE: this was an if(missing(types))elseif before 20160212: 
  ## thus setting types overruled cuttypes: was this necessary sometimes?
  if ( cuttypes ) 
    types <- types[types%in%feat[,columns["type"]]]

  ## ylimits without typord
  typey <- NULL
  ylim <- c(-.1,1.1)
  ## generate y-order and ylim, unless y is set by types
  ## TODO: generate a matrix of feature overlaps to select y-value and ylim
  if ( typord & (nrow(coor)>0 & length(types)>0) ) {
    offset <- 1/length(types)
    typey <- 1-offset*0:(length(types)-1)
    names(typey) <- types
    ylim <- c(min(typey)-offset, max(typey)+offset)
  }
  ## plot
  #old.par <- par(no.readonly = TRUE)
  #on.exit(par(old.par))
  plot(0,xlim=xlim,ylim=ylim,col=NA,axes=FALSE,ylab=ylab,xlab=NA)
  offset <- - 1/nrow(feat)
  y <- .5
  for ( i in order(rowMeans(coor)) ) {
      
    typ <- as.character(feat[i,columns["type"]])
    if ( !typ %in% types ) next
    yt <- ifelse(typord, typey[typ], y) 
    col <- 1
    if ( "color"%in%names(columns) ) # take feature color if present
        if ( columns["color"] %in% colnames(feat) )
            col <- as.character(feat[i,columns["color"]])
    
    
    x <- coor[i,]
    if ( x[1]!=x[2] )
      arrows(x0=x[1], y0=yt, col=col,y1=yt, x1=x[2], length=.075, ...)
    points(x=x[1], y=yt, col=col, pch=3, ...)

      ## plot names
      if ( names ) {
          if ( feat[i,"insd"] ) x <- mean(x)
          else if ( feat[i,"span"] ) x <- mean(xlim)
          else if ( feat[i,"lovl"] ) x <- max(x)
          else if ( feat[i,"rovl"] ) x <- min(x)
          shadowtext(x=x,y=yt,labels=feat[i,columns["name"]],
                     col=col,bg="white",r=.15,cex=tcx,pos=tpos)
      }
    if ( typord ) next
    y <- y + offset
    if ( y > 0.9 & offset> 0 ) { offset <- -offset; y <- 1 }
    if ( y < 0.1 & offset< 0 ) { offset <- -offset; y <- 0 + offset/2 }
  }
  if ( typord & (nrow(coor)>0 & length(types)>0) ) 
    axis(2,at=typey,labels=names(typey),las=2,line=line,tcl=.3, mgp=c(2,.5,0),
         cex.axis=ycx)
  if ( axis1 ) axis(1,line=-1)

  ## silently return y-position of types (list of plotted features)
  features <- typey #rownames(feat[feat[,columns["type"]] %in% types,])
}

## plot 3D genomeData as a heatmap; kept private since official version
## of this function is maintained in package genomeBrowser
segment.plotHeat <- function(data, coors, orig.idx, breaks, colors,
                             orig.col, ylab="", axes=FALSE, axis2=FALSE,
                             colnorm=FALSE, chrS, coor, args) {
  
  ## parse arguments, which will also override other 
  if ( !missing(args) ) {    
    for ( i in seq_along(args) )       
      assign(names(args)[i],args[[i]])
  }

  chr <- coors[1]
  xlim <- sort(coors[2:3])

  ## split off coordinate columns, if not separately supplied
  ## TODO: use as alternative to chrS
  ## TODO: clean this up a bit (chrS vs. coor vs. data)
  firstcol <- 1
  if ( sum(c("chr","coor")%in%colnames(data))==2 )
    firstcol <- 3

  ## cut to xlim - TODO: consolidate with coor in plotFeatures and plotBlocks
  if ( missing(chrS) ) { ## search indices - SLOW
    idx <- which(as.numeric(data[,"chr"]) == chr &
                 (as.numeric(data[,"coor"]) >= xlim[1] &
                  as.numeric(data[,"coor"]) <= xlim[2]))
  } else { ## take direct index - ONLY FOR EXPANDED and ORDERED DATA
    start <- xlim[1]+chrS[chr]
    end <- xlim[2]+chrS[chr]
    if ( start < 1 ) start <- 1
    if ( end > nrow(data) ) end <- nrow(data) ## TODO: cut by chromosome!
    idx <- start:end
  }
  
  dat <- as.matrix(data[idx,,drop=FALSE])
  if ( !missing(orig.idx) ) orig.idx <- orig.idx[idx] 
  if ( missing(chrS) )
    x <- dat[,"coor"]
  else
    x <- xlim[1]:xlim[2] # start:end - chrS[chr]
  ##  x <- coor[idx,"coor"]
  dat <- as.matrix(dat[,firstcol:ncol(data),drop=FALSE])
  y <- seq_len(ncol(dat))

  if ( length(x)==0 | !any(!is.na(dat)) ) {  ## empty?
    plot(1,xlim=xlim,col=NA,axes=FALSE,xlab="chromosome position",ylab=ylab)
  }else{

    if ( colnorm ) { # between 0 and 1
      dat <- t(apply(dat,1,function(x) {
        if ( sum(is.na(x))==length(x) ) return(x)
        mx <- max(x,na.rm=TRUE)
        mn <- min(x,na.rm=TRUE)
        if ( mx==mn ) rep(0.5,length(x)) else ((x-mn)/(mx-mn))}))
    }
    
    if ( missing(breaks) ) {
      if ( !missing(colors) ) nbreaks <- length(colors)+1
      else nbreaks <- 101
      breaks <- seq(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE),length.out=nbreaks)
    }
    nbreaks <- length(breaks)
    if ( missing(colors) )
      colors <- gray(seq(1,0,length.out=nbreaks-1))
    
    
    p1 <- image(x=x,y=y,z=dat, axes=FALSE, ylab=ylab, xlim=xlim,
                xlab="chromosome position",breaks=breaks,col=colors)
  }
  ## if present/requested, plot locations of original data
  ## (or arbitrary discrete info on positions)
  if ( !missing(orig.idx) ) {
    if ( missing(orig.col) ) orig.col <- 2
    # cat(paste("plotting original data locations\n"))
    axis(1,at=x[which(orig.idx)],labels=NA,tcl=.2,
         col=NA,col.ticks=orig.col,lwd.ticks=1)
  }
  if ( axes ) {
    axis(1)
    axis(2,las=2)
  } else if ( axis2 )
    axis(2,las=2)
}

### PLOTTING SEGMENTATION INTERNAL STRUCTURES 




## plot a matrix as seen in R; kept private since official version
## of this function is maintained in package segmenTools
image_matrix <- function(dat, text, text.col,
                         axis=1:2, axis1.col, axis2.col, ...) {

    ## reverse columns and transpose
    if ( nrow(dat)>1 )
        imgdat <- t(apply(dat, 2, rev))
    else
        imgdat <- t(dat)
    image(x=seq_len(ncol(dat)), y=seq_len(nrow(dat)), z=imgdat, axes=FALSE, ...)

    ## add text
    if ( !missing(text) ) {
        if ( missing(text.col) )
            text.col <- rep(1, length(c(text)))
        text(x=rep(seq_len(ncol(dat)),nrow(dat)),
             y=rep(rev(seq_len(nrow(dat))),each=ncol(dat)),
             paste(t(text)),col=t(text.col))
        }
    
    ## add axes
    ## TODO : handle axes=FALSE
    if ( !missing(axis) ) {
        if ( 1 %in% axis ) 
            if ( !missing(axis1.col) ) # colored ticks
                for ( i in seq_len(ncol(dat)) )
                    axis(1, at=i,colnames(dat)[i],
                         col.axis=axis1.col[i], col=axis1.col[i],
                         las=2, cex.axis=1.5, lwd=2)
            else
                axis(1, at=seq_len(ncol(dat)), labels=colnames(dat), las=2)
        if ( 2 %in% axis )
            if ( !missing(axis2.col) ) # colored ticks
                for ( i in seq_len(nrow(dat)) )
                        axis(2, at=nrow(dat)-i+1, rownames(dat)[i],
                             col.axis=axis2.col[i], col=axis2.col[i],
                             las=2, cex.axis=1.5, lwd=2)
            else
                axis(2, at=rev(seq_len(nrow(dat))), rownames(dat),las=2)        
    }
} 


#' Plot method for the "timeseries" object.
#' 
#' plot the processed time-series object returned from
#' \code{\link{processTimeseries}}.
#' @param x a time-series object as returned by
#' \code{\link{processTimeseries}}
#' @param plot a string vector indicating the values to be plotted;
#' "total": plot of the total signal, summed over
#' the time-points, and indicating the applied threshold \code{low.thresh};
#' note that the total levels may have been transformed (e.g. by
#' \code{\link{log_1}} or \code{\link{ash}}) depending on the arguments
#' \code{trafo} and \code{dc.trafo} in \code{\link{processTimeseries}};
#' "timeseries": plot the complete time-series as a heatmap, where time is
#' plotted bottom-up on the y-axis and segmentation coordinates on the x-axis;
## "clusters": this requires
## a clustering set from clusterTimeseries and plots average time-series
## for clusters, independent of their genomic coordinate
#' @param xaxis x-values to use as x-axis (e.g. to reflect absolute
#' chromosomal coordinates)
#' @param ylabh plot y-axis title horizontally
## @param cset a clustering object for the time series in \code{x} from
## \code{\link{clusterTimeseries}}, only required for option `clusters'
## in argument \code{plot}
## @param K a cluster number K to use in clustered time series plots,
## only required for option `clusters' in argument \code{plot},
#' @param ... additional arguments to plot of total signal
#'@export
plot.timeseries <- function(x, plot=c("total","timeseries"),
                            xaxis, ylabh=TRUE, ...) {
    
    tset <- x
    
    ## get time-series data
    ts <- tset$ts # incl. all trafos and zeros set to NA
    ts[tset$zero.vals,] <- NA
    
    tot <- tset$tot # total of the time-series

    ## mock "chromosome" coordinates
    rel.coors <- c(chr=1,start=1,end=nrow(ts)) 
    if ( missing(xaxis) )
      coors <- rel.coors
    else
      coors <- c(chr=1,start=min(xaxis), end=max(xaxis))

    ## timeseries heatmap colors
    colors0 <- gray(seq(1,0,length.out=100))

    ## settings
    low.thresh <- tset$settings$low.thresh
    logged <- tset$settings$trafo != "raw" # plot with log if not already done
    
    if ( "total" %in% plot ) {
        plot(coors["start"]:coors["end"],tot,log=ifelse(logged,"","y"),
             type="l",lwd=2,axes=FALSE,ylab=NA,xlab=NA, ...)
        graphics::polygon(x=c(coors["start"],coors["start"],
                              coors["end"],coors["end"]),
                          y=c(min(tot,na.rm=TRUE),rep(low.thresh,2),
                              min(tot,na.rm=TRUE)),col="#00000055",border=NA)
        graphics::abline(h=low.thresh,col="#000000BB")
        lines(coors["start"]:coors["end"],tot)
        axis(2, las=2)
        #axis(1)
        if ( ylabh )
            graphics::mtext("total signal", 2, 4, las=2, cex=1.25)
        else
            graphics::mtext("total signal", 2, 2)
    }
    if ( "timeseries" %in% plot ) {
        segment.plotHeat(ts,coors=rel.coors,chrS=0,colors=colors0, colnorm=TRUE)
        axis(2,at=seq_len(ncol(ts)))
        #axis(1)
        if ( ylabh )
            graphics::mtext("time points", 2, 4, las=2, cex=1.25)
        else
            graphics::mtext("time points", 2, 2)
    }
    ## TODO: provide this for function plot.clustering instead!
    #if ( "clusters" %in% plot ) {
    #    if ( missing(cset) )
    #      stop("Option `clusters' in plot requires a clustering;",
    #           "see argument `cset'")
    #    if ( missing(K) )
    #      K <- cset$selected
    #    if ( is.null(K) | ! K %in% colnames(cset$clusters) )
    #      stop("No cluster selection found; use argument `selected'")
    #    avg <- getClsAvg(x, cset$clusters[,K], q=.1)
    #    plotClsAvg(avg, col=cset$colors[[K]])
    #    
    #}
}


#' Plot method for the "clustering" object.
#' 
#' plot the clustering object returned by \code{\link{clusterTimeseries}}
#' @param x  a "clustering" object as returned by
#' \code{\link{clusterTimeseries}}
#' @param k a numeric or string vector indicating the clusterings to be plotted;
#' specifically the column numbers or names in the matrix of clusterings
#' in \code{cset$clusters}; if missing all columns will be plotted
#' and the calling code must take care of properly assigning \code{par(mfcol)}
#' or \code{layout} for the plot
#' @param sort if \code{TRUE} and the clustering is yet unsorted a cluster
#' sorting will be calculated based on "ccor" cluster-cluster similarity
#' matrix \code{x$Ccc}; see \code{\link{sortClusters}}
#' @param xaxis optionally x-values to use as x-axis (e.g. to reflect absolute
#' chromosomal coordinates)
#' @param axes list of axes to plot, numbers as used as first argument
#' in function \code{axis}
#' @param pch argument \code{pch} (symbol) for plot
#' @param ylabh plot "clustering" horizontally at y-axis
#' @param ... additional arguments to plot, eg. to set point \code{cex}
#' @return returns the input "clustering" object with (potentially new)
#' cluster sorting and colors as in shown in the plot
#'@export
plot.clustering <- function(x, k, sort=FALSE, xaxis,
                            axes=1:2, pch=16, ylabh=TRUE, ...) {

    cset <- x
    
    if ( missing(k) ) 
        k <- seq_len(ncol(cset$clusters))

    ## cluster sorting via Ccc (cluster-cluster correlation)
    if ( !"sorting" %in% names(cset) )
        cset <- sortClusters(cset, sort=sort, verb=1)

    ## cluster colors
    if ( !"colors" %in% names(cset) )
      cset <- colorClusters(cset)

    ## plotting all: layout or mfcol must be set from outside;
    ## or a specific k chosen to only plot the k'th clustering
    if ( length(k)>1 )
        par(mfcol=c(length(k),1))
    for ( i in k ) {
        seq <- as.character(cset$clusters[,i])
        cls.srt <- cset$sorting[[i]]
        if ( missing(xaxis) )
          xaxis <- seq_along(seq)
        y <- seq_along(cls.srt)
        names(y) <- cls.srt
        cols <- cset$colors[[i]]
        ## plot original clustering
        plot(xaxis,y[seq],axes=FALSE,xlab="",ylab=NA,
             col=cols[seq], pch=pch, ...)
        if ( 2%in%axes )
          axis(2, at=y, labels=names(y), las=2)
        if ( 1%in% axes )
          axis(1,labels=FALSE)
        if ( 3%in%axes )
          axis(3,labels=FALSE)
        if ( 4%in%axes)
          axis(4, at=y, labels=names(y), las=2)
        graphics::mtext(colnames(cset$clusters)[i], 2, 2)
        if ( ylabh )
            graphics::mtext("clustering", side=2, line=4.5, las=2, cex=1.25)
    }
    silent <- cset # silent return of cset with sorting and clustering
}

#' Plot method for the "segments" object.
#' 
#' plot the final segmentation objects returned by
#' \code{\link{segmentClusters}} and \code{\link{segmentCluster.batch}}
#' @param x a segmentation object as returned by
#' \code{\link{segmentClusters}} and \code{\link{segmentCluster.batch}}
#' @param params a named vector of parameter settings used in
#' \code{\link{segmentCluster.batch}} allows to filter plotted segment
#' types, e.g. params=c(S="icor") will only plot segments where
#' the scoring function (parameter S) "icor" was used.
#' @param types a string vector indicating segment types to plot (a subset of
#' \code{x$ids}; defaults to all in \code{x$ids})
#' @param xaxis optional x-values to use as x-axis (e.g. to reflect absolute
#' chromosomal coordinates)
#' @param show.fused show the fuse tag as a black x
#' @param plot string vector indicating which data should be plotted;
#' "segments": plot segments as arrows; "S1" plot the scoring vectors
#' \code{s(i,j,c} for all \code{c}; "S" plot the derivative of
#' matrix \code{S(i,c)} for all \code{c}
#' @param ... additional arguments forwarded to
#' \code{\link[graphics:arrows]{arrows}}, eg. to set \code{lwd} for
#' for \code{plot="segments"}, or to \code{\link[graphics:matplot]{matplot}}
#' for \code{plot="S"}
#'@export
plot.segments <- function(x, plot=c("S","segments"),  types, params,
                          xaxis, show.fused=FALSE, ...) {
    for ( pl in plot ) 
        plotSegments(x, plot=pl, types=types, params=params,
                     xaxis=xaxis, show.fused=show.fused, ...)
}
## above public wrapper allows to sort the plots by the order in plots
plotSegments <- function(x, plot=c("segments", "S", "S1"), types, params,
                         xaxis, show.fused=FALSE, ...) {
    sset <- x

    ## sub-types to be plotted
    if ( missing(types) ) types <- NULL
    if ( is.null(types) ) {
        if ( !missing(params) ) {
            if ( "settings"%in% names(sset) ) { 
                for ( i in seq_along(params) ) {
                    pid <- names(params)[i]
                    p <- params[i]
                    typ <- sset$settings[,pid] == p
                    types <- c(types, rownames(sset$settings)[typ])
                }
            }
        } else if ( "ids" %in% names(sset) ) 
          types <- sset$ids # rownames(sset$settings)
        else
          types <- "segments" # for direct function segmentClusters
    }
    ## type column not present from direct function segmentClusters
    if ( !"type"%in% colnames(sset$segments) ) 
        sset$segments <- cbind(sset$segments,
                               type=rep("segments",nrow(sset$segments)))
    
    
    ## one plot for all segments 
    if ( "segments" %in% plot ) {

        ## mock "chromosome" coordinates
        if ( missing(xaxis) )
          coors <- c(chr=1,start=1,end=sset$N) 
        else
          coors <- c(chr=1,start=min(xaxis), end=max(xaxis))

        ## get & process segments
        ## add colors!
        segs <- sset$segments
        if ( is.null(segs) ) {
            plot(1,1,col=NA,axes=FALSE,ylab=NA,xlab=NA)
            text(1,1,"no segments",cex=2)
        } else {        
            if ( !"color" %in% colnames(segs) ) {
                segs <- cbind(segs, color=rep(NA, nrow(segs)))
                for ( typ in unique(segs[,"type"]) ) {
                    cols <- sset$colors[[typ]]
                    segs[segs[,"type"]==typ,"color"] <-
                        cols[as.character(segs[segs[,"type"]==typ,"CL"])]
                }
            }

            ## column mapping required for segment.plotFeatures
            columns <- c(type="type", start="start", end="end",
                         color="color") # name="ID" required for names=T

            ## filter allsegs by segments for the current clustering
            ypos <- segment.plotFeatures(segs, types=types,
                                         coors=coors, typord=TRUE,cuttypes=TRUE,
                                         ylab="", names=FALSE,columns=columns,
                                         tcx=.5, ...)
            ##axis(1)
            ## plot fuse tag - only present in segments from batch function
            if ( show.fused & "fuse" %in% colnames(segs) ) {
                fuse <- segs[segs[,"fuse"],]
                points(fuse[,"start"], ypos[fuse[,"type"]], col="black",
                       pch=4, lwd=1, cex=1.5)
            }
        }
    }
    
    ## colors for S1 heatmap
    colors0 <- gray(seq(1,0,length.out=100))
    
    if ( any(c("S","S1") %in% plot) ) {

        if ( !"SK"%in%names(sset) )
          stop("ERROR: matrix plot was requested but matrices where",
               "not stored. Use `save.matrix=TRUE' in segmentation or",
               "omit 'S' and 'S1' from  argument `plot' in plot")
        
        ## get matrices, cluster sorting and colors
        SK <- sset$SK[types]
        if ( "sorting" %in% names(sset) )
            sk.srt <- sset$sorting[types]
        else
            sk.srt <- lapply(SK, function(x) colnames(x$S))
        sk.col <- sset$colors[types]

        ## one plot for each segmentation!
        for ( j in seq_along(SK) ) {

            ## cluster sorting
            srt <- as.character(sk.srt[[j]])
            if ( !"0" %in% srt ) # in generated sk.srt nuisance "0" is present
                srt <- c("0",srt)

            ## plot S1 as heatmap
            if ( "S1" %in% plot ) {
                S1 <- t(SK[[j]]$S1[,rev(srt)]) # sort by cluster sorting!
                S1 <- S1/apply(S1,1,mean) # TODO: is colMeans most informative?
                ## TODO: add from segmenTools
                image_matrix(S1,axis=2, col=colors0, ylab=expression(S[1](i,c)))
                graphics::mtext(names(SK)[j], side=2 , line=4.5, las=2)
            }
            if ( !"S" %in% plot ) next
            
            ## cluster coloring; add black for nuisance
            sgcols <- sk.col[[j]][srt]
            ## add alpha
            sgcols <- paste(sgcols,"EE",sep="") 
            names(sgcols) <- srt
            ## set NA: only show clusters that actually produced a segment
            tp <- sset$segments[,"type"]%in%names(SK)[j]
            sgcols[!names(sgcols)%in%as.character(sset$segments[tp,"CL"])] <- NA

            ## get matrix and sort according to cluster sorting
            ## TODO: plot by segment; highlight winning segment!!
            S <- SK[[j]]$S[,srt, drop=FALSE]
            dS <- apply(S,2,function(x) c(0,diff(x)))

            ## x-axis: x can be passed to use real coordinates
            if ( missing(xaxis) )
              xaxis <- seq_len(nrow(S))
            xlim <- range(xaxis)             

            xrng <- stats::quantile(xaxis,c(.05,.95), na.rm=TRUE)
            xidx <- which(xaxis>xrng[1]&xaxis<xrng[2]) #x%in%xrng[1]:xrng[2]
            ylim <- stats::quantile(ash(dS[xidx,]),c(0,1), na.rm=TRUE)
            plot(-1,col=NA,ylim=ylim,xlim=xlim,
                 ylab=expression(ash(Delta~S(i,c))))
            lines(xaxis,ash(dS[,1]),lwd=7,col="#00000015") # NUI: BACKGROUND 
            lines(xaxis,ash(dS[,1]),lwd=1,lty=3,col="#00000099") 
            graphics::matplot(xaxis, ash(dS), type="l",
                              add=TRUE, col=sgcols, lty=1, ...)
            graphics::mtext(names(SK)[j], side=2 , line=4.5, las=2)
        }
    }
}

#' Summary plot for the \code{segmenTier} pipeline.
#' 
#' Plot all objects from the segmentation pipeline, i.e. the processed
#' time-series, the clustering, the internal scoring matrices and
#' the final segments.
#' @param tset a time-series object as returned by
#' \code{\link{processTimeseries}}
#' @param cset a clusterings object as returned by
#' \code{\link{clusterTimeseries}}
#' @param sset a segmentation object as returned by
#' \code{\link{segmentClusters}} and \code{\link{segmentCluster.batch}}
#' @param split split segment plots by clustering plots
#' @param plot.matrix include the internal scoring matrices in the plot
#' @param mai margins of individual plots, see \code{par}
#' @param ... further arguments to
#' \code{\link{plot.clustering}} (\code{cset}) and
#' \code{\link{plot.segments}} (\code{sset}). Note: these may
#' conflict and cause errors, but eg. a combination of \code{cex=0.5, lwd=3}
#' works, affecting cluster point size and segment line width, respectively.
#'@export
plotSegmentation <- function(tset, cset, sset, split=FALSE, plot.matrix=FALSE,
                             mai=c(.01,2,.01,.01), ...) {

    nt <- ifelse(is.null(tset), 0, 2)
    nsg <- length(sset$ids)# total number of segmentations
    nk <- length(cset$ids) # number of clusterings
    spk <- nsg/nk # segmentations per clustering
    ## number of plots
    ## each clustering can have multiple segmentations; plot each
    ## 2 for time-series; and for each clustering 2 (clustering and segments),
    ## plus S1 & S 
    nplots <- nt + nk * (ifelse(split,2,1) + ifelse(plot.matrix, 2*spk, 0)) +
      ifelse(split,0,1)
    
    old.par <- par(c("mai","xaxs","mfcol")) # store plot params and set new
    par(mfcol=c(nplots,1), xaxs="i", mai=mai)

    ## TIME-SERIES PLOT UTILITY: plot both the total signal (optionally used
    ## for threshold) and a heatmap of the time-series
    if ( !is.null(tset) ) {
        plot(tset, plot=c("total","timeseries"))
        axis(1)
    }
    
    ## CLUSTERING PLOT UTILITY: 
    ## NOTE that clusterings are sorted (by their similarity matrix `Ccc`)
    ## and colored along a color-wheel
    if ( !split ) {
        for ( k in seq_len(ncol(cset$clusters)) )
            plot(cset, k, ...)
        plot <- "segments"
        if ( plot.matrix ) plot <- c("segments","S","S1")
        plot(sset, plot=plot, ...)
    } else
        for ( k in seq_len(ncol(cset$clusters)) ) {
            ## plot clustering
            plot(cset, k, ...)
            ## SEGMENTATION PLOT UTILITY
            ## plot all segments, S1(i,c), S(i,c) for this clustering
            kid <- cset$ids[k]
            types <- rownames(sset$settings)[sset$settings[,"K"] %in% kid]
            plot <- "segments"
            if ( plot.matrix ) plot <- c("segments","S","S1")
            plot(sset, plot=plot, types=types, ...) 
        }
    par(old.par)
}
