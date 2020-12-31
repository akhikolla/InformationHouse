# ---------------
# pgrid
# ---------------
# 
# Constructor
#
pgrid <- function(mass, boundary, gridtriple, generator, structure) {
	
  stopifnot(is.array(mass))   # matrices and 2-dim arrays seem to be exactly the same (incl. classes)
  dimension <- length(dim(mass))
  n <- dim(mass)
  N <- prod(n)
  if (missing(structure)) {
  	if (length(unique(dim(mass)) == 1)) {
  	  structure <- "square"
  	} else {
  	  structure <- "rectangular"
  	}
  }	
  nmiss <- missing(boundary) + missing(gridtriple) + missing(generator)
  if (nmiss == 3) {
  	  # changed 12/10/16: always assume square pixels and first boundaries c(0,1) if nothing is given
  	  boundary <- rep(0,2*dimension)
  	  boundary[seq(2,2*dimension,2)] <- n/n[1]
  	  gridtriple <- t(sapply(n, function(m) {c(1/(2 * n[1]), m/n[1] - 1/(2 * n[1]), 1/n[1])}))
      generator <- lapply(n, function(m) {seq(1/(2 * n[1]), m/n[1] - 1/(2 * n[1]), 1/n[1])})
  } else if (nmiss <= 1) {
  	  stop("only one of 'boundary', 'gridtriple', and 'generator' may be specified.")
  } else if (!missing(boundary)) {
  	  stopifnot(length(boundary) == 2*dimension)
  	  even <- seq(2,2*dimension,2)
      odd <- seq(1,2*dimension-1,2)
  	  halfinter <- (boundary[even]-boundary[odd])/(2*n)
  	  gridtriple <- cbind(boundary[odd]+halfinter, boundary[even]-halfinter, 2*halfinter)
  	  generator <- lapply(1:dimension, function(i) {seq(gridtriple[i,1], gridtriple[i,2], gridtriple[i,3])})
  } else if (!missing(gridtriple)) {
  	  if (is.vector(gridtriple)) {
  	  	gridtriple <- matrix(gridtriple, dimension, length(gridtriple), byrow=TRUE)
  	  }
  	  stopifnot(all(dim(gridtriple) == c(dimension, 3)))
  	  stopifnot(isTRUE(all.equal(gridtriple[,1]+(n-1)*gridtriple[,3], gridtriple[,2])))
  	  boundarymat <- cbind(gridtriple[,1] - gridtriple[,3]/2, gridtriple[,2] + gridtriple[,3]/2)
  	  boundary <- as.vector(t(boundarymat))
  	  generator <- lapply(1:dimension, function(i) {seq(gridtriple[i,1], gridtriple[i,2], gridtriple[i,3])})
  } else {  # (!missing(generator))
  	  if (!is.list(generator)) {
  	  	generator <- lapply(1:dimension, function(x) {return(generator)})
  	  }
      stopifnot(is.list(generator) && length(generator) == dimension && all(sapply(generator, length) == n))
      temp <- lapply(generator, function(x) {diff(diff(x))})
      stopifnot(all(sapply(1:dimension, function(i) {isTRUE(all.equal(temp[[i]], rep(0,n[i]-2)))})))
      gridtriple <- cbind(sapply(generator, function(x) {x[1]}),
                          sapply(1:dimension, function(i) {generator[[i]][n[i]]}),
                          sapply(generator, function(x) {mean(diff(x))}))     
  	  boundarymat <- cbind(gridtriple[,1] - gridtriple[,3]/2, gridtriple[,2] + gridtriple[,3]/2)
  	  boundary <- round(as.vector(t(boundarymat)), 12)
  }  	
  
  pixelarea <- prod(gridtriple[,3])
  totmass <- sum(mass)     # sum of masses at each pixel
  totcontmass <- totmass*pixelarea   # mass interpreted as uniformely distributed over pixel area
  res <- list(structure=structure, dimension=dimension, n=n, N=N, boundary=boundary, gridtriple=gridtriple,
              generator=generator, mass=mass, totmass=totmass, totcontmass=totcontmass)
  class(res) <- "pgrid"
  return(res)
}


# 
# plot method
# (this plots one pgrid object, two pgrid objects (next to each other or as diff), or two pgrid objects
# as diff and their transportation plan)
# it also calls plot_pgrid_wpp for plotting semi-discrete transports
# diffpic only does something for two pgrid objects without transportation plan
# rot=FALSE uses the usual R image-convention
# rot=TRUE plots mass-matrices the same way as they are displayed in numeric output and adapts the tplan-arrows accordingly
plot.pgrid <- function(x,y=NULL,tplan=NULL,mass=c("colour","thickness"),length=0.1,acol,bcol=4,lwd,rot=FALSE,overlay=FALSE,static.mass=TRUE,...) {
  stopifnot(class(x) == "pgrid")
  if (class(y) == "wpp") {
  	if (missing(acol)) { acol <- "#996699" }
  	if (missing(lwd)) { lwd <- 1.5 }
  	return(plot_pgrid_wpp(x,y,tplan,pmass=TRUE,cex=0.8,length=length,acol=acol,col=bcol,lwd=lwd,rot=TRUE,...))
  }
  a <- x
  if (a$dimension != 2) stop("plot.pgrid is currently only implemented for 2-d grids")
  xi <- a$generator[[1]]
  eta <- a$generator[[2]]  
  if (missing(y)) {
  	image2(xi, eta, a$mass, rot=rot, col=grey(0:200/200), asp=1, xlab="", ylab="")
  	invisible()
  } else {
  	stopifnot(class(y) == "pgrid")
    b <- y
    if (b$dimension != 2) stop("plot.pgrid is currently only implemented for 2-d grids")
    if (!(a$structure %in% c("square", "rectangular")) || !(b$structure %in% c("square", "rectangular")))
    stop("transport.pgrid is currently only implemented for rectangular pixel grids")
#    ngrid <- a$n
#    Ngrid <- a$N
    if (missing(tplan)) {
      if (!overlay) {
      	par(mfrow=c(1,2))
   	    image2(xi, eta, a$mass, rot=rot, col=grey(0:200/200), asp=1, xlab="", ylab="", ...)
   	    image2(b$generator[[1]], b$generator[[2]], b$mass, rot=rot, col=grey(0:200/200), asp=1, xlab="", ylab="", ...)
   	    par(mfrow=c(1,1))      	
      } else {
      	stopifnot(compatible(a,b))
      	zeta <- a$mass-b$mass
      	image2(xi,eta,zeta, rot=rot, col=grey(0:200/200), asp=1, xlab="", ylab="", ...)
      }
      invisible()
    } else {
      stopifnot(compatible(a,b))
      mass <- match.arg(mass)
      if (missing(acol))
        acol <- 2
      if (missing(lwd))
        lwd <- 3    
      gg <- expand.grid(xi,eta)
      maxmass <- max(tplan[,3])	
      zeta <- a$mass-b$mass
      
      if (rot) {
      	yco <- rev(gg[,1])
        xco <- gg[,2]
      } else {
        xco <- gg[,1]
        yco <- gg[,2]
      }
      image2(xi,eta,zeta, rot=rot, col=grey(0:200/200), asp=1, axes=FALSE, xlab="", ylab="", ...)
      if (mass == "colour") {
        stplan <- tplan[order(tplan[,3]),]	
        cc <- heat.colors(128)
        arrcols <- cc[as.numeric( as.character(cut(stplan[,3], breaks=seq(0,maxmass,length.out=129), labels=1:128)) )]
        wh <- which(stplan$from != stplan$to)
        arrows(xco[stplan$from[wh]],yco[stplan$from[wh]],xco[stplan$to[wh]],yco[stplan$to[wh]],
           angle=5, length, col=arrcols[wh], lwd=lwd)
        if (static.mass == TRUE) {   
          nwh <- (1:dim(stplan)[1])[-wh]
          points(xco[stplan$from[nwh]],yco[stplan$from[nwh]], pch=16, cex=0.5+lwd*0.1, col=arrcols[nwh])
          points(xco[stplan$from[nwh]],yco[stplan$from[nwh]], cex=0.5+lwd*0.1)
        }
      } else {
      	wh <- which(tplan$from != tplan$to)
        arrows(xco[tplan$from[wh]],yco[tplan$from[wh]],xco[tplan$to[wh]],yco[tplan$to[wh]],
           angle=5, length, col=acol, lwd=(8*tplan[,3]/maxmass)[wh])
        if (static.mass == TRUE) {   
          nwh <- (1:dim(tplan)[1])[-wh]
          points(xco[tplan$from[nwh]],yco[tplan$from[nwh]], pch=16, cex=0.5 + (8*tplan[,3]/maxmass)[nwh] * 0.1, col=acol)
          points(xco[tplan$from[nwh]],yco[tplan$from[nwh]], cex=0.5 + (8*tplan[,3]/maxmass)[nwh] * 0.1)
        }
      }
      invisible()
    }
  }
}


image2 <- function(x,y,z,rot=FALSE,...) {
  rotclock <- function(m) t(m)[,nrow(m):1]	
  if (rot) {
  	image(y,x,rotclock(z),...)
  } else {
  	image(x,y,z,...)
  }
}


image3 <- function(z,x=1:dim(z)[1],y=1:dim(z)[2],rot=TRUE,...) {
  rotclock <- function(m) t(m)[,nrow(m):1]	
  if (rot) {
  	image(y,x,rotclock(z),...)
  } else {
    image(x,y,z,...)
  }
}


print.pgrid <- function(x, ...) {
  stopifnot(class(x) == "pgrid")
  if (x$dimension == 2) {
  	cat("Regularly spaced ",x$n[1],"x",x$n[2]," pixel grid on [",
  	  x$boundary[1],",",x$boundary[2],"] x [",x$boundary[3],",",x$boundary[4],"].\n",sep="")
  	cat("x-gridtriple:", x$gridtriple[1,], "\n")
  	cat("y-gridtriple:", x$gridtriple[2,], "\n")
  	cat("pixel masses range from ", min(x$mass), " to ", max(x$mass), "\n", sep="")
  	cat("total pixel mass: ", x$totmass, "\n", sep="")
  	cat("total continuum mass: ", x$totcontmass, "\n", sep="")    	
  } else
  if (x$dimension == 3) {
  	cat("Regularly spaced ",x$n[1],"x",x$n[2],"x",x$n[3]," pixel grid on [",
  	  x$boundary[1],",",x$boundary[2],"] x [", x$boundary[3],",",x$boundary[4],
  	  "] x [",x$boundary[5],",",x$boundary[6],"].\n",sep="")
  	cat("x-gridtriple:", x$gridtriple[1,], "\n")
  	cat("y-gridtriple:", x$gridtriple[2,], "\n")
  	cat("z-gridtriple:", x$gridtriple[3,], "\n")
  	cat("pixel masses range from ", min(x$mass), " to ", max(x$mass), "\n", sep="")
  	cat("total pixel mass: ", x$totmass, "\n", sep="")
  	cat("total continuum mass: ", x$totcontmass, "\n", sep="")    	
  } else {
  	cat("Regularly spaced pixel grid in ", x$dimension, " dimensions.\n", sep="")
    cat("number of pixels in each dimension:", x$n, "\n")
    cat("bounding box:", x$boundary, "\n")
    cat("gridtriples:\n")
      for (i in 1:x$dimension) cat("", x$gridtriple[i,], "\n")
  	cat("pixel masses range from ", min(x$mass), " to ", max(x$mass), "\n", sep="")
  	cat("total pixel mass: ", x$totmass, "\n", sep="")
  	cat("total continuum mass: ", x$totcontmass, "\n", sep="")    	
  }
  invisible(x)
}


summary.pgrid <- function(object, ...) {
  print.pgrid(object)
}


# ---------------
# pp  (this class is not *really* useful so far)
# Later: include different masses at points
# ---------------
# 
# Constructor
#
pp <- function(coordinates) {
  coordinates <- as.matrix(coordinates)
  dimension <- dim(coordinates)[2]
  N <- dim(coordinates)[1]
  res <- list(dimension=dimension, N=N, coordinates=coordinates)
  class(res) <- "pp"
  return(res)  
}

# 
# plot method
#
plot.pp <- function(x,y=NULL,tplan=NULL,cols=c(4,2),cex=0.8,acol=grey(0.3),lwd=1,overlay=TRUE,...) {
  stopifnot(class(x) == "pp" && class(x) == "pp")  
  if (x$dimension != 2) stop("plot.pp is currently only implemented for 2-d point patterns")  
  oldpars <- par(xpd=TRUE, xaxs="i", yaxs="i")
  if (missing(y)) {
  	# if (missing(cols)) { cols <- 1 }
  	plot(x$coordinates, axes=FALSE, xlab="", ylab="", col=cols[1], pch=16, cex=cex, asp=1, ...)
  	invisible()
  } else {
  	if (y$dimension != 2) stop("plot.pp is currently only implemented for 2-d point patterns")
  	min1 <- min(x$coordinates[,1],y$coordinates[,1])
    max1 <- max(x$coordinates[,1],y$coordinates[,1])
    min2 <- min(x$coordinates[,2],y$coordinates[,2])
    max2 <- max(x$coordinates[,2],y$coordinates[,2])
  	if (missing(tplan)) {
      if (!overlay) {
      	par(mfrow=c(1,2))
      	plot(x$coordinates, xlim=c(min1,max1), ylim=c(min2,max2), axes=FALSE, xlab="", ylab="", col=cols[1], pch=16, cex=cex, asp=1, ...)
      	plot(y$coordinates, xlim=c(min1,max1), ylim=c(min2,max2), axes=FALSE, xlab="", ylab="", col=cols[1], pch=16, cex=cex, asp=1, ...)
   	    par(mfrow=c(1,1))      	
      } else {
      	plot(x$coordinates, xlim=c(min1,max1), ylim=c(min2,max2), axes=FALSE, xlab="", ylab="", col=cols[1], pch=16, cex=cex, asp=1, ...)
        points(y$coordinates, col=cols[2], pch=16, cex=cex)
      }
      invisible()
  	} else {
  	  # first empty plot because we want transport arrows be covered by the points
  	  plot(NULL,type="n", xlim=c(min1,max1), ylim=c(min2,max2), axes=FALSE, xlab="", ylab="", xaxs="i", asp=1, ...)
      segments(x$coordinates[tplan$from,1], x$coordinates[tplan$from,2], y$coordinates[tplan$to,1], y$coordinates[tplan$to,2], col=acol, lwd=lwd)
      points(x$coordinates, col=cols[1], pch=16, cex=cex)
      points(y$coordinates, col=cols[2], pch=16, cex=cex)
      invisible()
    }
  }
  par(oldpars)
}


print.pp <- function(x, ...) {
  stopifnot(class(x) == "pp")
  cat("Pattern of ",x$N," points in ",x$dimension, " dimensions.\n", sep="")
  cat("Minimal coordinates:", apply(x$coordinates,2,min), "\n")
  cat("Maximal coordinates:", apply(x$coordinates,2,max), "\n")
  invisible(x)
}


summary.pp <- function(object, ...) {
  print.pp(object)
}



# ---------------
# wpp
# ---------------
# 
# Constructor
#
wpp <- function(coordinates, mass){
  coordinates <- as.matrix(coordinates)
  NN <- dim(coordinates)[1]
  stopifnot(length(mass) == NN)
  dimension <- dim(coordinates)[2]
  totmass <- sum(mass)
  stopifnot(all(mass >= rep(0,NN)))
  stopifnot(totmass > 0)
  massnonzero <- (mass != 0)
  mass <- mass[massnonzero]
  zeropoints <- coordinates[!massnonzero,]
  coordinates <- coordinates[massnonzero,]
  N <- dim(coordinates)[1]
  res <- list(dimension = dimension, N=N, coordinates = coordinates, mass = mass, totmass = totmass)
  class(res) <- "wpp"
  if (N < NN) { 
    attr(res, "zeropoints") <- zeropoints 
  }
  return(res)
}


plot.wpp <- function(x,y=NULL,tplan=NULL,pmass=TRUE,tmass=TRUE,cols=c(4,2),cex=0.8,aglevel=0.4,acol=grey(0.3),lwd=1,overlay=TRUE,...) {
  stopifnot(class(x) == "wpp")  
  if (x$dimension != 2) stop("plot.wpp is currently only implemented for 2-d point patterns")  
  oldpars <- par(xpd=TRUE, xaxs="i", yaxs="i")
#  col1 <- heat.colors(128)[1:60]
#  xorder <- order(x$mass)
#  col1 <- col1[as.numeric(as.character(cut(x$mass,breaks=seq(0,max(x$mass),length.out =61),labels=1:60)))]
  if (missing(y)) {
  	if (pmass) {
  	  massmean <- mean(x$mass)
  	  cexfac <- sqrt(x$mass/massmean)
  	  toosmall <- (cexfac < 0.25)
  	  toolarge <- (cexfac > 5)
  	  justright <- !(toosmall | toolarge)	
  	  cexfac <- cexfac[justright]  
  	  plot(x$coordinates, type="n", xlab="", ylab="", asp=1, ...)	
  	  # drop im folgenden wichtig, sonst werden bei nur einer Zeile Punkte (c_11, 1) und (c_12, 2) geplottet 
  	  points(x$coordinates[toolarge,,drop=FALSE], col=cols[1], pch=16, cex=5*cex) 
  	  points(x$coordinates[toolarge,,drop=FALSE], col=grey(1), pch=10, cex=5*cex, lwd=2)
  	  # black contours vor better visibility if points overlapp:
  	  # points(x$coordinates[toolarge,], col=1, pch=1, cex=0.92*5*cex)
  	  points(x$coordinates[justright,,drop=FALSE], col=cols[1], pch=16, cex=cexfac*cex)
  	  # points(x$coordinates[justright,], col=1, pch=1, cex=0.92*cexfac*cex)
  	  points(x$coordinates[toosmall,,drop=FALSE], col=cols[1], pch="+", cex=0.4*cex)  	  
    } else {
      plot(x$coordinates, xlab="", ylab="", col=cols[1], pch=16, cex=cex, asp=1, ...)
    }
    invisible()
  } else {
    stopifnot(class(y) == "wpp")
    if (y$dimension != 2) stop("plot.wpp is currently only implemented for 2-d point patterns")
    
    if (missing(tplan)) {
      if (!overlay) {
        par(mfrow=c(1,2))
        plot(x=x,tplan=NULL,pmass=pmass,cols=cols[1],cex=0.8,lwd=lwd,overlay=FALSE,...)
        plot(x=y,tplan=NULL,pmass=pmass,cols=cols[2],cex=0.8,lwd=lwd,overlay=FALSE,...)
        par(mfrow=c(1,1))      	
      } else {
      	allcoords <- rbind(x$coordinates,y$coordinates)  
  	    plot(allcoords, type="n", xlab="", ylab="", asp=1, ...)	
  	    if (pmass) {
  	      massmean <- mean(c(x$mass,y$mass))
  	      # First pp:
  	      cexfac <- sqrt(x$mass/massmean)
  	      toosmall <- (cexfac < 0.25)
  	      toolarge <- (cexfac > 5)
  	      justright <- !(toosmall | toolarge)	
  	      cexfac <- cexfac[justright]  
  	      plot(x$coordinates, type="n", xlab="", ylab="", asp=1, ...)	
  	      points(x$coordinates[toolarge,,drop=FALSE], col=cols[1], pch=16, cex=5*cex) 
  	      points(x$coordinates[toolarge,,drop=FALSE], col=grey(1), pch=10, cex=5*cex, lwd=2)
  	      # black contours vor better visibility if points overlapp:
  	      # points(x$coordinates[toolarge,], col=1, pch=1, cex=0.92*5*cex)
  	      points(x$coordinates[justright,,drop=FALSE], col=cols[1], pch=16, cex=cexfac*cex)
  	      # points(x$coordinates[justright,], col=1, pch=1, cex=0.92*cexfac*cex)
  	      points(x$coordinates[toosmall,,drop=FALSE], col=cols[1], pch="+", cex=0.4*cex)  	  
  	      # Second pp:  
  	      cexfac <- sqrt(y$mass/massmean)
  	      toosmall <- (cexfac < 0.25)
  	      toolarge <- (cexfac > 5)
  	      justright <- !(toosmall | toolarge)	
  	      cexfac <- cexfac[justright]  
  	      # plot(y$coordinates, type="n", xlab="", ylab="", asp=1, ...)	
  	      points(y$coordinates[toolarge,,drop=FALSE], col=cols[2], pch=16, cex=5*cex) 
  	      points(y$coordinates[toolarge,,drop=FALSE], col=grey(1), pch=10, cex=5*cex, lwd=2)
  	      # black contours vor better visibility if points overlapp:
  	      # points(y$coordinates[toolarge,], col=1, pch=1, cex=0.92*5*cex)
  	      points(y$coordinates[justright,,drop=FALSE], col=cols[2], pch=16, cex=cexfac*cex)
  	      # points(x$coordinates[justright,], col=1, pch=1, cex=0.92*cexfac*cex)
  	      points(y$coordinates[toosmall,,drop=FALSE], col=cols[2], pch="+", cex=0.4*cex)  	  
        } else {
          points(x$coordinates, col=cols[1], pch=16, cex=cex)
          points(y$coordinates, col=cols[2], pch=16, cex=cex)
        }
        invisible()
      }
    } else {
      allcoords <- rbind(x$coordinates,y$coordinates)  
  	  plot(allcoords, type="n", xlab="", ylab="", asp=1, ...)
  	  if (tmass) {
  	    tmassrat <- sum(tplan$mass>0)*(tplan$mass/x$totmass) 
  	    # is ratio transported mass / mean transported mass (over non-negative transports)
  	    # 1 if same amount of mass on each non-zero transport
  	    # The following solution is probably a bit extreme if point patterns are big and general
  	    # (everything that is not within a factor 10 of mean mass is cut off)
  	    colfac <- 1-log10(tmassrat)
  	    toosmall <- (colfac > 2)  # mass to small i.e. grey value too large
  	    toolarge <- (colfac < 0)    
  	    justright <- !(toosmall | toolarge)	
  	    colfac <- colfac[justright]
  	    # im folgenden kein drop nötig, da tplan ein data.frame ist
  	    tplansmall <- tplan[toosmall,]
  	    tplanlarge <- tplan[toolarge,]
  	    tplanright <- tplan[justright,]
        segments(x$coordinates[tplanlarge$from,1], x$coordinates[tplanlarge$from,2],
          y$coordinates[tplanlarge$to,1], y$coordinates[tplanlarge$to,2], col=1, lwd=1.5*lwd)
        segments(x$coordinates[tplanright$from,1], x$coordinates[tplanright$from,2],
          y$coordinates[tplanright$to,1], y$coordinates[tplanright$to,2], col=grey(colfac*aglevel), lwd=lwd)
        segments(x$coordinates[tplansmall$from,1], x$coordinates[tplansmall$from,2],
          y$coordinates[tplansmall$to,1], y$coordinates[tplansmall$to,2], col=grey(0.8), lwd=lwd, lty=2)
      } else {
      	segments(x$coordinates[tplan$from,1], x$coordinates[tplan$from,2], y$coordinates[tplan$to,1],
      	  y$coordinates[tplan$to,2], col=acol, lwd=lwd)
      }
  	  if (pmass) {
  	    massmean <- mean(c(x$mass,y$mass))
  	    # First pp:
  	    cexfac <- sqrt(x$mass/massmean)
  	    toosmall <- (cexfac < 0.25)
  	    toolarge <- (cexfac > 5)
  	    justright <- !(toosmall | toolarge)	
  	    cexfac <- cexfac[justright]  
  	    # plot(x$coordinates, type="n", xlab="", ylab="", asp=1, ...)	
  	    points(x$coordinates[toolarge,,drop=FALSE], col=cols[1], pch=16, cex=5*cex) 
  	    points(x$coordinates[toolarge,,drop=FALSE], col=grey(1), pch=10, cex=5*cex, lwd=2)
  	    # black contours vor better visibility if points overlapp:
  	    # points(x$coordinates[toolarge,], col=1, pch=1, cex=0.92*5*cex)
  	    points(x$coordinates[justright,,drop=FALSE], col=cols[1], pch=16, cex=cexfac*cex)
  	    # points(x$coordinates[justright,], col=1, pch=1, cex=0.92*cexfac*cex)
  	    points(x$coordinates[toosmall,,drop=FALSE], col=cols[1], pch="+", cex=0.4*cex)  	  
  	    # Second pp:  
  	    cexfac <- sqrt(y$mass/massmean)
  	    toosmall <- (cexfac < 0.25)
  	    toolarge <- (cexfac > 5)
  	    justright <- !(toosmall | toolarge)	
  	    cexfac <- cexfac[justright]  
  	    # plot(y$coordinates, type="n", xlab="", ylab="", asp=1, ...)	
  	    points(y$coordinates[toolarge,,drop=FALSE], col=cols[2], pch=16, cex=5*cex) 
  	    points(y$coordinates[toolarge,,drop=FALSE], col=grey(1), pch=10, cex=5*cex, lwd=2)
  	    # black contours vor better visibility if points overlapp:
  	    # points(x$coordinates[toolarge,], col=1, pch=1, cex=0.92*5*cex)
  	    points(y$coordinates[justright,,drop=FALSE], col=cols[2], pch=16, cex=cexfac*cex)
  	    # points(y$coordinates[justright,], col=1, pch=1, cex=0.92*cexfac*cex)
  	    points(y$coordinates[toosmall,,drop=FALSE], col=cols[2], pch="+", cex=0.4*cex)  	  
      } else {
        points(y$coordinates, col=cols[1], pch=16, cex=cex)
        points(y$coordinates, col=cols[2], pch=16, cex=cex)
      }
      invisible()
    }
  }
  par(oldpars)
}


print.wpp <- function(x, ...) {
  stopifnot(class(x) == "wpp")
  cat("Pattern of ",x$N," points in ",x$dimension, " dimensions with total mass ", x$totmass,".\n", sep="")
  cat("Minimal coordinates:", apply(x$coordinates,2,min), "\n")
  cat("Maximal coordinates:", apply(x$coordinates,2,max), "\n")
  invisible(x)
}


summary.wpp <- function(object, ...) {
  print.wpp(object)
}


# ---------------
# pgrid and wpp
# ---------------
# 
# plot semidiscrete transport maps
# 
plot_pgrid_wpp <- function(x,y,tplan,pmass=TRUE,cex=0.8,length=0.1,acol="#996699",col=4,lwd=1.5,rot=TRUE,...) {
  # korrekterweise müsste bei rot=FALSE wpp gedreht werden
  stopifnot(class(x) == "pgrid" && class(y) == "wpp")
  stopifnot(class(tplan) %in% c("apollonius_diagram", "power_diagram"))
  
  tplansites <- as.matrix(tplan$sites[,1:2])
  if (!isTRUE(all.equal(y$coordinates[order(y$coordinates[,1]),],tplansites[order(tplansites[,1]),],check.attributes = FALSE))) {
  	stop("y and target sites of transport tesselation do not match")
  }
  a <- x
  if (a$dimension != 2) stop("plotting of transport from pgrid to wpp is currently only supported in 2-d")
  xi <- a$generator[[1]]
  eta <- a$generator[[2]]  
  # b <- y
  # N <- b$N
  # if (b$dimension != 2) stop("plotting of transport from pgrid to wpp is currently only supported in 2-d")
  image2(xi, eta, a$mass, rot=rot, col=grey(0:200/200), asp=1, xlab="", ylab="")
  
  if (class(tplan) == "apollonius_diagram") {
    plot_apollonius(y$coordinates, tplan$weights, show_points = TRUE,
                    show_weights = FALSE, add_to_weights = 0, add = TRUE, col=col, lwd=lwd, ...)
  } else {  
    plot(tplan, weights=FALSE, add=TRUE, col=col, lwd=lwd, ...)
  }
  
  if (pmass) {
  	massmean <- mean(y$mass)
  	cexfac <- sqrt(y$mass/massmean)
  	toosmall <- (cexfac < 0.25)
  	toolarge <- (cexfac > 5)
  	justright <- !(toosmall | toolarge)	
  	cexfac <- cexfac[justright]  
  	# drop im folgenden wichtig, sonst werden bei nur einer Zeile Punkte (c_11, 1) und (c_12, 2) geplottet 
  	points(y$coordinates[toolarge,,drop=FALSE], col="#CC0000", pch=16, cex=5*cex) 
  	points(y$coordinates[toolarge,,drop=FALSE], col=grey(1), pch=10, cex=5*cex, lwd=2)
  	# black contours vor better visibility if points overlapp:
  	# points(x$coordinates[toolarge,], col=1, pch=1, cex=0.92*5*cex)
  	points(y$coordinates[justright,,drop=FALSE], col="#CC0000", pch=16, cex=cexfac*cex)
  	# points(x$coordinates[justright,], col=1, pch=1, cex=0.92*cexfac*cex)
  	points(y$coordinates[toosmall,,drop=FALSE], col="#CC0000", pch="+", cex=0.4*cex)  	  
  }
  
  # Draw the arrows (not needed for apollonius diagram):
  if(class(tplan) == "power_diagram") {
    cells <- tplan$cells[!sapply(tplan$cells,function(cc){all(is.na(cc))})]  # remove na cells
    # centroid formula from https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
    # area <- 0.5 * sum(x[-n]*y[-1]-x[-1]*y[-n])
    # 1-d: centroidx <- sum((x[-n]+x[-1])*(x[-n]*y[-1]-x[-1]*y[-n]))/(6*area)
    # 1-d: centroidy <- sum((y[-n]+y[-1])*(x[-n]*y[-1]-x[-1]*y[-n]))/(6*area)
    centroid <- function(xx,yy) {
  	  xx <- c(xx,xx[1])
  	  yy <- c(yy,yy[1])
  	  n <- length(xx)
  	  area <- 0.5 * sum(xx[-n]*yy[-1]-xx[-1]*yy[-n])
  	  cx <- sum((xx[-n]+xx[-1])*(xx[-n]*yy[-1]-xx[-1]*yy[-n]))/(6*area)
  	  cy <- sum((yy[-n]+yy[-1])*(xx[-n]*yy[-1]-xx[-1]*yy[-n]))/(6*area)
  	  return(c(cx,cy))
    }
    centroids <- sapply(cells,function(cc) {centroid(cc[,1],cc[,2])})  # compute centroids (gives 2 x {no. of cells} matrix)
    arrows(centroids[1,],centroids[2,],tplan$sites[,1],tplan$sites[,2],lwd=lwd*1,col=acol,angle=20,length=length)
  } 
  
  return(invisible())
}



# ---------------
# Minor methods for classes pgrid and pp
# ---------------
# 

# new method for R-generic all.equal 
all.equal.pgrid <- function(target, current, ...) {
  class(target) <- "list"
  class(current) <- "list"
  NextMethod("all.equal")
}

all.equal.pp <- function(target, current, ...) {
  class(target) <- "list"
  class(current) <- "list"
  NextMethod("all.equal")
}

all.equal.wpp <- function(target, current, ...) {
  class(target) <- "list"
  class(current) <- "list"
  NextMethod("all.equal")
}



# the following is very minimalistic
# --> include inner consistency tests of the two objects in later versions
compatible <- function(target, current, ...) {
  stopifnot(class(target) == class(current))	
  UseMethod("compatible")
}

compatible.pgrid <- function(target, current, ...) {
  return(all(target$n == current$n) && isTRUE(all.equal(target$generator, current$generator)))
}

compatible.pp <- function(target, current, ...) {
  return(target$N == current$N && target$dimension == current$dimension)
}

compatible.wpp <- function(target, current, ...) {
  return(target$dimension == current$dimension && isTRUE(all.equal(target$totmass, current$totmass)))
}



# transforms integer measure-vectors to integer measure-vectors of fixed total mass
# N is the target sum of the measure, afaics it's not per se a problem to go beyond
# .Machine$integer.max if it's not due to individual entries
fudge <- function(temp,N=1e9) {
  n <- length(temp)
  if (sum(temp) < N) {
     sel <- sample(1:n,N-sum(temp),replace=TRUE)
     for (i in sel) {
   	   temp[i] <- temp[i]+1
     }
  } else while (sum(temp) > N) {
     sel <- sample(which(temp>1),sum(temp)-N,replace=TRUE)
     temp[sel] <- temp[sel]-1
  }
  return(temp)
}

# Generate Positions of a uniform 2-d grid on [0,1]^2
grid_positions<-function(n,m){
  G1<-expand.grid(n:1,1:m)
  G2<-cbind(G1[,2],G1[,1])
  G2<-G2-0.5
  G2[,1]<-G2[,1]/m
  G2[,2]<-G2[,2]/n
  return(G2)
}