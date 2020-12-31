oo <- function(p,layout) {
    g <- layout[1]
    blk <- g*layout[2]
    oo <- NULL
    start <- 1
    end <- blk
    while (start<=p) {
        oo <- c(oo,end:start)
        start <- end+1
        end <- min(end+blk,p)  
    }
    if (g>1) {
        i1 <- 1
        i2 <- g
        while (i1<=p) {
            oo[i1:i2] <- rev(oo[i1:i2])
            i1 <- i2+1
            i2 <- min(p,i2+g)
        }
    }
    oo
}

phase1Plot <- function(x) {
    d <- dim(x)
    if (is.null(d)) {
        dim(x) <- d <- c(1,1,length(x))
    } else {
        if (length(d)!=2) stop("x must be a nxm matrix")
        dim(x) <- d <- c(1,d)
    }
    mphase1Plot(x)
}

mphase1Plot <- function(x,layout=c(1,p)) {
    d <- dim(x)
    if ((length(d)<2) || (length(d)>3)) stop("x must be a pxnxm array")
    if (length(d)==2) dim(x) <- d <- c(d[1],1,d[2])
    p <- d[1]
    n <- d[2]
    m <- d[3]
    nm <- dimnames(x)[[1]]
    if (is.null(nm)) nm <- paste("X",1:p,sep="")
    a <- expand.grid(nm,1:n,1:m)
    xf <- 1:m
    print(lattice::xyplot(
        x~a[,3]|reorder(a[,1],rep(oo(p,layout),m*n)),
        layout=layout, scales=list(y="free"),
        xlab="Time",ylab="Observations",
        panel = function(x,y,...) {
            lattice::panel.xyplot(x,y,...)
            f <- tapply(y,x,mean)
            lattice::llines(xf,f,lty="solid",col=trellis.par.get()$superpose.line$col[2])
         }
        ))
}

