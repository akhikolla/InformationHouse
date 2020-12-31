tri.mesh <- function(x,y=NULL,duplicate="error"){
    if(is.null(x))
        stop("argument x missing.")
    if(is.null(y)){
        x1<-x$x
        y1<-x$y
        if (is.null(x1) || is.null(y1))
            stop("argument y missing and x contains no $x or $y component.")
    } else {
        x1<-x
        y1<-y
    }

    n <- length(x1)
    if(length(y1)!=n)
        stop("length of x and y differ.")
    ## handle duplicate points:
    xy <- paste(x1, y1, sep =",")
    i <- match(xy, xy)
    if(duplicate!="error")
        {
            if(duplicate!="remove" & duplicate!="error" & duplicate!="strip"){
                stop("possible values for \'duplicate\' are \"error\", \"strip\" and \"remove\"")
            }
            else{
                if(duplicate=="remove")
                    ord <- !duplicated(xy)
                if(duplicate=="strip")
                    ord <- (hist(i,plot=FALSE,freq=TRUE,breaks=seq(0.5,max(i)+0.5,1))$counts==1)
                x1 <- x1[ord]
                y1 <- y1[ord]
                n <- length(x1)
            }
        }
    else
        if(any(duplicated(xy)))
            stop("duplicate data points")

    ans <- shull.deltri(x1,y1)
    nt <- length(ans$i1)

    ## note: triangles are enumerated in c++ starting with 0, so add 1 here
    ## points are enumerated started with 1
    tri.obj<-list(n=ans$n,x=ans$x,y=ans$y,
                  nt=ans$nt,
                  trlist=ans$trlist,
                  cclist=ans$cclist,
                  nchull=ans$nch,
                  chull=ans$ch,
                  narcs=ans$na,
                  arcs=cbind(ans$a1,ans$a2),
                  call=match.call())

    colnames(tri.obj$trlist) <- c("i1","i2","i3","j1","j2","j3","k1","k2","k3")
    colnames(tri.obj$cclist) <- c("x","y","r","area","ratio")
    colnames(tri.obj$arcs) <- c("from","to")
    class(tri.obj)<-"triSht"
    invisible(tri.obj)
}


