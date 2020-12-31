cc.plot <- function(tm,mark,scale,type,a) {
    l <- length(a)
    m <- length(a[[1]])
    pan <- rep(names(a),rep(m,l))
    o <- rep(1:l,rep(m,l))
    x <- rep(tm,l)
    y <- unlist(a)
    g <- rep("stat",length(y))
    for (i in 1:l) {
        cl <- attr(a[[i]],"CL")
        if (!is.null(cl)) {
            for (j in 1:3) {
                    pan <- c(pan,names(a)[i],names(a)[i])
                    o <- c(o,i,i)
                    x <- c(x,min(tm)-1,max(tm)+1)
                    y <- c(y,cl[j],cl[j])
                    g <- c(g,rep(paste("vcl",j),2))
            }
        }
    }
    pan <- reorder(pan,-o)
    j <- l+1
    ggpanel <- function(x,y,...) {
        lattice::panel.xyplot(x,y,...)
        if (m<=100) lattice::panel.points(x[1:m],y[1:m],pch=20)
        j <<- j-1
        cl <- attr(a[[j]],"CL")
        if (!is.null(cl)) {
            if (mark=="all") {
                out <- integer(0)
                if (!is.na(cl[2])) out <- which(y[1:m]<cl[2])
                if (!is.na(cl[3])) out <- c(out,which(y[1:m]>cl[3]))
                for (i in out) lattice::panel.text(x[i],y[i],x[i],cex=1.2)
            }
            if (mark=="one") {
                out <- NA
                if (!is.na(cl[2])) {
                    i <- which.min(y[1:m])
                    if (y[i]<cl[2]) out <- i
                }
                if (!is.na(cl[3])) {
                    i <- which.max(y[1:m])
                    if ((y[i]>cl[3]) &
                        (is.na(out) || ((y[i]-cl[3])>=(cl[2]-y[out]))))
                        out <- i
                }
                if (!is.na(out)) lattice::panel.text(x[out],y[out],x[out],cex=1.2)
            }
        }
        ff <- attr(a[[j]],"fit")
        if (!is.null(ff)) {
            lattice::llines(c(x-0.5,x[m]+0.5),c(ff,ff[m]),lty="dashed",
                   type="s",col=lattice::trellis.par.get()$superpose.line$col[2])
            i <- 2
            if (mark=="all") {
                while (i<=m) {
                    if (ff[i]!=ff[i-1]){
                        if (i==2) {
                            lattice::panel.text(x[1],ff[1],x[1]); 
                        } else if ((i==m)||((i<m)&&(ff[i+1]==ff[i-1]))) {
                            lattice::panel.text(x[i],ff[i],x[i]); i<-i+1;
                        } else {
                            lattice::panel.text(x[i]-0.5,ff[i],x[i])
                        }
                    }
                    i <- i+1
                }
            }
        }
    }
    col <- c(lattice::trellis.par.get()$superpose.line$col[1],
             rep(lattice::trellis.par.get()$superpose.line$col[2],3))
    print(lattice::xyplot(y~x|pan,type=type,distribute.type=TRUE,
                          xlab="Time",ylab="",lty=c(1,2,3,3),col=col,
                          group=g,scales=list(y=scale),layout=c(1,l),panel=ggpanel))
}


dmodq <- function(stat,alpha) {
    eps <- .Machine$double.eps
    if (alpha<=eps) {
        q <- max(stat)+eps
    } else if (alpha>=1-eps) {
        q <- min(stat)-eps
    } else {
        q <- quantile(stat,1-alpha)
    }
    q
}

dbalance2 <- function(stat,alpha) {
    obj <- function(beta) {
        q <- apply(stat,1,dmodq,beta)
        mean((stat[1,]>q[1])|(stat[2,]>q[2]))-alpha
    }
    apply(stat,1,dmodq,uniroot(obj,c(0,1))$root)
}

dbalance3 <- function(stat,alpha) {
    obj <- function(beta) {
        q <- apply(stat[2:3,],1,dmodq,beta)
        l <- c(dmodq(stat[1,],mean((stat[2,]>q[1])|(stat[3,]>q[2]))),q)
        mean((stat[1,]>l[1])|(stat[2,]>l[2])|(stat[3,]>l[3]))-alpha
    }
    beta <- uniroot(obj,c(0,1))$root
    q <- apply(stat[2:3,],1,dmodq,beta)
    c(dmodq(stat[1,],mean((stat[2,]>q[1])|(stat[3,]>q[2]))),q)
}



check <- function(a) {
    check1 <- function(x) {
        cl <- attr(x,"CL")
        (!is.na(cl[2]) && any(x<cl[2]))||(!is.na(cl[3]) && any(x>cl[3]))
    }
    any(sapply(a$value,check1))
}
