mphase1 <- function(x, plot=TRUE, post.signal=TRUE,
                    isolated=dim(x)[2]>1, step=TRUE,
                    alpha=0.05, gamma=0.5,
                    K=min(50,round(sqrt(dim(x)[3]))), lmin=5, L=1000,
                    seed=11642257) {
    d <- dim(x)
    if ((length(d)<2) || (length(d)>3)) stop("x must be a pxnxm array")
    if (length(d)==2) dim(x) <- c(d[1],1,d[2])
    if (L<100) stop("The number of permutation L is too small")
    if (K<1) stop("The number of potential shifts K is too small")
    if (lmin<1) stop("The minimum length of a step shift lmin is too small")
    if (isolated && (dim(x)[2]<=1)) stop("Detection of isolated shifts requires subgrouped data")
    if (!isolated && !step) stop("No shift detection required")
    if (!is.na(seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
        }
        set.seed(seed)
    }
    u <- MPHASE1(x,isolated,step,K,lmin,L)
    u$x <- x
    if (is.null(dimnames(x))) {
        nm <- paste("X",1:dim(x)[1],sep="")
        dimnames(u$x) <- list(nm,NULL,NULL)
    } else if (is.null(dimnames(x)[[1]])) {
        nm <- paste("X",1:dim(x)[1],sep="")
        dnm <- dimnames(x)
        dnm[[1]] <- nm
        dimnames(u$x) <- dnm
    } else {
        nm <- dimnames(x)[[1]]
    }
    dimnames(u$signed.ranks) <- dimnames(x)
    if (!post.signal) plot <- FALSE
    if (plot) post.signal <- TRUE
    u$call <- match.call()
    u$call$plot <- plot
    u$call$post.signal <- post.signal
    u$call$isolated <- isolated
    u$call$step <- step
    u$call$alpha <- alpha
    u$call$gamma <- gamma
    u$call$K <- K
    u$call$lmin <- 5
    u$call$L <- L
    u$call$seed <- seed
    u$center <- as.numeric(crossprod(u$r,u$center))
    names(u$center) <- nm
    u$scatter <- crossprod(u$r)
    dimnames(u$scatter) <- list(nm,nm)
    class(u) <- "mphase1"
    if (post.signal) postsignal.mphase1(u,plot,alpha,gamma) else u
}

postsignal <- function(x,...) UseMethod("postsignal")


postsignal.mphase1 <- function(x,plot=TRUE,alpha=0.05,gamma=0.5,...) {
    x$call$plot <- plot
    x$call$post.signal <- TRUE
    x$call$alpha <- alpha
    x$call$gamma <- gamma
    dd <- dim(x$x)
    p <- dd[1]
    n <- dd[2]
    m <- dd[3]
    if (x$p.value<alpha) {
        mn <- m*n
        n1 <- n-1
        ncp <- NROW(x$forward)
        X <- matrix(0,mn,ncp+1)
        X[,1] <- 1
        for (i in 1:ncp) {
            a <- (x$forward$time[i]-1)*n+1
            X[if (x$forward$type[i]=="Isolated") a:(a+n1) else a:mn,i+1] <- 1
        }
        Si <- solve(t(x$r))
        idx <- alasso(XX <- kronecker(X,Si), as.numeric(x$signed.ranks),
                      gamma, if (n==1) p*m-p else 2*p*m-p)
        idx[1:p] <- TRUE
        y <- as.numeric(Si%*%matrix(x$x,p))
        x$fitted <- crossprod(x$r,matrix(qr.fitted(qr(XX[,idx]),y),p))
        dim(x$fitted) <- dd
        x$residuals <- x$x-x$fitted
        jj <- seq(p+1,2*p)
        x$alasso <- data.frame(type=character(0),time=integer(0),variables=character(0)) 
        for (i in 1:ncp) {
            a <- which(idx[jj])
            if (length(a)>0) {
                x$alasso <- rbind(x$alasso,
                                  data.frame(type=x$forward$type[i],
                                             time=x$forward$time[i],
                                             variables=paste(a,collapse=",")))
            }
            jj <- jj+p
        }
    } else {
        ff <- apply(x$x,1,mean) 
        x$fitted <- array(ff,dd)
        x$residuals <- x$x-x$fitted
        x$alasso <- "None"
    }
    dimnames(x$fitted) <- dimnames(x$residuals) <- dimnames(x$x) 
    if (plot) plot.mphase1(x)
    x
}

print.mphase1 <- function(x,...) {
    cat("Call:\n")
    print(x$call)
    cat("\np-value",if (x$p.value<0.001) "< 0.001" else paste("=",x$p.value),"\n")
    if (!is.null(x$alasso)) {
        cat("\nLocation Shifts:\n")
        print(x$alasso)
    }
    invisible(x)
}

plot.mphase1 <- function(x,layout=c(1,p),...) {
    if (is.null(x$alasso))
        object <- postsignal.mphase1(x,FALSE,alpha=x$call$alpha,gamma=x$call$gamma)
    dd <- dim(x$x)
    p <- dd[1]
    n <- dd[2]
    m <- dd[3]
    y <- as.numeric(apply(x$x,c(1,3),mean))
    nm <- dimnames(x$x)[[1]]
    a <- expand.grid(nm,1:m)
    j <- jj <- 0
    uu <- oo(p,layout)
    fitted <- x$fitted[uu,1,]
    eps <- sqrt(.Machine$double.eps)
    col <- c(trellis.par.get()$superpose.line$col[1],
             rep(trellis.par.get()$superpose.line$col[2],3))
    print(xyplot(
        y~a[,2]|reorder(a[,1],rep(uu,m)), type="l", col=col,
        scales=list(y="free"),layout=layout, 
        main=if (x$p.value<0.001) "p-value<0.001" else
                                        paste("p-value=",round(x$p.value,3),sep=""),
        ylab=if (n==1) "Observations" else "Subgroup Means", xlab="Time",
        prepanel = function(x,y,...) {
            j <<- j+1
            list(ylim=range(y,fitted[j,]))
        },        
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            if (m<=100) panel.points(x[1:m],y[1:m],pch=20)
            jj <<- jj+1
            ff <- fitted[jj,]
            llines(c(x-0.5,x[m]+0.5),c(ff,ff[m]),lty="dashed",
                   type="s",col=trellis.par.get()$superpose.line$col[2])
            i <- 2
            while (i<=m) {
                if (abs(ff[i]-ff[i-1])>eps) {
                    if (i==2) {
                        panel.text(1,ff[1],1); 
                    } else if ((i==m)||((i<m)&&(abs(ff[i+1]-ff[i-1])<eps))) {
                        panel.text(i,ff[i],i); i<-i+1;
                    } else {
                        panel.text(i-0.5,ff[i],i)
                    }
                }
                i <- i+1
            }
        },
    ))
    invisible(x)
}


