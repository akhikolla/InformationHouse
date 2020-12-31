rsp <- function(y, plot=TRUE, L=1000, seed=11642257, alpha=0.05,
                maxsteps=min(50,round(NROW(y)/15)),
                lmin=max(5,min(10,round(NROW(y)/10)))) {
    d <- dim(y)
    if (is.null(d)) {
        dim(y) <- c(length(y),1)
    } else {
        if (length(d)!=2) stop("y must be a nxm matrix")
        y <- t(y)        
    }
    opt <- list(m=NROW(y),n=NCOL(y),
                L=L,seed=seed,maxsteps=maxsteps,lmin=lmin,alpha=alpha)
    if (!is.na(opt$seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
        }
        set.seed(opt$seed)
    }
    m <- .rsp(y,"level",opt)
    s <- .rsp(y-m$fit,"scale",opt)
    p <- p.adjust(c(m$p,s$p))
    if (p[1]>=alpha) m$fit <- rep(mean(m$yi),length(m$yi))
    if (p[2]>=alpha) s$fit <- rep(mean(s$yi),length(s$yi))
    if (plot) {
        a <- list(m$yi,s$yi)
        attr(a[[1]],"fit") <- m$fit
        attr(a[[2]],"fit") <- s$fit
        pp <- function(p) if (p<0.001) "(p<0.001)" else paste("(p=",round(p,3),")",sep="")  
        names(a) <- c(paste("Level",pp(p[1])),
                      paste("Scale",pp(p[2])))
        cc.plot(1:opt$m,"all","free","l",a)
    }
    ans <- list(p=p,stat=cbind(m$yi,s$yi),fit=cbind(m$fit,s$fit))
    names(ans$p) <- colnames(ans$stat) <- colnames(ans$fit) <- c("level","scale")
    if (plot) invisible(ans) else ans
}

.rsp <- function(y, type, opt) {
    nsteps <- opt$maxsteps
    ipar <- c(opt$m,opt$n,nsteps,opt$lmin,if (type=="level") 1 else 2,opt$L)
    nstat <- if (opt$n==1) nsteps else nsteps+1
    mod <- .C("ggdotrsp",as.integer(ipar),as.double(y),
              steps=integer(1+2*(nsteps+1)),stat=double(nstat),
              perm=double(nstat*opt$L),PACKAGE="dfphase1")
    perm <- matrix(mod$perm,nstat)
    a <- apply(perm,1,mean)
    b <- pmax(apply(perm,1,sd),.Machine$double.eps)
    r <- (cbind(mod$stat,perm)-a)/b
    w <- apply(r,2,max)
    p <- (sum(w>=w[1])-1)/opt$L
    yi <- if (type=="level") rowMeans(y) else sqrt(rowMeans(y*y))
    if (p>opt$alpha) {
        fit <- rep(mean(yi),length(yi))
    } else {
        j <- which.max(r[,1])
        nstat <- NROW(r)
        if ((NCOL(y)==1) || (j<nstat)) {
            g <- factor(.C("ggstepfactor",as.integer(j+1),
                           as.integer(mod$steps),integer(NROW(y)),
                           PACKAGE="dfphase1")[[3]])
            fit <- tapply(yi,g,mean)[g]
        } else {
            q <- quantile(w,1-opt$alpha)
            out <- if (type=="level") yi else yi*yi
            out <- (abs(out-mean(out))-a[nstat])/b[nstat]
            fit <- ifelse(out>=q,yi,mean(yi[out<q]))
        }
    }
    list(p=p,yi=as.double(yi),fit=as.double(fit))
}

