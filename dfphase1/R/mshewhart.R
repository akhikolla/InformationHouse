mshewhart <- function(x, subset, stat=c("T2Var","T2","Var","Depth Ranks"),
                      score=c("Identity","Signed Ranks","Spatial Signs",
                              "Spatial Ranks","Marginal Ranks"),
                      loc.scatter = c("Classic","MCD"),
                      plot=TRUE, FAP=0.05, seed=11642257, L=1000, limits=NA) {
    d <- dim(x)
    if (length(d)!=3) stop("x must be a pxnxm array")
    if (d[2]<2) stop("Size of each subgroup must be greater than 1")
    tm <- 1:d[3]
    if (!missing(subset)) {
        x <- x[,,subset,drop=FALSE]
        d <- dim(x)
        tm <- tm[subset]
    }
    if (length(tm)<1)
        stop("Number of subgroups/times must be greater than 1 (after subsetting)")
    stat <- match.arg(stat)
    score <- match.arg(score)
    if (is.na(limits)) {
        if (L<100) stop("The number of permutation is too low")
        if (!is.na(seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
            }
            set.seed(seed)
        }
    } else {
        if (switch(stat,T2Var=2,1)!=length(limits))
            stop("Uncorrected number of control limits")
        L <- 0
        seed <- NA
    }
    if (stat=="Depth Ranks") {
        v <- ggdepthranks(x,L)
        score <- "Depth Ranks"
    } else if (score=="Identity") {
        loc.scatter <- match.arg(loc.scatter)
        if (loc.scatter=="MCD") {
            dim(x) <- c(d[1],d[2]*d[3])
            u <- robustbase::covMcd(t(x))
            sc <- chol(u$cov)
            x <- solve(t(sc),x-u$center)
            dim(x) <- d
            v <- ggscore2mshewhart(x,stat,L)
            v$center <- u$center
            v$scatter <- sc
        } else {
            v <- ggclassicmshewhart(x,stat,L)
        }
    } else  if (score %in% c("Signed Ranks","Spatial Signs","Spatial Ranks")) {
        u <- ggrscore(x,score)
        v <- ggscore2mshewhart(u$score,stat,L)
        v$center <- u$center
        v$scatter <- u$scatter
    } else if (score=="Marginal Ranks") {
        N <- d[2]*d[3]
        dim(x) <- c(d[1],N)
        for (i in 1:d[1]) x[i,] <- rank(x[i,])
        x <- x-(N+1)/2
        sc <- chol(tcrossprod(x/sqrt(N)))
        x <- solve(t(sc),x)
        dim(x) <- d
        v <- ggscore2mshewhart(x,stat,L)
        v$center <- rep(N/2,d[1])
        v$scatter <- sc
    }
    if (is.na(limits)) {
        if (stat=="T2Var") {
            v$limits <- dbalance2(v$sp,FAP)
        } else {
            v$limits <- dmodq(v$sp,FAP)
        }
    } else {
        v$limits <- limits
    }
    v$stat <- stat
    v$sp <- NULL
    v$score <- score
    v$L <- L
    v$seed <- seed
    v$FAP <- FAP
    if (plot) {
        if (stat=="T2Var") {
            u <- list(T2=v$statistic[1,],Var=v$statistic[2,])
            attr(u$T2,"CL") <- c(0,NA,v$limits[1])
            attr(u$Var,"CL") <- c(0,NA,v$limits[2])
            cc.plot(tm,"all","free",c("h","l","l","l"),u)
        } else {
            u <- list(v$statistic)
            names(u) <- stat
            attr(u[[1]],"CL") <- c(0,NA,v$limits)
            cc.plot(tm,"all","same",c("h","l","l","l"),u)
        }
    }
    if (stat == "T2Var") {
        v$T2 <- v$statistic[1,]
        v$Var <- v$statistic[2,]
    } else if (stat == "T2") {
        v$T2 <- v$statistic
    } else if (stat == "Var") {
        v$Var <- v$statistic
    } else {
        v$DepthRanks <- v$statistic
    }
    for (i in c("T2","Var","DepthRanks"))
        if (!is.null(v[[i]])) names(v[[i]]) <- tm
    v$statistic <- NULL
    invisible(v)
}

mshewhart.normal.limits <- function(p, n, m, stat=c("T2Var","T2","Var", "Depth Ranks"),
                                    score=c("Identity","Signed Ranks","Spatial Signs",
                                            "Spatial Ranks","Marginal Ranks"),
                                    loc.scatter = c("Classic","MCD"),
                                    FAP=0.05, seed=11642257, L=100000) {
    stat <- match.arg(stat)
    score <- match.arg(score)
    loc.scatter <- match.arg(loc.scatter)
    if (!is.na(seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
        }
        set.seed(seed)
    }
    if (stat=="T2Var") {
        l <- c(0,0)
        v <- replicate(L,{
            x <- array(rnorm(p*n*m),c(p,n,m))
            u <- mshewhart(x,stat=stat,score=score,loc.scatter=loc.scatter,plot=FALSE,limits=l)
            apply(u$statistic,1,max)
        })
        dbalance2(v,FAP)
    } else {
        v <- replicate(L,{
            x <- array(rnorm(p*n*m),c(p,n,m))
            u <- mshewhart(x,stat=stat,score=score,loc.scatter=loc.scatter,plot=FALSE,limits=0)
            max(if (stat=="T2") u$T2 else if (stat=="Var") u$Var else u$DepthRank)
        })
        dmodq(v,FAP)
    }
}
