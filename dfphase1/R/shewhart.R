shewhart <- function(x, subset, stat=c("XbarS","Xbar","S","Rank","lRank","sRank"),
                     aggregation=c("median","mean"), plot=TRUE,
                     FAP=0.05, seed=11642257,
                     L = if(stat %in% c("XbarS","Xbar","S")) 1000 else 100000,
                     limits=NA) {
    if (!is.matrix(x) || !is.numeric(x))
        stop("x must be a numeric nxm matrix")
    tm <- 1:NCOL(x)
    if (!missing(subset)) {
        x <- x[,subset,drop=FALSE]
        tm <- tm[subset]
    }
    if ((NCOL(x)<2) || (NROW(x)<2))
        stop("n and m must be greater than 1 (after subsetting)")
    stat <- match.arg(stat)
    if (is.na(limits)) {
        if (L<100) stop("L is too low") 
        if (!is.na(seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
            }
            set.seed(seed)
        }
    } else {
        if (switch(stat,XbarS=3,Rank=2,S=2,1)!=length(limits))
            stop("Uncorrected number of control limits")
        L <- 0
        seed <- NA
    }
    if (stat %in% c("XbarS","Xbar","S")) {
        aggregation <- match.arg(aggregation)
        u <- ggxbars(x,aggregation=="mean",L)
        if (stat=="Xbar") u$S <- NULL
        if (stat=="S") u$Xbar <- NULL
        if (is.na(limits)) limits <- .shewhart.limits(u$perm,stat,FAP)
        u$perm <- NULL
    } else {
        u <- ggrank(x)
        if (stat=="lRank") u$sRank <- NULL
        if (stat=="sRank") u$lRank <- NULL
        if (is.na(limits))
            limits <- shewhart.normal.limits(NROW(x),NCOL(x),stat,FAP=FAP,seed=seed,L=L)
    }
    u$limits <- limits
    u$stat <- stat
    u$L <- L
    u$aggregation <- aggregation
    u$FAP <- FAP
    u$seed <- seed
    if (plot) {
        if (stat=="XbarS") {
            v <- list(Xbar=u$Xbar,S=u$S)
            dl <- u$limits[1]*u$scale/sqrt(NROW(x))
            attr(v[["Xbar"]],"CL") <- c(u$center,u$center-dl,u$center+dl)
            attr(v[["S"]],"CL") <- c(u$scale,u$scale*u$limits[2],u$scale*u$limits[3]) 
        } else if (stat=="Xbar") {
            v <- list(Xbar=u$Xbar)
            dl <- u$limits*u$scale/sqrt(NROW(x))
            attr(v[["Xbar"]],"CL") <- c(u$center,u$center-dl,u$center+dl)
        } else if (stat=="S") {
            v <- list(S=u$S)
            attr(v[["S"]],"CL") <- c(u$scale,u$scale*u$limits[1],u$scale*u$limits[2]) 
        } else if (stat=="Rank") {
            v <- list(lRank=u$lRank,sRank=u$sRank)
            attr(v[["lRank"]],"CL") <- c(0,-u$limits[1],u$limits[1])
            attr(v[["sRank"]],"CL") <- c(0,-u$limits[2],u$limits[2])
        } else if (stat=="lRank") {
            v <- list(lRank=u$lRank)
            attr(v[["lRank"]],"CL") <- c(0,-u$limits,u$limits)
        } else if (stat=="sRank") {
            v <- list(sRank=u$sRank)
            attr(v[["sRank"]],"CL") <- c(0,-u$limits,u$limits)
        }
        cc.plot(tm,"all","free","l",v)
    }
    for (i in c("Xbar","S","lRank","sRank"))
        if (!is.null(u[[i]])) names(u[[i]]) <- tm
    invisible(u)
}

.shewhart.limits <- function(v,stat,FAP) {
    switch(stat,
           XbarS={l<-dbalance3(v,FAP);l[2]<- (-l[2]);l},
           Xbar=dmodq(v[1,],FAP),
           S={l<-dbalance2(v[2:3,],FAP);l[1] <- -l[1]; l},
           Rank=dbalance2(v,FAP),
           lRank=dmodq(v[1,],FAP),
           sRank=dmodq(v[2,],FAP)
           )
}

shewhart.normal.limits <- function(n, m, stat=c("XbarS","Xbar","S","Rank","lRank","sRank"),
                                   aggregation=c("median","mean"),
                                   FAP=0.05, seed=11642257, L=100000) {
    if ((n<2) || (m<2)) stop("n and m must be greater than 1")
    stat <- match.arg(stat)
    aggregation <- match.arg(aggregation)
    if (!is.na(seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
        }
        set.seed(seed)
    }
    .shewhart.limits(if (stat %in% c("XbarS","Xbar","S"))
                         ggxbarsall(n,m,aggregation=="mean",L)
                     else
                        ggrankall(n,m,L), 
                     stat,FAP)
}
