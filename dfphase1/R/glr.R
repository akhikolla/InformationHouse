changepoint <- function(x, subset, score=c("Identity","Ranks"), only.mean=FALSE, plot=TRUE,
                        FAP=0.05, seed=11642257, L=10000, limits=NA) {
    d <- dim(x)
    if (is.null(d)) {
        dim(x) <- d <- c(1,1,length(x))
    } else {
        if (length(d)!=2) stop("x must be a nxm matrix")
        dim(x) <- d <- c(1,d)
    }
    score <- match.arg(score)
    if (score=="Ranks") score <- "Marginal Ranks"
    u <- mchangepoint(x,subset,score,only.mean,plot,FAP,seed,L,limits)
    u$score <- score
    invisible(u)
}


mchangepoint <- function(x, subset, score=c("Identity","Signed Ranks","Spatial Signs",
                                            "Spatial Ranks","Marginal Ranks"),
                         only.mean=FALSE, plot=TRUE,
                         FAP=0.05, seed=11642257, L=10000, limits=NA) {
    d <- dim(x)
    if ((length(d)<2) || (length(d)>3)) stop("x must be a pxnxm array")
    if (length(d)==2) dim(x) <- d <- c(d[1],1,d[2])
    tm <- 1:d[3]
    if (!missing(subset)) {
        x <- x[,,subset,drop=FALSE]
        d <- dim(x)
        tm <- tm[subset]
    }
    if (length(tm)<5)
        stop("Number of subgroups/times must be greater than 5")
    score <- match.arg(score)
    if (score %in% c("Signed Ranks","Spatial Signs","Spatial Ranks")) {
        x <- ggrscore(x,score)$score
    } else if (score=="Marginal Ranks") {
        for (i in 1:d[1]) x[i,,] <- rank(x[i,,])
    }
    if (is.na(limits)) {
        if (L<1000) stop("The number of permutation is too low")
        if (!is.na(seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                kept <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                on.exit(assign(".Random.seed", kept, envir = .GlobalEnv))
            }
            set.seed(seed)
        }
    } else {
        if (length(limits)!=dim(x)[3]) stop("limits's length is not correct")
        L <- 0
        seed <- NA
    }
    a <- ggglrchart(x,only.mean,L)
    if (is.na(limits)) {
        u <- rowMeans(a$glr.perm)
        a$limits <- u*quantile(apply(a$glr.perm/u,2,max,na.rm=TRUE),1-FAP)
        a$glr.perm <- NULL
    } else {
        a$limits <- limits
    }
    if (!only.mean) {
        u <- ggglrchart(x,TRUE,0)
        d <- dim(x)
        N <- d[2]*d[3]
        a$mean <- N*log(1+u$glr/N)
        a$dispersion <- a$glr-a$mean
    }
    if (plot) {
        if (only.mean) {
            v <- list("Adjusted Likelihood Ratio"=a$glr/a$limits)
            attr(v[[1]],"CL") <- c(NA,0,1)
        } else {    
            v <- list("Adjusted Likelihood Ratio"=a$glr/a$limits,
                      "Diagnostic / Level"=a$mean/a$limits,
                      "Diagnostic / Scale"=a$dispersion/a$limits) 
            attr(v[[1]],"CL") <- c(NA,0,1)
            attr(v[[2]],"CL") <- attr(v[[3]],"CL") <- c(NA,0,NA)
        }
        cc.plot(tm,"one","same","l",v)
    }
    a$score <- score
    a$seed <- seed
    a$L <- L
    a$only.mean <- only.mean
    for (i in c("glr","limits","mean","dispersion"))
        if (!is.null(a[[i]])) names(a[[i]]) <- tm
    invisible(a)
}

changepoint.normal.limits <- function(n, m, score=c("Identity","Ranks"), only.mean=FALSE, 
                                      FAP=0.05, seed=11642257, L=100000) {
    score <- match.arg(score)
    if (score=="Ranks") score <- "Marginal Ranks"
    mchangepoint.normal.limits(1,n,m,score,only.mean,FAP,seed,L)
}

mchangepoint.normal.limits <- function(p, n, m,
                                       score=c("Identity","Signed Ranks","Spatial Signs",
                                               "Spatial Ranks","Marginal Ranks"),
                                       only.mean=FALSE, FAP=0.05, seed=11642257, L=100000) {
    limits <- rep(1,m)
    v <- replicate(L,
                   mchangepoint(array(rnorm(p*n*m),c(p,n,m)),
                                score=score,only.mean=only.mean,
                                plot=FALSE,FAP=FAP,seed=NA,L=0,limits=limits)$glr)
    u <- rowMeans(v)
    u*quantile(apply(v/u,2,max,na.rm=TRUE),1-FAP)    
}



