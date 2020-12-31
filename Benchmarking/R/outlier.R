# $Id: outlier.R 226 2020-06-27 21:48:14Z lao $
# Undersoeg om der er outliere i datasaet.
# En outlier er en eller flere enheder, der adskiller sig kraftigt 
# fra alle andre enheder. 
# Her defineres outlier som enheder som aendrer arealet/rumfanget af 
# den maengde som alle enheder udspaender. Dvs. det undersoeges om 
# arealet/rumfanget aendres markant hvis en eller flere enheder fjernes.
# Areal/rumfang beregnes som determinant af t(XY) %*% XY, en
# (m+n)x(n+m) matrix; af hensyn til hastighed bruges crossprod(XY). 

# C:\prg\R\bin\x64\Rscript.exe --vanilla test_ap.R


outlier.ap <- function(X, Y, NDEL = 3L, NLEN = 25L, TRANSPOSE = FALSE)  
{
    DEBUG=FALSE
    
    if (TRANSPOSE) {
        X <- t(X)
        Y <- t(Y)
    }
    X <- tjek_data(X)
    Y <- tjek_data(Y)
    if (dim(X)[1] != dim(Y)[1]) stop("Number of firms differ in 'X' and 'Y'")

    xy <- cbind(X,Y)
    XY <- crossprod(xy)
    S <- det(XY)  # det(t(xy)%*%xy)
    R <- NDEL        # hoejeste antal udeladte units
    K <- dim(X)[1]   # number of firms/observations
    last = min(NLEN, K)   # gem 25 mindste vaerdier af RX
    ratio <- array(Inf, c(last,R) )
    imat <- matrix(NA, nrow=R, ncol=R)
    rmin <- array(Inf, dim=R)
    if (DEBUG) { lbnr <- 1 }
    for ( r in 1L:R )  {
        if (DEBUG) print( paste("Remove ",r," observations",sep=""), quote=FALSE )
        # remove r observationer
        RX <- rep(Inf,last)
        rrlast <- Inf
   
        if (K < r) 
            stop("Number of units less than number to be removed; K < NDEL")
        e <- 0
        h <- r     # antal udeladte units
        del <- 1L:r  # Maengde af udeladte; foerste maengde er de 'r' foerste units
        if (DEBUG) {
            cat(format(lbnr, width=5), "*: ", sep=""); lbnr <- lbnr + 1
            print(del)
        }
        Dxy = XY - crossprod(xy[del,,drop=FALSE])
        rr <- det(Dxy)/S  ## det( t(Dxy) %*% Dxy )/S
        RX[1] <- rr
        rrlast <- RX[last]
   
        # count <- as.integer(round(choose(K, r)))
        # print( paste("Number of combinations", count), quote=FALSE )
        nmmp1 <- K - r + 1L   # Indeks for unit een efter sidst mulige medtagne unit 
        while (del[1L] != nmmp1) {
            if (e < K - h) {
                # der er flere der mangler at blive udeladt
                h <- 1L
                e <- del[r]
                j <- 1L
            }
            else {
                if (DEBUG) print("Skift at base")
                e <- del[r - h]
                h <- h + 1L
                j <- 1L:h
            }
            del[r - h + j] <- e + j
            if (DEBUG) {
                cat(format(lbnr, width=5), " : ", sep=""); lbnr <- lbnr + 1
                print(del)
            }
            Dxy = XY - crossprod(xy[del,,drop=FALSE])
            rr <- det(Dxy)/S   ## det( t(Dxy) %*% Dxy )/S
            if ( rr < rrlast ) {  # gem de |last| mindste
                if ( rr < min(RX) ) {
                    imat[r,1:r] <- del
                }
                RX[last] <- rr
                RX <- sort(RX)
                rrlast <- RX[last]
            }
        } ## while
        rmin[r] <- min(RX)
        Rratio <- log(RX/min(RX))
        ratio[,r] <- Rratio
    } ## for (r) 
    return( list(ratio=ratio, imat=imat, r0=rmin) )
}  # outlier.ap



outlier.ap.plot <- function(ratio, NLEN = 25L, xlab="Number of firms deleted", 
                            ylab="Log ratio", ..., ylim)  
{
    nlen <- min(NLEN, dim(ratio)[1])
    R <- dim(ratio)[2]
    if (missing(ylim))  ylim <- range(ratio[1:nlen,])
    ry <- matrix(1:R, nrow=nlen, ncol=R, byrow=TRUE)
    plot(1:R, rep(0,R), type="p", ylim=ylim, xaxt="n", ylab=ylab, xlab=xlab, ...)
    axis(1, at=1:R, labels=c(1:R))
    points(ry[-1,], ratio[2:nlen,])
    lines(1:R, ratio[2,], lty="dashed")
}  # outlier.ap.plot



#----------------------------------------------------------------------
#----------------------------------------------------------------------

outlierC.ap <- function(X, Y, NDEL = 3, NLEN = 25, TRANSPOSE = FALSE)  
{
    # Tjek af input
    if (TRANSPOSE) {
        X <- t(X)
        Y <- t(Y)
    }
    X <- tjek_data(X)
    Y <- tjek_data(Y)
    if (dim(X)[1] != dim(Y)[1]) stop("Number of firms differ in 'X' and 'Y'")

    # Data til outlierCpp
    xy <- cbind(X,Y)
    K <- dim(X)[1]   # number of firms/observations
    R <- min(K, NDEL)   # hoejeste antal udeladte units
    nlen = min(NLEN, K) # gem 25 mindste vaerdier af RX

    # Variabler til at holde output
    ratio <- array(Inf, dim = c(nlen, R))
    imat <- array(Inf, dim = c(R, R))   # Rcpp kan ikke haandtere 'NA' saa
                                        # 'Inf' bruges i stedet.
    rmin <- array(Inf, dim=R)

    outlierCpp(K, R, xy, ratio, imat, rmin)
    imat[is.infinite(imat)] <- NA

    return(list(ratio=ratio, imat=imat, r0=rmin))
}
    
#----------------------------------------------------------------------
#----------------------------------------------------------------------
