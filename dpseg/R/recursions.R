

### SCORING FUNCTIONS
## for generic recursion S[i-jumps] + score(i,j)
get_scoreij <- function(type) 
    get(paste0("scoreij","_",type), mode="function")
scoreij_var <- function(i,j,x,y,...) -var(residuals(lm(y[i:j]~x[i:j], ...)))
scoreij_cor <- function(i,j,x,y,...) abs(cor(y[i:j], x[i:j], ...))-1
scoreij_r2 <- function(i,j,x,y,...) summary(lm(y[i:j]~x[i:j], ...))$r.squared-1
scoreij_r2adj <- function(i,j,x,y,...)
    summary(lm(y[i:j]~x[i:j],...))$adj.r.squared-1
## pre-calculated matrix - used as "y"
scoreij_matrix <- function(i,j,x,y,...) y[i,j]



## generic dynamic programing recursion 
## 
## Generic dynamic programing recursion, solving
## \deqn{S_j = max_{\substack{\\i \le j-\Ell_{\min}\\i \ge j-\Ell_{\max}}} S_{i-jumps} + score(i,j) - P}{S[j] = max(S[i-jumps] + score(i,j) -P)}.
## 
## This implementation is highly inefficent when used with actual scoring
## functions, but fast when used with a pre-calculated scoring matrix,
## passed via argument \code{y}. The function is mainly used to (a) define
## input and output of recursions to be compatible with the \code{dpseg}
## wrapper (TODO: use setGeneric/setMethod), (b) to quickly test
## new scoring functionss, or (c) use with pre-calculated scoring matrices.
recursion_generic <- function(x, y, maxl, jumps=0, P=0, minl=3, S0=1,
                              type="var", storev=FALSE, storem=FALSE,
                              scoref, ...) {

    ## use pre-calculated scoring matrix
    if ( "matrix" %in% class(y) ) {
        type <- "matrix"
        scoref <- scoreij_matrix
        x <- seq_len(nrow(y))
    }
    
    ## get scoring function
    if ( missing(scoref) )
        scoref <- get_scoreij(type)

    ## ALLOW DISCRETE JUMPS
    ## jumps=TRUE -> S[i-1] + score(i,j)
    ## jumps=FALSE -> S[i] + score(i,j)
    jumps <- as.numeric(jumps)
    
    ## initialize
    N <- length(x)
    S <- imax <- numeric(N) # intialized to 0
    S[] <- -P

    if ( missing(maxl) )
        maxl <- N 
            
    ## score matrix: store for educational movies or debugging,
    ## or re-use with type="matrix"
    SCR <- NULL
    if ( storem ) 
        SCR <- matrix(NA, nrow=N, ncol=N)
    
    ## RECURSION
    for ( j in 2:N) {

        irng <- (j-maxl):(j-minl) +1 
        irng <- irng[irng>0] # skip irng<1 at start 
        si <- rep(-Inf, maxl-minl+1) # TODO: restrict to irng

        for ( i in irng ) { 

            idx <- i - (j-maxl) 
            ## get score and S[i] vector
            scr <- scoref(i,j,x,y,...)
            si[idx] <- ifelse(i==1&jumps==1, S0, S[i-jumps]) + scr - P
            
            if ( storem ) SCR[i,j] <- scr
        }
        S[j] <- max(si)
        idx <- which.max(si)
        imax[j] <- idx + (j-maxl) 
    }
    list(imax=imax, S=S, SCR=SCR)
}

## efficient lin.reg. recursion
## note: scoref ignored, only for consistency with recursion_generic
## S[i-jumps] + score(i,j) 
recursion_linreg <- function(x, y, maxl, jumps=0, P=0, minl=3, S0=1,
                             type="var", storev=TRUE, storem=FALSE,
                             scoref, ...) {
    
    ## ALLOW DISCRETE JUMPS
    ## jumps=TRUE -> S[i-1] + score(i,j)
    ## jumps=FALSE -> S[i] + score(i,j)
    jumps <- as.numeric(jumps)

    ## initialize
    N <- length(x)
    S <- imax <- numeric(N) # init to 0
    S[] <- -P #*(1:minl)
    
    if ( missing(maxl) )
        maxl <- N 

    ## store matrix: only required for educational movies!
    SCR <- NULL
    if ( storem ) 
        SCR <- matrix(NA, nrow=N, ncol=N)

    ## initialize x2, y2, xy
    sy <- sx <- sy2 <- sx2 <- sxy <- rep(0, N)
    sx[1] <- x[1]
    sy[1] <- y[1]
    sx2[1] <- x[1]^2
    sy2[1] <- y[1]^2
    sxy[1] <- y[1]*x[1]
    
    if ( storev )
        islp <- icpt <- ivar <- irsq <- rep(NA, N)

    ## RECURSION
    for ( j in 2:N ) {

        ## incremental sums 
        sx[j] <-  sx[j-1]  + x[j]        # \sum x
        sy[j] <-  sy[j-1]  + y[j]        # \sum y
        sx2[j] <- sx2[j-1] + x[j]^2      # \sum x^2
        sy2[j] <- sy2[j-1] + y[j]^2      # \sum y^2
        sxy[j] <- sxy[j-1] + x[j]*y[j]   # \sum x*y

        ## sum of squares for i=1 
        n <- j
        Syy <- (sy2[j]) - (sy[j])^2/n
        ## Sxx = \sum x^2 - (\sum x)^2/n
        Sxx <- (sx2[j]) - (sx[j])^2/n
        ## Sxy = \sum x*y - sum(x)*sum(y)/n
        Sxy <- (sxy[j]) - (sx[j])*(sy[j])/n

        ## range of i values tested
        irng <- (j-maxl):(j-minl) + 1
        irng <- irng[irng>0] # skip negative ranges at start

        ## scoring vector
        si <- rep(-Inf, maxl+1-minl) # rep(-Inf,length(irng)) #
        if ( storev ) slpi <- cpti <- vari <- r2i <- cori <- si

        ## calculate all scores i->j
        for ( i in irng ) {

            idx <- i - (j-maxl) # si vector index
            
            ## sum of squares: variance & covariance
            ## subtract previous
            ij <- i - 1 # sum of squares subtraction index
            n <- j-ij   # segment length i:j
            
            ## ij==0 (at i=1 with jumps==FALSE) handled above
            if( ij>0 ) {
                ## Syy = \sum y^2 - (\sum y)^2/n
                Syy <- (sy2[j]-sy2[ij]) - (sy[j]-sy[ij])^2/n
                ## Sxx = \sum x^2 - (\sum x)^2/n
                Sxx <- (sx2[j]-sx2[ij]) - (sx[j]-sx[ij])^2/n
                ## Sxy = \sum x*y - sum(x)*sum(y)/n
                Sxy <- (sxy[j]-sxy[ij]) - (sx[j]-sx[ij])*(sy[j]-sy[ij])/n
            } 
            
            ## only for debug, can later be removed
            position <- paste0("i=",i,", j=",j, ", idx=",idx)
            if ( i<1 ) stop("wrong i: ", position)
            if ( idx==0 ) stop("wrong index: ", position)
            if ( n<minl ) warning("too short segment: ", position)


            ## slope and intercept from i+1 to j:
            slp <- Sxy/Sxx
            
            ## scoring function: -var, r-squared-1 or Pearson-1 
            var <- -(Syy - Sxy*slp) / (n-1)
            r2 <- slp*slp*Sxx/Syy -1
            cor <- slp*sqrt(Sxx/Syy) -1 # NOTE: abs(cor)

            ## get score(i,j): vari or r2i!
            scr <- get(paste0(type), mode="numeric")

            ## S[i-jumps] + score(i,j) - P
            if ( Sxx!=0.0 ) 
                si[idx] <- ifelse((i-jumps)==0, S0, S[i-jumps]) + scr - P
            else warning("numeric cancellation: Sxx==0.0")

            ## store values
            if ( storev ) {
                ## intercept
                if ( ij>0 )
                    cpti[idx] <- ((sy[j]-sy[ij]) - slp*(sx[j]-sx[ij]))/n
                else 
                    cpti[idx] <- ((sy[j]) - slp*(sx[j]))/n
                slpi[idx] <- slp
                vari[idx] <- var
                r2i[idx] <- r2 +1
                cori[idx] <- cor +1
            }
            ## store matrix
            if ( storem ) SCR[i,j] <- scr
        }
        
        ## GET MAXIMUM SCORE
        ## S[j] = max S[i-jumps] + score(i,j)
        S[j] <- max(si)
        idx <- which.max(si)

        ## reverse max?
        ridx <- length(si) - which.max(rev(si)) + 1
        if ( idx!=ridx & j >= minl) 
            warning("multiple max at j: ", j, ",", idx,", ", ridx)

        ## ... store which i yielded maximum
        imax[j] <- idx + (j-maxl)
        
        ## store recorded values
        if ( storev ) {
            ## correct sign for variances
            vari <- -vari
            ivar[j] <- vari[idx]
            islp[j] <- slpi[idx]
            icpt[j] <- cpti[idx]
            irsq[j] <- r2i[idx]
        }
    }

    ## store recorded values
    values <- NULL
    if ( storev ) values <- list(icpt=icpt, islp=islp, ivar=ivar, irsq=irsq)

    list(imax=imax, S=S, SCR=SCR, values=values)
}


### NOTE: outdated first version, kept for demonstration purposes
## for recursion_old S[i] + score[i+jumps][j]

get_scorexy <- function(type) 
    get(paste0("scorexy","_",type), mode="function")
scorexy_var <- function(x, y, ...) -var(residuals(lm(y~x, ...)))
scorexy_cor <- function(x, y, ...) abs(cor(y, x, ...))-1
scorexy_r2 <- function(x, y, ...) summary(lm(y~x, ...))$r.squared-1
scorexy_r2adj <- function(x, y, ...) summary(lm(y~x,...))$adj.r.squared-1

## NOTE: S0 has no effect, kept for consistency
## S[i] + score(i+jumps,j) 
recursion_old <- function(x, y, maxl, jumps=0, P=0, minl=3, S0, type="var",
                          storev=FALSE, storem=FALSE, scoref, ...) {

    if ( missing(scoref) )
        scoref <- get_scorexy(type)

    ## ALLOW DISCRETE JUMPS
    ## jumps=TRUE -> score(i+1,j)
    ## jumps=FALSE -> score(i,j)
    jumps <- as.numeric(jumps)
    
    ## initialize
    N <- length(x)
    S <- imax <- numeric(N) # init to 0
    imax[] <- 1 # TODO: required?
    S[] <- -P

    if ( missing(maxl) )
        maxl <- N - jumps
            
    ## store matrix: only required for educational movies!
    SCR <- NULL
    if ( storem ) 
        SCR <- matrix(NA, nrow=N, ncol=N)

    ## RECURSION
    for ( j in 2:N) {

        irng <- (j-maxl):(j-minl) +1-jumps # shift range
        irng <- irng[irng>0] # skip irng<1 at start 
        si <- rep(-Inf, maxl-minl+1) # TODO: restrict to irng

        for ( i in irng ) { 

            idx <- i - (j-maxl) + jumps  # vector index
            scr <- scoref(x[(i + jumps):j], y[(i + jumps):j], ...)
            si[idx] <- S[i] + scr - P
            
            if ( storem ) SCR[i,j] <- scr
        }
        S[j] <- max(si)
        idx <- which.max(si)
        imax[j] <- idx + (j-maxl) - jumps
    }
    list(imax=imax, S=S, SCR=SCR)
}
