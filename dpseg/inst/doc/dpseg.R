## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, global.par=TRUE,
                      fig.path = ".", fig.pos = 'h',
                      fig.height = 2.7, fig.width = 4, fig.align = "center")
knitr::opts_knit$set(global.par = TRUE)

## redo benchmarking or use existing figure?
REDOLONG <- FALSE
times <- 100

## ---- include=FALSE-----------------------------------------------------------
par(mai=c(.5,.5,.3,.5),mgp=c(1.2,.3,0),tcl=-.25, cex.main=.75)

## ---- eval=FALSE--------------------------------------------------------------
#  simple_backtrace <- function (imax, jumps = 0) {
#      end <- length(imax) # end of the first segment = end of the data
#      ends <- c(end)      # initiate vector
#      while (end > 1) {
#          end <- imax[end] - jumps # end of previous segment
#          ends <- c(end, ends)     # note the order, new end is prepended
#      }
#      ends
#  }

## ----dpsegdemo, fig.cap="Segmentation of a growth curve (*E. coli* in M9 minimal medium) by `dpseg`. Vertical lines indicate segment starts (dashed red) and ends (solid black). The `predict` method returns the piecewise linear model (dashed yellow).\\label{fig:dpseg}"----
library(dpseg)

type <- "var"  # use the (default) scoring, -var(residuals(lm(y~x)))
jumps <- FALSE # allow discrete jumps between segments?
P <- 1e-4      # break-point penalty, use higher P for longer segments

# get example data `oddata` - bacterial growth measured as optical density OD
x <- oddata$Time
y <- log(oddata[,"A2"]) # note: exponential growth -> log(y) is linear

# NOTE: the scoring function results are stored as a matrix for re-use below
segs <- dpseg(x=x, y=y, jumps=jumps, P=P, type=type, store.matrix=TRUE)
print(segs)
plot(segs)
lines(predict(segs),lty=2, lwd=3, col="yellow")

## ---- eval=FALSE--------------------------------------------------------------
#  # use options frames, and res to control file size
#  movie(segs, format="gif", file.name="movie", path="pkg/vignettes/figures",
#        frames=seq(1,length(x),3), delay=.3,res=50)

## ----movie, echo=FALSE, fig.cap="Animation of the progress of the algorithm through the data (gray circles and right y-axis). The black vertical line is the current position $j$, and the black circles are the values of the scoring function $\\text{score}(i,j)$ for all $i<j$. The blue line was the $i_\\text{max}$ where the maximal value of the recursion $S_j$ was found. The thin red line indicates $j-\\lmin$. Colors indicate the final segmentation after backtracing. \\label{fig:movie}"----

# NOTE: results pre-calculated with above command
knitr::include_graphics(c("figures/movie.gif"))

## ---- include=FALSE-----------------------------------------------------------
par(mai=c(.5,.5,.3,.5),mgp=c(1.2,.3,0),tcl=-.25, cex.main=.75)

## -----------------------------------------------------------------------------
p <- estimateP(x=x, y=y, plot=TRUE)
plot(dpseg(x=x, y=y, jumps=jumps, P=round(p,3)))

## ----estimatesource-----------------------------------------------------------
simple_estimateP <- function (x, y, ...) {
    var(smooth.spline(x, y, ...)$y - y)
}

## ---- include=FALSE-----------------------------------------------------------
par(mai=c(.5,.5,.3,.5),mgp=c(1.2,.3,0),tcl=-.25, cex.main=.75)

## ----scanp, eval=TRUE, fig.cap="A higher $P$ will yield fewer and longer segments. $P$ should be chosen close to the optimal value of the scoring function. \\label{fig:pscan}"----

## NOTE: dpseg is slower for many segments!
sp <- scanP(x=x, y=y, P=seq(-.01,.1,length.out=50), plot=TRUE)

## ----precalcmatrix------------------------------------------------------------
## use the scoring matrix from previous run for generic recursion,
## with store.matrix=TRUE, and test different parameters.
segm <- dpseg(y=segs$SCR, jumps=jumps, P=2*P, minl=5, maxl=50) 
print(segm)

## ---- include=FALSE-----------------------------------------------------------
par(mai=c(.5,.5,.3,.5),mgp=c(1.2,.3,0),tcl=-.25, cex.main=.75)

## ----addlm, fig.cap="Compared to the first run above both, the larger $\\lmin=5$ and $P=0.0002$, contributed to get fewer segments, while the lower $\\lmax=50$ split the last segment in two."----

## add lm-based regression coefficients
segm <- addLm(segm, x=x, y=y)
plot(segm)

## ----rcppdynprogdata, warning=FALSE, message=FALSE, fig.cap="Segmentation of the test data from the RcppDynProg package with the default scoring function, negative variance of residuals (left), and with Pearson correlation (right).", fig.height = 5.4, fig.width = 8----

## example from https://winvector.github.io/RcppDynProg/articles/SegmentationL.html
set.seed(2018)
d <- data.frame(x = 0.05*(1:(300))) # ordered in x
d$y_real <- sin((0.3*d$x)^2) 
d$y <- d$y_real + 0.25*rnorm(length(d$y_real))

par(mfrow=c(2,2))
scanP(x=d$x, y=d$y, P=seq(-0.005,0,length.out=50), verb=0) # note: no messages by verb=0
scanP(x=d$x, y=d$y, P=seq(-0.05,0,length.out=50), type="cor", verb=0)
plot(dpseg(x=d$x, y=d$y, P=-.0025, verb=0)) 
plot(dpseg(x=d$x, y=d$y, type="cor", verb=0))

## ---- include=FALSE-----------------------------------------------------------
par(mfcol=c(1,1), mai=c(.5,.5,.3,.5),mgp=c(1.2,.3,0),tcl=-.25, cex.main=.75)

## ----customscoring, eval=FALSE------------------------------------------------
#  score_rsq <- function(i, j, x, y,...) summary(lm(y[i:j]~x[i:j]))$r.squared
#  segn <- dpseg(x=x, y=y, jumps=jumps, P=.99, scoref=score_rsq, add.lm=TRUE)
#  plot(segn)

## ----custom, echo=FALSE, fig.cap="Note, that as above a meaningful penalty term $P$ should be close to the optimal value of the scoring function, here $R^2=1$."----

# NOTE: since running custom scoring functions takes long,
# the results were pre-calculated for the vignette
knitr::include_graphics(c("figures/customscoring-1.png"))

## ----rcppdynprog, echo=FALSE, eval=FALSE--------------------------------------
#  yp <- RcppDynProg::piecewise_linear(x=x,y=y)
#  plot(x,y)
#  lines(x,yp, col=2, lwd=2)

## ----benchmark, echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Benchmarking of R & Rcpp implementations."----

if (!REDOLONG) {
    knitr::include_graphics("figures/benchmark-1_100.png")
}else{
    ## time recursion_generic, recursion_linreg and recursion_linreg_c
    library(microbenchmark)
    library(ggplot2)

    ## USE check_results to compare different implementations
    ## currently deactivated, since RcppDynProg is also tested
    check_results <- function(values) {
        error <- FALSE
        for ( i in 2:length(values) )
            error <- error|any(values[[1]]!=values[[i]])
        !error
    }
    
    
    ## for pre-calculated matrix
    Sj<- dpseg(x=x, y=y,  jumps=jumps, P=P, store.matrix=TRUE, verb=0)
    SCR <- Sj$SCR
    
    ## NOTE `times` defined in knitr pre-amble to speed up compilation
    ## generic_r takes very long, only for final figure
    mbm <- microbenchmark::microbenchmark(times=times,
                          "generic_r" = {
                              Sjg<- dpseg(x=x, y=y, jumps=jumps, P=P, verb=0, 
                                          recursion=dpseg:::recursion_generic,
                                          store.values=FALSE)
                              Sjg$traceback},
                          "incremental_r" = {
                              Sji<- dpseg(x=x, y=y, jumps=jumps, P=P, verb=0,
                                          recursion=dpseg:::recursion_linreg) 
                              Sji$traceback},
                          "precalculated_r" = {                          
                              Sjp <- dpseg(y=SCR, jumps=jumps, P=P, verb=0,
                                           recursion=dpseg:::recursion_generic)
                              Sjp$traceback},
                          "incremental_rcpp" = {
                              Sjc<- dpseg(x=x, y=y, jumps=jumps, P=P, verb=0,
                                          recursion=dpseg:::recursion_linreg_c) 
                              Sjc$traceback},
                          "precalculated_rcpp" = {                          
                              Sjp <- dpseg(y=SCR, jumps=jumps, P=P, verb=0)
                              Sjp$traceback},
			  "RcppDynProg" = {
			      Sjp <- RcppDynProg::piecewise_linear(x=x,y=y)
			      Sjp
			      },    
                          check = NULL) #check_results)
    
    #mbm
    ggplot2::autoplot(mbm)
    ##boxplot(mbm$time ~ mbm$expr, horizontal=TRUE, axes=FALSE)
    ##axis(2, las=2)
}

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
## RECURSION
recursion <- function(x, y,  maxl, jumps=0, P=0, minl=3, S0=1) {
  N = length(x)
  S = numeric(N) # init to 0
  imax = numeric(N)
  S[1] = -P
  for ( j in 2:N) {
    si = rep(-Inf, maxl-minl+1)
    irng = (j-maxl):(j-minl) +1
    irng = irng[irng>0]
    for ( i in irng ) { 
      idx = i-(j-maxl) 
      sij = ifelse(i==1&jumps==1, S0, S[i-jumps])
      si[idx] = sij + score(x, y, i, j) - P
    }
    S[j] = max(si)
    idx = which.max(si)
    imax[j] = idx + (j-maxl)
  }
  imax
}

## SCORING FUNCTION
score <- function(x, y, k, l) -var(residuals(lm(y[k:l]~x[k:l])))

## backtracing
backtrace <- function(imax, jumps=0) {
  end = length(imax) # end of last segment
  ends = end
  while( end>1 ) {
    end = imax[end] - jumps
    ends = c(end, ends) 
  }
  ends[1] <- ends[1] + jumps # note: start of first segment
  ends
}

## ---- include=FALSE-----------------------------------------------------------
par(mai=c(.5,.5,.3,.5),mgp=c(1.2,.3,0),tcl=-.25, cex.main=.75)

## ---- eval=TRUE, echo=TRUE, fig.cap="Correction segmentation (horizontal lines) for a primitive test case."----
# simple test case 
k1=1
k2=.05
k3=-.5
y1 <- k1*1:5
y2 <- k2*1:5 + k1*5
y3 <- k3*1:5 + k2*5 + k1*5
set.seed(1)
ym <- c(y1, y2, y3)
nsd <- .25 # noise, standard deviation
y <- ym + rnorm(length(ym), 0, nsd) # add noise
x <- 1:length(y)

## run recursion
JUMPS = 0
imax = recursion(x, y, maxl=length(x), jumps=JUMPS, P=0, minl=3, S0=1)

## backtrace
ends = backtrace(imax, jumps=JUMPS)
print(ends)

## plot
plot(x,y)
lines(x,ym)
legend("bottom", title="slopes:", legend=paste(k1,k2,k3,sep=", "), bty="n")
abline(v=ends)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
## for development only: silently reload data to test code
jumps <- FALSE # allow discrete jumps between segments?
P <- 1e-4      # break-point penalty, use higher P for longer segments

# get example data `oddata` - bacterial growth measured as optical density OD
x <- oddata$Time
y <- log(oddata[,"A2"]) # note: exponential growth -> log(y) is linear

