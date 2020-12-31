#' Plot method for GeDS objects. Plots GeDS fits.
#'
#' @aliases plot.GeDS
#' @param x a \code{\link{GeDS-Class}} object from which the GeDS fit(s) should be extracted.
#' @param which a numeric vector specifying the iterations of stage A for which the corresponding GeDS fits should be plotted.
#' It has to be a subset of  \code{1:nrow(x$stored)}. See details.
#' @param DEV logical variable specifying whether a plot representing the deviance at each iteration of stage A should be produced or not.
#' @param ask logical variable specifying whether the user should be prompted before changing the plot page.
#' @param main optional character string to be used as a title of the plot.
#' @param legend.pos the position of the legend within the panel. See \link[graphics]{legend} for details.
#' @param new.window logical variable specifying whether the plot should be shown in a new window or in the active one.
#' @param wait time, in seconds, the system should wait before plotting a new page.
#' Ignored if \code{ask = TRUE}.
#' @param ... further arguments to be passed to the \code{\link[graphics]{plot.default}} function.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit that should be plotted.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param type character string specifying the type of plot required.
#' Should be set either to "\code{Polygon}" if the
#' user wants to get also the control polygon of the GeDS fit,  \code{"NCI"} or  \code{"ACI"} if 95\%
#' confidence bands for the predictions should be plotted (see details)
#' or \code{"none"} if only the fitted GeDS curve should be plotted.
#' Applies only when plotting a univariate spline regression.
#'
#' @details This method is provided in order to allow the user to plot the GeDS  fits contained
#' in the \code{\link{GeDS-Class}} objects.
#'
#' Since in Stage A of the GeDS algorithm the knots of a linear spline fit are sequentially located, one at a time, the user may wish to visually
#' inspect this process using the argument \code{which}.
#' The latter specifies a particular iteration number (or a vector of such numbers) for which the corresponding
#' linear fit(s) should be plotted.
#' The \code{ask} and \code{wait} arguments can help the user to manage these pages.
#'
#' By means of \code{ask} the user can determine for how long each page should appear on the screen.
#' Pages are sequentially replaced by pressing the enter button.
#'
#' Note that, in order to ensure stability, if the object was produced by the function \code{\link{GGeDS}},
#' plotting intermediate fits of stage A is allowed  only if \code{n = 2}, in contrast to objects produced
#'  by  \code{\link{NGeDS}} for which plotting intermediate results is allowed also for \code{n = }2 or 3 results.
#'
#' The confidence intervals obtained by setting \code{type = "NCI"} are approximate local
#' bands obtained considering the knots as fixed constants.
#' Hence the columns of the design matrix are seen as covariates and standard
#' methodology relying on the \code{se.fit} option of \code{predict.lm} or \code{predict.glm} is applied.
#'
#' Setting \code{type = "ACI"}, asymptotic confidence intervals are plotted. This option is
#' applicable only if the canonical link function has been used in the fitting procedure.
#'
#' @seealso \code{\link{NGeDS}} and \code{\link{GGeDS}}; \code{\link[graphics]{plot}}.
#'
#'
#' @examples
#' ###################################################
#' # Generate a data sample for the response variable
#' # Y and the single covariate X, assuming Normal noise
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression using NGeDS
#' (Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#' # Plot the final quadratic GeDS fit (red solid line)
#' # with its control polygon (blue dashed line)
#' plot(Gmod)
#'
#' # Plot the quadratic fit obtained from the linear fit at the 10th
#' # iteration of stage A i.e. after 9 internal knots have been inserted
#' # by the GeDS procedure
#' plot(Gmod, which=10)
#'
#' # Generate plots of all the intermediate fits obtained
#' # by running the GeDS procedure
#' \dontrun{
#' plot(Gmod, which=1:16)
#' }
#'
#' ###################################################
#' # Generate a data sample for the response variable Y and the covariate
#' # X assuming Poisson distributed error and a log link function
#'
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- exp(f_1(X))
#' # Generate Poisson distributed Y according to the mean model
#' Y <- rpois(N,means)
#'
#' # Fit a Poisson GeDS regression model using GGeDS
#' (Gmod2 <- GGeDS(Y ~ f(X), beta = 0.2, phi = 0.995, family = poisson(),
#'                 Xextr = c(-2,2)))
#'
#' # similar plots as before, but for the linear fit
#' plot(Gmod2, n = 2)
#' plot(Gmod2, which = 10, n = 2)
#' \dontrun{
#' plot(Gmod2, which = 1:16, n = 2)
#' plot(Gmod2, which = 1:16, n = 2, ask = T)
#' }
#'
#' @export
setMethod("plot", signature(x = "GeDS"),  function(x, which, DEV = FALSE, ask = FALSE, main,
                                                   legend.pos = "topright", new.window = FALSE, wait = 0.5,n=3L,
                                                   type=c("Polygon","NCI","ACI","none"),...){

  if(length(DEV)!= 1 || length(ask)!= 1 || length(new.window)!= 1 || length(wait)!= 1  )
    stop("Please, check the length of the parameters")

  draw.legend <- !(is.na(legend.pos))

  results <- list()
  results$terms <- x$terms

  DEV <- as.logical(DEV)
  ask <- as.logical(ask)
  new.window <- as.logical(new.window)
  wait <- as.numeric(wait)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
    slp <- 0
  } else {
    slp <- wait
  }
  X <- x$Args$X
  Y <- x$Args$Y
  Z <- x$Args$Z
  q <- x$Args$q
  offset <- x$Args$offset
  weights <- x$Args$weights
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' wrongly specified. Set to 3.")}
  if(n==2L) toprint = "Linear"
  if(n==3L) toprint = "Quadratic"
  if(n==4L) toprint = "Cubic"

  maxim <- nrow(x$Stored)
  others <- list(...)



  #dare opzione linear quadratic cubic e stamparle assieme al ctrl polygon
  #bisogna vedere come la vede il persp...
  if(x$Type == "LM - biv") {
    W <- x$Args$W
    newX <- seq(from=range(X)[1],to=range(X)[2],length=100)
    newY <- seq(from=range(Y)[1],to=range(Y)[2],length=100)
    grid.data <- expand.grid(newX,newY)#expand.grid(X,Y)#

    if(n == 2L) obj <- x$Linear
    if(n == 3L) obj <- x$Quadratic
    if(n == 4L) obj <- x$Cubic
    matriceX <- splineDesign(knots=obj$Xknots,derivs=rep(0,length(grid.data[,1])),
                             x=grid.data[,1],ord=n,outer.ok = T)
    matriceY <- splineDesign(knots=obj$Yknots,derivs=rep(0,length(grid.data[,2])),
                             x=grid.data[,2],ord=n,outer.ok = T)
    matricebiv <- GeDS:::tensorProd(matriceX,matriceY)
    if(!is.null(W)){
      Z <- Z - W%*%x$Linear$Theta[-(1:dim(matricebiv)[2])]
    }
    z <- matricebiv%*%x$Quadratic$Theta[1:dim(matricebiv)[2]]
    z <- matrix(z,100,100)
    res<-persp(newX, newY, z,
               theta = 60, phi = 30, expand = 0.5,border=NULL,xlab="X",ylab="Y",zlab="Z", main=paste0( x$Nintknots$X, " X internal knots and ", x$Nintknots$Y ," Y internal knots" ),
               ...)
    tmp <- (obj$Predicted-Z>0)
    points(trans3d(X[tmp],Y[tmp] ,Z[tmp], pmat = res), col = "red", pch = 16, cex=.1)
    points(trans3d(X[!tmp],Y[!tmp] ,Z[!tmp], pmat = res), col = "green", pch = 16, cex=.1)
  } else {


    col <- if("col" %in% names(others)) others$col else "red"
    others$col = NULL

    if(missing(which)) which <- x$Nintknots+1
    if(length(which)==1) if(which== "all") which <- 1:(maxim)

    extr <- x$Args$extr
    if(x$Type == "LM - Univ"){
      type <- match.arg(type)
      xname <- attr(x$terms,"specials")$f-1
      xname <- attr(x$terms,"term.labels")[xname]
      xname <- substr(xname,3,(nchar(xname)-1))
      yname <- rownames(attr(x$terms,"factors"))[1]
      maxim <- nrow(x$Stored)
      if (!(is.numeric(which) || is.null(which)) || any(which < 1) || any(which > maxim)) {
        stop(sprintf("'which' must be between 1 and %d", maxim),domain=NA)}
      if(!is.null(which)){
        which <- as.integer(which)
        last <- max(which)} else {
          last <-  0
        }


      if(any(which > maxim - q)) {
        warning("Plotting also iterations discarded by the algorithm")
      }

      for(i in which){
        if(new.window) dev.new()
        if(slp >0) {
          Sys.sleep(slp)
        }
        main0 <- if(missing(main)) paste0( i-1, " internal knots") else main
        plot(X,Y,main = main0,xlab=xname,ylab=yname,...)
        ik <- x$Stored[i,-c(1,2,(i+2),(i+3))]
        knt <- makenewknots(ik,n)

        temp<-SplineReg_LM(X=X,Y=Y,Z=Z,offset = offset,weights=weights,extr=NULL,
                           InterKnots=sort(c(knt,rep(extr,n))),n=n)
        if(length(unique(X))<15*length(knt) || length(X)< 100) {
          step <- (range(extr)[2]-range(extr)[1])/(8*length(knt))
          step <- rep(step,(15*length(knt)))
          XtoPlot <- min(extr)+c(0,cumsum(step))
          matrice <- splineDesign(knots=sort(c(knt,rep(extr,n))),derivs=rep(0,length(XtoPlot)),x=XtoPlot,ord=n,outer.ok = T)

          temp$Predicted <- cbind(matrice,Z)%*%temp$Theta


          X <- XtoPlot
        }
        results$X <- X
        results$predicted <- temp$Predicted


        lines(X,temp$Predicted,col=col)
        rug(c(knt,rep(extr,n)))
        if(type=="Polygon"){
          results$Polykn <- temp$Poly$Kn
          results$Polyth <- temp$Poly$Thetas

          lines(temp$Poly$Kn,temp$Poly$Thetas,col="blue",lty=2)
          points(temp$Poly$Kn,temp$Poly$Thetas,col="blue")
          if(draw.legend) legend(legend.pos,c("Data",toprint,"Polygon"),lty=c(NA,1,2),col=c("black","red","blue"),pch=c(1,NA,1))
        } else {
          if(type=="NCI"){
            lines(X,temp$NCI$Upp,col="grey",lty=2)
            lines(X,temp$NCI$Low,col="grey",lty=2)
            results$CIupp <- temp$NCI$Upp
            results$CIlow <- temp$NCI$Low
            if(draw.legend) legend(legend.pos,c("Data",toprint,"CI"),lty=c(NA,1,2),col=c("black","red","grey"),pch=c(1,NA,NA))
          } else {
            if (type=="ACI"){
              lines(X,temp$ACI$Upp,col="grey",lty=2)
              lines(X,temp$ACI$Low,col="grey",lty=2)
              results$CIupp <- temp$ACI$Upp
              results$CIlow <- temp$ACI$Low
              if(draw.legend) legend(legend.pos,c("Data",toprint,"CI"),lty=c(NA,1,2),col=c("black","red","grey"),pch=c(1,NA,NA))
            } else
              if(draw.legend) legend(legend.pos,c("Data",toprint),lty=c(NA,1),col=c("black","red"),pch=c(1,NA))
          }
        }
      }

      if(DEV){
        main0 <- if(missing(main)) "RSS" else main
        plot(1:x$iters+3,x$RSS/length(X), main = main0,xlab="Knots", ylab = "",...)
        if(slp >0) {
          Sys.sleep(slp)
        }
        main0 <- if(missing(main)) expression(sqrt(RSS)) else main
        plot(1:x$iters+3,(x$RSS/length(X))^.5, main = main0,xlab="Knots", ylab = "",...)
        if(slp >0) {
          Sys.sleep(slp)

        }
        main0 <- if(missing(main)) expression(phi) else main
        plot(1:(x$iters-q)+3,(x$RSS[(1+q):x$iters]/x$RSS[(1):(x$iters-q)]), main = main0,xlab="Knots", ylab = "",...)
        if(slp >0) {
          Sys.sleep(slp)

        }
        main0 <- if(missing(main)) expression(sqrt(phi)) else main
        plot(1:(x$iters-q)+3,(x$RSS[(1+q):x$iters]/x$RSS[(1):(x$iters-q)])^.5, main = main0,xlab="Knots", ylab = "",...)
        if(slp >0) {
          Sys.sleep(slp)

        }
      }
    } else {
      if(x$Type == "GLM - Univ") {
        family <- x$Args$family

        xname <- attr(x$terms,"specials")$f-1
        xname <- attr(x$terms,"term.labels")[xname]
        xname <- substr(xname,3,(nchar(xname)-1))

        type <- match.arg(type)
        maxim <- nrow(x$Stored)
        if (!is.numeric(which) || any(which < 1) || any(which > maxim)) {
          stop(sprintf("'which' must be between 1 and %d", maxim),domain=NA)}
        which <- as.integer(which)
        yylim <- if(length(which)>1) range(x$Linear$Predicted) else NULL
        last <- max(which)
        if(any(which > maxim - q)) {
          warning("Plotting also iterations discarded by the algorithm")
        }

        if(n!=2L && which != maxim - q) {
          which <- maxim - q
          warning("Stage A iterations can be plotted only for the linear spline")
        }
        for(i in which){
          if(new.window) dev.new()
          if(slp >0) {
            Sys.sleep(slp)
          }
          ik <- x$Stored[i,-c(1,2,(i+2),(i+3))]
          knt <- makenewknots(ik,n)#quad


          if(n!=2L){
            temp <- if(n==3L) x$Quadratic else x$Cubic
            temp2 <- temp

            if((length(unique(X))<8*length(knt) || length(X)< 100) & is.null(offset) & is.null(weights)) {
              step <- (range(extr)[2]-range(extr)[1])/(8*length(knt))
              step <- rep(step,(8*length(knt)))
              XtoPlot <- min(extr)+c(0,cumsum(step))


              X <- XtoPlot
            }
            matrice <- splineDesign(knots=sort(c(knt,rep(extr,n))),
                                    derivs=rep(0,length(X)),x=X,ord=n,outer.ok = T)

            temp$Predicted <- matrice%*%temp$Theta[1:ncol(matrice)]

          } else {

            temp <- list("Poly"=list())
            matrice <- splineDesign(knots=sort(c(knt,rep(extr,2))),
                                    derivs=rep(0,length(X)),x=X,ord=2,outer.ok = T)
            temp$Predicted <- matrice %*% x$Coefficients[i,1:(i+1)]
            temp$Poly$Kn <- sort(c(knt,extr))
            temp$Poly$Thetas <- x$Coefficients[i,1:(i+1)]
            temp2 <- temp

          }
          main0 <- if(missing(main)) paste0( i-1, " internal knots") else main

          results$X <- X
          results$pred <- temp$Predicted
          others2 <- c(list(x=X,"y"=temp$Predicted,type="l",main = main0,
                            col=col,xlab=xname,ylab="Predictor",ylim = yylim),others)
          do.call(plot.default, others2)
          rug(x$Stored[i,1:(i+3)])
          if(type=="Polygon"){
            results$Polykn <- temp$Poly$Kn
            results$Polyth <- temp$Poly$Thetas
            lines(temp$Poly$Kn,temp$Poly$Thetas,col="black",lty=2)
            points(temp$Poly$Kn,temp$Poly$Thetas,col="black")
          } else {
            if(type=="NCI"){
              yy <- Y
              xx <- if(n!=2L) temp$Basis else {
                x$Linear$Basis
              }
              temp3 <- glm(yy ~ -1+xx,offset=offset,weights=weights,family=x$Args$family, start= temp$Thetas)
              pred <- predict.glm(temp3,newdata=data.frame(xx=matrice,offset=0),se.fit=T,type="link")
              CIupp <- pred$fit+qnorm(.975)*pred$se.fit
              CIlow <- pred$fit-qnorm(.975)*pred$se.fit
              results$CIupp <- CIupp
              results$CIlow <- CIlow

              lines(X,CIupp,col="grey",lty=2)

              lines(X,CIlow,col="grey",lty=2)

            } else {
              if (type=="ACI"){

                CI <- confint.GeDS(object=x,n=n)
                CIupp <- CI[,2]
                CIlow <- CI[,1]
                results$CIupp <- CIupp
                results$CIlow <- CIlow

                lines(X,CIupp,col="grey",lty=2)

                lines(X,CIlow,col="grey",lty=2)
                #                  stop("'ACI' type confidence bands not available in the GLM framework")
              } else

                if(draw.legend) legend(legend.pos,c(toprint),lty=c(1),col=col,pch=c(NA))
            }
          }
        }

        if(DEV){
          main0 <- if(missing(main)) "DEV" else main
          plot(1:x$iters+3,x$RSS/length(X), main = main0,xlab="Knots", ylab = "",...)
          if(slp >0) {
            Sys.sleep(slp)
          }
          main0 <- if(missing(main)) expression(sqrt(DEV)) else main
          plot(1:x$iters+3,(x$RSS/length(X))^.5, main = main0,xlab="Knots", ylab = "",...)
          if(slp >0) {
            Sys.sleep(slp)
          }
          main0 <- if(missing(main)) expression(phi) else main
          plot(1:(x$iters-q)+3,(x$RSS[(1+q):x$iters]/x$RSS[(1):(x$iters-q)]), main = main0,xlab="Knots", ylab = "",...)
          if(slp >0) {
            Sys.sleep(slp)
          }
          main0 <- if(missing(main)) expression(sqrt(phi)) else main
          plot(1:(x$iters-q)+3,(x$RSS[(1+q):x$iters]/x$RSS[(1):(x$iters-q)])^.5, main = main0,xlab="Knots", ylab = "",...)
          if(slp >0) {
            Sys.sleep(slp)
          }
        }
      } else stop("Type not recognized")}
  }
  x$Plotinfo <- results
  invisible(x)
}
)




