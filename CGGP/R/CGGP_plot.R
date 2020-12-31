#' Heatmap of SG design depth
#' 
#' The values on the diagonal are largest design depth for that dimension.
#' The off-diagonal values are the largest design depth that both dimensions
#' have been measured at simultaneously.
#' A greater depth means that more points have been measured along that
#' dimension or two-dimensional subspace.
#'
#' @param CGGP CGGP object
#'
#' @return A heat map made from ggplot2
#' @export
#' @family CGGP plot functions
#' @references https://stackoverflow.com/questions/14290364/heatmap-with-values-ggplot2
#'
#' @examples
#' \donttest{
#' # All dimensions should look similar
#' d <- 8
#' SG = CGGPcreate(d,201)
#' CGGPplotheat(SG)
#' 
#' # The first and fourth dimensions are most active and will have greater depth
#' SG <- CGGPcreate(d=5, batchsize=50)
#' f <- function(x) {cos(2*pi*x[1]*3) + exp(4*x[4])}
#' for (i in 1:1) {
#'   SG <- CGGPfit(SG, Y=apply(SG$design, 1, f))
#'   SG <- CGGPappend(CGGP=SG, batchsize=200)
#' }
#' # SG <- CGGPfit(SG, Y=apply(SG$design, 1, f))
#' CGGPplotheat(SG)
#' }
CGGPplotheat <- function(CGGP) {
  skinny <- NULL
  for (i in 1:CGGP$d) {
    skinny <- rbind(skinny, c(i, i, max(CGGP$uo[,i])))
  }
  for (i in 1:(CGGP$d-1)) {
    for (j in (i+1):CGGP$d) {
      skinny <- rbind(skinny,
                      c(i, j, max(apply(CGGP$uo[,c(i,j)], 1, min))),
                      c(j, i, max(apply(CGGP$uo[,c(i,j)], 1, min)))
      )
    }
  }
  
  skdf <- data.frame(skinny)
  names(skdf) <- c('Var1', 'Var2', 'value')
  ggplot2::ggplot(skdf, ggplot2::aes_string('Var1', 'Var2')) +
    ggplot2::geom_tile(ggplot2::aes_string(fill = 'value')) + 
    ggplot2::geom_text(ggplot2::aes_string(label = 'round(value, 1)')) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    # ggplot2::scale_fill_gradient(low = "yellow", high = "red") +
    ggplot2::scale_x_continuous(breaks = 1:CGGP$d)  +
    ggplot2::scale_y_continuous(breaks = 1:CGGP$d) # labels=c() to set names
}


#' Histogram of measurements at each design depth of each input dimension
#' 
#' A greater design depth signifies a more important dimension.
#' Thus a larger right tail on the histogram are more important variables.
#'
#' @param CGGP CGGP object
#' @param ylog Should the y axis be put on a log scale?
#'
#' @return Histogram plot made using ggplot2
#' @export
#' @family CGGP plot functions
#'
#' @examples
#' \donttest{
#' # All dimensions should look similar
#' d <- 8
#' SG = CGGPcreate(d,201)
#' CGGPplothist(SG)
#' CGGPplothist(SG, ylog=FALSE)
#' 
#' # The first dimension is more active and will have greater depth
#' f <- function(x) {sin(x[1]^.6*5)}
#' SG <- CGGPcreate(d=5, batchsize=100)
#' SG <- CGGPfit(SG, apply(SG$design, 1, f))
#' SG <- CGGPappend(CGGP=SG, batchsize=1000)
#' CGGPplothist(SG)
#' }
CGGPplothist <- function(CGGP, ylog=TRUE) {
  p <- ggplot2::ggplot(reshape2::melt(data.frame(CGGP$uo), id.vars=NULL),
                       ggplot2::aes_string(x='value'))
  p <- p +ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .)
  if (ylog) {
    p <- p + ggplot2::scale_y_log10() #limits=c(.9999, NA))
  }
  p <- p + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) # rotates labels from vert to hor
  
  p
}


#' CGGP block plot
#' 
#' Plot the 2D projections of the blocks of an CGGP object.
#'
#' @param CGGP CGGP object
#' @param singleplot If only two dimensions, should a single plot be made?
#'
#' @return ggplot2 plot
#' @export
#' @family CGGP plot functions
#'
#' @examples
#' \donttest{
#' # The first and fourth dimensions are most active and will have greater depth
#' ss <- CGGPcreate(d=5, batchsize=50)
#' f <- function(x) {cos(2*pi*x[1]*3) + x[3]*exp(4*x[4])}
#' ss <- CGGPfit(ss, Y=apply(ss$design, 1, f))
#' ss <- CGGPappend(CGGP=ss, batchsize=100)
#' CGGPplotblocks(ss)
#' 
#' mat <- matrix(c(1,1,1,2,2,1,2,2,1,3), ncol=2, byrow=TRUE)
#' CGGPplotblocks(mat)
#' }
CGGPplotblocks <- function(CGGP, singleplot=TRUE) {
  
  if (inherits(CGGP, "CGGP")) {
    d <- CGGP$d
    uo <- CGGP$uo[1:CGGP$uoCOUNT,]
  } else if (is.matrix(CGGP)) {
    d <- ncol(CGGP)
    uo <- CGGP
  } else {stop("CGGPplotblocks only works on CGGP or matrix")}
  if (d>2) {singleplot <- FALSE} # Can only use singleplot in 2D
  
  alldf <- NULL
  # ds <- c(1,2)
  d1set <- if (singleplot) {1} else {1:d}
  for (d1 in d1set) {
    d2set <- setdiff(1:d, d1)
    for (d2 in d2set) {
      ds <- c(d1, d2)
      uods <- uo[,ds]
      uodsdf <- data.frame(uods)
      uods.unique <- plyr::ddply(uodsdf, c("X1", "X2"), function(x) data.frame(count=nrow(x)))
      
      gdf <- uods.unique
      gdf$xmin <- gdf$X1-1
      gdf$xmax <- gdf$X1
      gdf$ymin <- gdf$X2-1
      gdf$ymax <- gdf$X2
      gdf$d1 <- d1
      gdf$d2 <- d2
      alldf <- if (is.null(alldf)) gdf else rbind(alldf, gdf)
    }
  }
  p <- ggplot2::ggplot(alldf) +
    ggplot2::geom_rect(ggplot2::aes_string(xmin="xmin", xmax="xmax",
                                           ymin="ymin", ymax="ymax",
                                           fill="count"),
                       color="black")+
    ggplot2::scale_fill_gradient(low = "yellow", high = "red")
  
  if (!singleplot) {p <- p + ggplot2::facet_grid(d1 ~ d2)}
  else {p <- p+ggplot2::xlab("X1") + ggplot2::ylab("X2")}
  p
}



#' Plot validation prediction errors
#'
#' @param predmean Predicted mean
#' @param predvar Predicted variance
#' @param Yval Y validation data
#' @param d If output is multivariate, which column to use. Will do all if
#' left as NULL.
#'
#' @return None, makes a plot
#' @export
#' @importFrom graphics plot points polygon
#'
#' @examples
#' x <- matrix(runif(100*3), ncol=3)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(x, 1, f1)
#' # Create a linear model on the data
#' mod <- lm(y ~ ., data.frame(x))
#' # Predict at validation data
#' Xval <- matrix(runif(3*100), ncol=3)
#' mod.pred <- predict.lm(mod, data.frame(Xval), se.fit=TRUE)
#' # Compare to true results
#' Yval <- apply(Xval, 1, f1)
#' valplot(mod.pred$fit, mod.pred$se.fit^2, Yval=Yval)
valplot <- function(predmean, predvar, Yval, d=NULL) {
  
  if (!is.null(d)) {
    predmean <- predmean[,d]
    predvar  <- predvar[,d]
    Yval <- Yval[,d]
  }
  errmax <- max(sqrt(predvar), abs(predmean - Yval))
  
  errs <- unname(predmean-Yval)
  psds <- unname(sqrt(predvar))
  if (!is.matrix(errs) || ncol(errs)==1) { # vector, single output
    tdf <- data.frame(err=errs, psd=psds)
  } else { # multiple outputs, need to melt
    tdf <- cbind(reshape2::melt(errs), reshape2::melt(psds))[,c(3,6,2)]
    tdf <- data.frame(tdf)
    names(tdf) <- c("err", "psd", "outputdim")
  }
  # ggplot(tdf, aes(x=err, y=psd)) + geom_point()
  values <- data.frame(id=factor(c(2,1)), Error=factor(c('68%','95%')))
  positions <- data.frame(id=rep(values$id, each=3),
                          x=1.1*c(0,errmax*2,-errmax*2, 0,errmax,-errmax),
                          y=1.1*c(0,errmax,errmax,0,errmax,errmax))
  # Currently we need to manually merge the two together
  datapoly <- merge(values, positions, by = c("id"))
  
  # ggplot(datapoly, aes(x = x, y = y)) +
  # geom_polygon(aes(fill = value, group = id))
  # ggplot(tdf, aes(x=err, y=psd)) + geom_polygon(aes(fill = value, group = id, x=x, y=y), 
  #              datapoly, alpha=.2) + geom_point() +
  # xlab("Predicted - Actual") + ylab("Predicted error") + 
  #   coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,errmax))
  p <- ggplot2::ggplot(tdf, ggplot2::aes_string(x='err', y='psd')) + 
    ggplot2::geom_polygon(ggplot2::aes_string(fill = 'Error', group = 'id', x='x', y='y'),
                          datapoly[datapoly$id==2,], alpha=.2) + # Separate these so it does the right order
    ggplot2::geom_polygon(ggplot2::aes_string(fill = 'Error', group = 'id', x='x', y='y'),
                          datapoly[datapoly$id==1,], alpha=.2) + 
    ggplot2::geom_point() +
    ggplot2::xlab("Predicted - Actual") + ggplot2::ylab("Predicted error") + 
    ggplot2::coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,max(psds))) #errmax))
  if (is.matrix(errs) && ncol(errs) > 1) {
    p <- p + ggplot2::facet_wrap("outputdim")
  }
  p
  
}



#' Plot validation prediction errors for CGGP object
#'
#' @param CGGP CGGP object that has been fitted
#' @param Xval X validation data
#' @param Yval Y validation data
#' @param d If output is multivariate, which column to use. Will do all if
#' left as NULL.
#'
#' @return None, makes a plot
#' @export
#' @family CGGP plot functions
#' @importFrom graphics plot points polygon
#'
#' @examples
#' SG <- CGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(SG$design, 1, f1)
#' SG <- CGGPfit(SG, y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' CGGPvalplot(CGGP=SG, Xval=Xval, Yval=Yval)
CGGPvalplot <- function(CGGP, Xval, Yval, d=NULL) {
  ypred <- CGGPpred(CGGP=CGGP, xp=Xval)
  valplot(ypred$mean, ypred$var, Yval, d=d)
}


#' Calculate stats for prediction on validation data
#'
#' @param predmean Predicted mean
#' @param predvar Predicted variance
#' @param Yval Y validation data
#' @param bydim If multiple outputs, should it be done separately by dimension?
#' @param metrics Optional additional metrics to be calculated. Should have
#'  same first three parameters as this function.
#' @param min_var Minimum value of the predicted variance.
#' Negative or zero variances can cause errors.
#' @param RMSE Should root mean squared error (RMSE) be included?
#' @param score Should score be included?
#' @param CRPscore Should CRP score be included?
#' @param coverage Should coverage be included?
#' @param R2 Should R^2 be included?
#' @param corr Should correlation between predicted and true mean be included?
#' @param MAE Should mean absolute error (MAE) be included?
#' @param MIS90 Should mean interval score for 90\% confidence be included?
#' See Gneiting and Raftery (2007).
#'
#' @return data frame
#' @export
#' @importFrom stats pnorm dnorm cor
#' @references Gneiting, Tilmann, and Adrian E. Raftery.
#' "Strictly proper scoring rules, prediction, and estimation."
#' Journal of the American Statistical Association 102.477 (2007): 359-378.
#'
#' @examples
#' valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9))
#' valstats(cbind(c(0,1,2), c(1,2,3)),
#'          cbind(c(.01,.01,.01),c(.1,.1,.1)),
#'          cbind(c(0,1.1,1.9),c(1,2,3)))
#' valstats(cbind(c(0,1,2), c(8,12,34)),
#'          cbind(c(.01,.01,.01),c(1.1,.81,1.1)),
#'          cbind(c(0,1.1,1.9),c(10,20,30)), bydim=FALSE)
#' valstats(cbind(c(.8,1.2,3.4), c(8,12,34)),
#'          cbind(c(.01,.01,.01),c(1.1,.81,1.1)),
#'          cbind(c(1,2,3),c(10,20,30)), bydim=FALSE)
valstats <- function(predmean, predvar, Yval, bydim=TRUE,
                     RMSE=TRUE, score=TRUE, CRPscore=TRUE, coverage=TRUE,
                     corr=TRUE, R2=TRUE, MAE=FALSE, MIS90=FALSE,
                     metrics,
                     min_var=.Machine$double.eps) {
  if (missing(predvar) || is.null(predvar)) {
    predvar <- NaN * predmean
  }
  # Boost up all variances to at least min_var, should be tiny number
  if (is.numeric(min_var) && !is.na(min_var) && min_var>=0) {
    predvar <- pmax(predvar, min_var)
  }
  if ((is.matrix(predmean) && 
       (any(dim(predmean)!=dim(predvar)) || any(dim(predmean)!=dim(Yval)))
  ) ||
  (length(predmean)!=length(predvar) || length(predmean)!=length(Yval))
  ) {
    stop("Shape of predmean, predvar, and Yval must match in valstats")
  }
  if (bydim && is.matrix(predmean) && ncol(predmean) > 1) {
    return(do.call("rbind",
                   lapply(1:ncol(predmean),
                          function(i) {
                            valstats(predmean=predmean[,i], predvar=predvar[,i], Yval=Yval[,i])
                          })))
  }
  
  m <- predmean
  v <- pmax(predvar, 0)
  s <- sqrt(v)
  z <- (Yval - m) / s
  
  # Create out df, add desired elements
  out <- data.frame(DELETETHISELEMENT=NaN)
  if (RMSE) out$RMSE <- sqrt(mean((predmean - Yval)^2))
  if (score) out$score <- mean((Yval-predmean)^2/predvar+log(predvar))
  if (CRPscore) out$CRPscore <- - mean(s * (1/sqrt(pi) - 2*dnorm(z) - z * (2*pnorm(z) - 1)))
  if (coverage) out$coverage <- mean((Yval<= predmean+1.96*sqrt(predvar)) & 
                                       (Yval>= predmean-1.96*sqrt(predvar)))
  if (corr) out$corr <- cor(c(predmean), c(Yval))
  if (R2) out$R2 <- 1 - (sum((Yval - predmean)^2) / sum((Yval - mean(Yval))^2))
  if (MAE) out$MAE <- mean(abs(predmean - Yval))
  if (MIS90) {
    out$MIS90 <- mean(3.28 * s + 20 * pmax(0, m - Yval - 1.64 * s) + 20 * pmax(0, -m + Yval - 1.64 * s))
  }
  
  # Remove initial element
  out$DELETETHISELEMENT <- NULL
  
  # Add normalized RMSE, dividing MSE in each dimension by variance
  if (is.matrix(predmean) && ncol(predmean) > 1) {
    # out$RMSEnorm <- sqrt(mean(sweep((predmean - Yval)^2, 2, apply(Yval, 2, var), "/")))
    out$RMSEnorm <- sqrt(mean(colMeans((predmean - Yval)^2) / apply(Yval, 2, var)))
  }
  
  # Check if user wants more metrics
  if (!missing(metrics)) {
    if (is.function(metrics)) {
      out$metric <- metrics(predmean, predvar, Yval)
    } else if (is.list(metrics)) {
      metrics.names <- names(metrics)
      if (is.null(metrics.names)) {
        metrics.names <- paste0("metrics", 1:length(metrics))
      }
      for (iii in 1:length(metrics)) {
        out[[metrics.names[iii]]] <- metrics[[iii]](predmean, predvar, Yval)
      }
    }
  }
  
  out
}

#' Calculate stats for CGGP prediction on validation data
#'
#' @param CGGP CGGP object
#' @param Xval X validation matrix
#' @param Yval Y validation data
#' @param bydim If multiple outputs, should it be done separately by dimension?
#' @param ... Passed to valstats, such as which stats to calculate.
#'
#' @return data frame
#' @export
#'
#' @examples
#' \donttest{
#' SG <- CGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(SG$design, 1, f1)
#' SG <- CGGPfit(SG, y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' CGGPvalstats(CGGP=SG, Xval=Xval, Yval=Yval)
#' 
#' # Multiple outputs
#' SG <- CGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
#' y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
#' y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
#' y <- cbind(y1, y2)
#' SG <- CGGPfit(SG, Y=y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- cbind(apply(Xval, 1, f1),
#'               apply(Xval, 1, f2))
#' CGGPvalstats(SG, Xval, Yval)
#' CGGPvalstats(SG, Xval, Yval, bydim=FALSE)
#' }
CGGPvalstats <- function(CGGP, Xval, Yval, bydim=TRUE, ...) {
  # Make predictions
  ypred <- CGGPpred(CGGP=CGGP, xp=Xval)
  # Use valstats to get df with values
  valstats(ypred$mean, ypred$var, Yval=Yval, bydim=bydim, ...)
}

#' Plot correlation samples
#' 
#' Plot samples for a given correlation function and parameters.
#' Useful for getting an idea of what the correlation parameters mean
#' in terms of smoothness.
#'
#' @param Corr Correlation function or CGGP object.
#' If CGGP object, it will make plots for thetaMAP,
#' the max a posteriori theta.
#' @param theta Parameters for Corr
#' @param numlines Number of sample paths to draw
#' @param zero Should the sample paths start at y=0?
#' @param outdims Which output dimensions should be used?
#'
#' @return Plot
#' @export
#' @family CGGP plot functions
#' @importFrom graphics par
#'
#' @examples
#' \donttest{
#' CGGPplotcorr()
#' CGGPplotcorr(theta=c(-2,-1,0,1))
#' 
#' SG <- CGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+sin(2*pi*x[2]*3)}
#' y <- apply(SG$design, 1, f)
#' SG <- CGGPfit(SG, Y=y)
#' CGGPplotcorr(SG)
#' }
CGGPplotcorr <- function(Corr=CGGP_internal_CorrMatGaussian, theta=NULL,
                         numlines=20,
                         outdims=NULL,
                         zero=TRUE) {
  # Points along x axis
  n <- 100
  xl <- seq(0,1,l=n)
  
  if (inherits(Corr, "CGGP")) {
    if (is.null(theta)) {theta <- Corr$thetaMAP}
    Corr <- Corr$CorrMat
  }
  
  nparam <- Corr(return_numpara=TRUE)
  if (is.null(theta)) {theta <- rep(0, nparam)}
  
  thetadims <- if (is.matrix(theta)) {ncol(theta)} else {1}
  if (is.null(outdims)) {
    outdims <- 1:thetadims
  } else {
    theta <- theta[,outdims, drop=FALSE]
  }
  numoutdims <- length(outdims)
  
  thetadims <- if (is.matrix(theta)) {ncol(theta)} else {1}
  
  numindims <- length(theta) / nparam / numoutdims
  indims <- 1:numindims
  ggdf <- NULL
  for (indim in indims) {
    for (outdim in outdims) {
      theta.thisloop <- if (thetadims>1) {theta[1:nparam + nparam*(indim-1), outdim]} 
      else {theta[1:nparam + nparam*(indim-1)]}
      # Can change px.mean and px.cov to be conditional on other data
      px.mean <- rep(0,n)
      px.cov <- Corr(xl, xl, theta=theta.thisloop)
      # Generate sample paths
      samplepaths <- newy <- MASS::mvrnorm(n=numlines, mu=px.mean, Sigma=px.cov)
      
      # Set so all start at 0.
      if (zero) {
        samplepathsplot <- sweep(samplepaths, 1, samplepaths[,1])
      }
      
      # Make data correct shape and add
      newdf <- cbind(reshape2::melt(data.frame(t(samplepathsplot)), id.vars=c()),
                     x=rep(xl, numlines), d=paste0("X",indim), outdim=paste0("Y",outdim))
      if (is.null(ggdf)) {ggdf <- newdf}
      else {ggdf <- rbind(ggdf, newdf)}
      
    }
  }
  
  # Return plot
  p <- ggplot2::ggplot(ggdf, 
                       ggplot2::aes_string(x="x", y="value", color="variable")) + 
    ggplot2::geom_line() + ggplot2::theme(legend.position="none")
  if (numindims > 1 && numoutdims==1) {
    p <- p + ggplot2::facet_grid(d ~ .)
  } else if (numindims > 1 && numoutdims > 1) {
    p <- p + ggplot2::facet_grid(d ~ outdim)
  }
  p
}


#' CGGP slice plot
#' 
#' Show prediction plots when varying over only one dimension.
#' Most useful when setting all values to 0.5 because it will
#' have the most points.
#'
#' @param CGGP  CGGP object
#' @param proj Point to project onto
#' @param np Number of points to use along each dimension
#' @param color Color to make error region
#' @param outdims If multiple outputs, which of them should be plotted?
#' @param scales Parameter passed to ggplot2::facet_grid()
#' @param facet If "grid", will use ggplot2::facet_grid(), if "wrap" will
#' use ggplot2::facet_wrap(). Only applicable for a single output dimension.
#'
#' @return ggplot2 object
#' @export
#' @family CGGP plot functions
#'
#' @examples
#' \donttest{
#' d <- 5
#' f1 <- function(x){x[1]+x[2]^2 + cos(x[3]^2*2*pi*4) - 3.3}
#' s1 <- CGGPcreate(d, 200)
#' s1 <- CGGPfit(s1, apply(s1$design, 1, f1))
#' #s1 <- CGGPappend(s1, 200)
#' #s1 <- CGGPfit(s1, apply(s1$design, 1, f1))
#' CGGPplotslice(s1)
#' CGGPplotslice(s1, 0.)
#' CGGPplotslice(s1, s1$design[nrow(s1$design),])
#' }
CGGPplotslice <- function(CGGP, proj=.5, np=300, color="pink", outdims, scales="free_y", facet="grid") {
  if (!is.null(CGGP$design_unevaluated)) {stop("CGGP must be updated with all data")}
  if (length(proj) == 1) {proj <- rep(proj, CGGP$d)}
  if (length(proj) != CGGP$d) {stop("proj should be of length CGGP$d or 1")}
  d <- CGGP$d
  
  tdfall <- NULL
  pointdfall <- NULL
  
  Y <- as.matrix(CGGP$Y)
  if (missing(outdims)) {
    numoutdims <- ncol(Y)
    outdims <- 1:numoutdims
  } else {
    numoutdims <- length(outdims)
  }
  
  for (d5 in 1:d) {
    xl <- seq(0,1,l=np)
    m <- matrix(proj,np,d, byrow=T)
    m[,d5] <- xl
    # Won't work if it used PCA or shared parameters
    
    p5 <- try(CGGPpred(CGGP, m, outdims=outdims), silent = TRUE)
    if (inherits(p5, "try-error")) {
      p5 <- CGGPpred(CGGP, m)
    }
    # p5 %>% str
    # poly <- cbind(c(rep(xl,each=2))[-c(1,2*np)])
    # If multiple outputs, get all of them
    for (outdim in outdims) {
      if (is.matrix(p5$mean)) { #numoutdims > 1) {
        tdf <- data.frame(mean=p5$mean[,outdim], var=p5$var[,outdim])
      } else { # Single out dim
        tdf <- as.data.frame(p5)
      }
      tdf$outdim <- outdim
      tdf$var <- pmax(0, tdf$var) # No negative values
      tdf$sd <- sqrt(tdf$var)
      tdf$meanp2sd <- tdf$mean + 2*tdf$sd
      tdf$meanm2sd <- tdf$mean - 2*tdf$sd
      tdf$x <- xl
      tdf$d <- d5
      tdfall <- rbind(tdfall, tdf)
    }
    
    w2.5 <- apply(CGGP$design[,-d5, drop=FALSE], 1, function(x) all(abs(x - proj[-d5]) < 1e-8))
    x2.5 <- CGGP$design[w2.5,, drop=FALSE]
    y2.5 <- Y[w2.5,, drop=FALSE]
    # plot(x2.5[,d5], y2.5)
    pointdf <- NULL
    if (length(y2.5) > 0) {
      for (outdim in outdims) {
        pointdf <- rbind(pointdf, data.frame(x=x2.5[,d5], y=y2.5[,outdim], d=d5, outdim=outdim))
      }
    }
    pointdfall <- rbind(pointdfall, pointdf)
  }
  
  tdfall$d <- paste0("X", tdfall$d)
  tdfall$outdim <- paste0("Y", tdfall$outdim)
  if (!is.null(pointdfall)) {
    pointdfall$d <- paste0("X", pointdfall$d)
    pointdfall$outdim <- paste0("Y", pointdfall$outdim)
  }
  p <- ggplot2::ggplot(tdfall, ggplot2::aes_string(x='x')) + 
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin='meanm2sd', ymax='meanp2sd'), color=color, fill=color) +
    ggplot2::geom_line(ggplot2::aes_string(y='mean'))
  if (!is.null(pointdfall)) {
    p <- p + ggplot2::geom_point(ggplot2::aes_string(x='x', y='y'), data=pointdfall)
  }
  if (numoutdims == 1) {
    if (facet == "grid") {
      p <- p + ggplot2::facet_grid(d ~ .)
    } else if (facet == "wrap") {
      p <- p + ggplot2::facet_wrap(d ~ .)
    } else {
      stop("Not valid facet argument in CGGPplotslice")
    }
  } else {
    p <- p +ggplot2::facet_grid(outdim ~ d, scales=scales) #"free_y")
  }
  p
}


#' Plot something similar to a semivariogram
#' 
#' It's not actually a variogram or semivariogram.
#' It shows how the correlation function falls off as distance increases.
#'
#' @param CGGP CGGP object
#' @param facet How should the plots be faceted? If 1, in a row,
#' if 2, in a column, if 3, wrapped around.
#' @param outdims Which output dimensions should be shown.
#'
#' @return ggplot2 object
#' @export
#' @family CGGP plot functions
#'
#' @examples 
#' SG <- CGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+x[3]^.4*sin(2*pi*x[2]^2*3) + .1*exp(3*x[3])}
#' y <- apply(SG$design, 1, f)
#' SG <- CGGPfit(SG, Y=y)
#' CGGPplotvariogram(SG)
CGGPplotvariogram <- function(CGGP, facet=1, outdims=NULL) {
  vdf <- NULL
  xl <- seq(0,1,l=101)
  if (is.null(outdims)) {
    outdims <- if (is.matrix(CGGP$Y)) {1:ncol(CGGP$Y)} else {1}
  }
  for (di in 1:CGGP$d) {
    for (outdim in outdims) {
      theta.thisiter <- if (is.matrix(CGGP$thetaMAP)) {CGGP$thetaMAP[1:CGGP$numpara + (di-1)*CGGP$numpara, outdim]
      } else {
        CGGP$thetaMAP[1:CGGP$numpara + (di-1)*CGGP$numpara]
      }
      tmp <- c(CGGP$CorrMat(c(0), xl, theta.thisiter))
      tmpdf <- data.frame(x=xl, y=tmp, d=paste0("X",di), outdim=paste0("Y", outdim))
      vdf <- if (is.null(vdf)) tmpdf else rbind(vdf, tmpdf)
    }
  }
  p <- ggplot2::ggplot(vdf, ggplot2::aes_string(x='x', y='y')) + ggplot2::geom_line() + ggplot2::ylim(c(0,1))
  # p <- p + facet_grid(d ~ .)
  if (facet==1) {
    p <- p + ggplot2::facet_grid(outdim ~ d)
    p <- p + ggplot2::scale_x_continuous(breaks=c(0,1))
  } else if (facet==2) {
    p <- p + ggplot2::facet_grid(d ~ outdim)
  } else if (facet==3) {
    p <- p + ggplot2::facet_wrap(d ~ outdim)
  } else {
    stop("facet is not 1, 2, or 3")
  }
  p
}


#' Plot theta samples
#'
#' @param CGGP CGGP object
#'
#' @return ggplot2 object
#' @export
#' @family CGGP plot functions
#'
#' @examples
#' gs <- CGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+x[3]^.4*sin(2*pi*x[2]^2*3) + .1*exp(3*x[3])}
#' y <- apply(gs$design, 1, f)
#' gs <- CGGPfit(gs, Y=y)
#' CGGPplottheta(gs)
CGGPplottheta <- function(CGGP) {
  # stripchart(data.frame(t(CGGP$thetaPostSamples)), xlim=c(-1,1))
  # stripchart(data.frame(t(CGGP$thetaMAP)), add=T, col=2, pch=17)
  if (is.matrix(CGGP$thetaMAP)) {
    tsamp <- NULL
    tmap <- NULL
    for (i in 1:dim(CGGP$thetaPostSamples)[3]) {
      tsamp <- rbind(tsamp,
                     reshape2::melt(CGGP$thetaPostSamples[,,i]), pdim=i)
      tmap <- rbind(tmap,
                    data.frame(value=CGGP$thetaMAP[,i], Var1=1:nrow(CGGP$thetaMAP), pdim=i))
    }
  } else { # single output
    tsamp <- reshape2::melt(CGGP$thetaPostSamples)
    tmap <- data.frame(value=CGGP$thetaMAP, Var1=1:length(CGGP$thetaMAP))
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data=tmap, mapping=ggplot2::aes_string("value", "Var1"), color="green", size=5) + 
    ggplot2::geom_point(data=tsamp, mapping=ggplot2::aes_string("value", "Var1"))
  if (is.matrix(CGGP$thetaMAP)) {
    p <- p + ggplot2::facet_grid(. ~ pdim)
  }
  p
}


#' Plot negative log posterior likelihood of samples
#'
#' @param CGGP CGGP object
#'
#' @return ggplot2 object
#' @export
#' @family CGGP plot functions
#'
#' @examples
#' gs <- CGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+x[3]^.4*sin(2*pi*x[2]^2*3) + .1*exp(3*x[3])}
#' y <- apply(gs$design, 1, f)
#' gs <- CGGPfit(gs, Y=y)
#' CGGPplotsamplesneglogpost(gs)
CGGPplotsamplesneglogpost <- function(CGGP) {
  # Number of output parameter dimensions, one plot for each
  nopd <- if (is.matrix(CGGP$thetaPostSamples)) {1}
  else {dim(CGGP$thetaPostSamples)[3]}
  
  if (nopd==1) {
    neglogpost_thetaMAP <- CGGP_internal_neglogpost(CGGP$thetaMAP,
                                                    CGGP, y=CGGP$y,
                                                    Xs=CGGP$Xs, ys=CGGP$ys)
  } else {
    neglogpost_thetaMAP <- apply(CGGP$thetaMAP, 2, CGGP_internal_neglogpost,
                                 CGGP, y=CGGP$y,
                                 Xs=CGGP$Xs, ys=CGGP$ys)
  }
  
  over_dim <- if (nopd == 1) {2} else {2:3}
  neglogpost_samples <- apply(CGGP$thetaPostSamples, over_dim,
                              CGGP_internal_neglogpost,
                              CGGP=CGGP, y=CGGP$y,
                              Xs=CGGP$Xs, ys=CGGP$ys)
  num_Inf <- sum(is.infinite(neglogpost_samples))
  if (num_Inf > 0) {
    warning(paste(num_Inf, "neglogpost samples are Inf"))
  }
  nlps_melt <- reshape2::melt(neglogpost_samples)
  nlps_melt$Y_var2 <- paste0('Y', nlps_melt$Var2)
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_histogram(ggplot2::aes_string(x='value'),
                            nlps_melt, bins=30) + 
    
    ggplot2::xlab("neglogpost") +
    ggplot2::ggtitle("neglogpost of theta samples (blue is MAP)")
  
  if (nopd > 1) {
    vl <- data.frame(nlp=neglogpost_thetaMAP,
                     Y_var2=paste0('Y', 1:length(neglogpost_thetaMAP)))
    p <- p +
      ggplot2::geom_vline(data=vl,
                          mapping=ggplot2::aes_string(xintercept="nlp"),
                          color="blue", size=2) +
      ggplot2::facet_wrap(. ~ Y_var2)
  } else {
    vl <- data.frame(nlp=neglogpost_thetaMAP)
    p <- p +
      ggplot2::geom_vline(data=vl,
                          mapping=ggplot2::aes_string(xintercept="nlp"),
                                 color="blue", size=2)
  }
  
  p
}

#' Plot CGGP block selection over time
#' 
#' Shows the order in which blocks were selected
#' for each dimension.
#' Gives an idea of how the selections change over time.
#'
#' @param CGGP CGGP object
#' @param indims Which input dimensions should be shown?
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' gs <- CGGPcreate(d=3, batchsize=100)
#' # All dimensions will look similar
#' CGGPplotblockselection(gs)
#' \donttest{
#' # You need to append with CGGPappend after fitting to see a difference
#' f <- function(x){x[1]^1.2}
#' y <- apply(gs$design, 1, f)
#' gs <- CGGPfit(gs, Y=y)
#' gs <- CGGPappend(gs, 100)
#' # Now you will see higher for X1 from 100 to 200 while others remain low.
#' CGGPplotblockselection(gs)
#' }
CGGPplotblockselection <- function(CGGP, indims) {
  uodf <- CGGP$uo[1:CGGP$uoCOUNT,]
  if (!missing(indims) && !is.null(indims)) {
    uodf <- uodf[,indims]
  }
  tdf2 <- reshape2::melt(data.frame(uodf,
                                    ninblock=CGGP$gridsize,
                                    ncumsum=cumsum(CGGP$gridsize),
                                    ind=1:CGGP$uoCOUNT),
                         id.vars=c("ind", "ninblock", "ncumsum"))
  tdf2$Var2 <- as.integer(substr(tdf2$variable, 2, 3))
  ggplot2::ggplot(data=tdf2, mapping=ggplot2::aes_string("ncumsum", "value", weight="ninblock")) + ggplot2::geom_point() +
    ggplot2::facet_grid(Var2 ~ .) +
    ggplot2::stat_smooth(color="green", method="loess", formula = y ~ x) +
    ggplot2::xlab("uo") + ggplot2::ylab("Block level")
}