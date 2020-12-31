######################### SEMIVARIOGRAM FUNCTION ###########################
# LAST UPDATE: 11/15/2020

#' Semivariogram for Geospatial Data
#'
#' This function creates semivariogram plots. It creates empirical semivariogram 
#' for raw data and \code{lm} object or parametric exponential semivariogram based 
#' on the estimation from \code{metropolis.krige}. Based on the user's chosen level 
#' of coarsening, the semivariogram is presented for various distances.
#'
#' @param x An object for which a semivariogram is desired. The object can 
#'   be a \code{krige} object, a \code{semivariance} object, a \code{lm} object, 
#'   or a vector of variables (or variable names in the \code{data}).
#' @param coords A matrix of coordinates for all observations or a vector of variable 
#'   names indicating the coordinates variables in the data. Alternatively, the 
#'   coordinates can also be specified separately using \code{east} and \code{north}.
#' @param data If object is a variable name, a data frame must be provided.
#' @param bins Number of bins into which distances should be divided. The observed 
#'   distances will be split at equal intervals, and the semivariance will be computed 
#'   within each interval. Defaults to 13 intervals.
#' @param terms A vector of strings specifies for which the semivariogram is created. 
#'   Options are "raw" (the semivariogram for raw data), "residual" (the semivariogram
#'   for residuals from linear regression). "parametric" (the parametric powered 
#'   exponential semivariogram based on the estimation from \code{metropolis.krige}). 
#'   Defaults will create all the applicable plots. 
#' @param east Alternative specification for the vector of eastings for all observations.
#' @param north Alternative specification for the vector of northing for all observations.
#' @param type A vector specify the type of plots for each term. Options are "l" 
#'   (line plot) and "p" (scatter plot). Defaults to \code{c(raw = "p", residual = "p", 
#'   parametric = "l")}
#' @param legend A logical argument for whether legend should be presented. Defaults to \code{TRUE}.
#' @param col A vector specify the color for each term. Defaults to \code{c(raw = "black", 
#'   residual = "blue", parametric = "red")}
#' @param pch A vector specify the points symbols for scatter plot. Suppressed for 
#'   line plot.
#' @param lty A vector specify the line type for line plot. Suppressed for 
#'   scatter plot.
#' @param \dots Additional arguments to be passed to \code{semivariogram} methods. 
#'   Further arguments that can passed to \code{plot()} function can be specified
#'   here.
#'   
#' @return An semivariogram plot. For \code{krige} object, it will return empirical 
#'   semivariograms for raw data and residuals of linear regression as well as a 
#'   parametric powered exponential semivariogram given the values of the estimates 
#'   from \code{metropolis.krige} as default.
#'   
#' @details The function creates semivariograms for diagnosing the spatial relationship
#'   that best describes the data and for diagnosing the degree to which the model 
#'   fits the spatial relationship. With a view of the empirical semivariogram, 
#'   a user can consult images of parametric semivariograms to determine whether 
#'   an exponential, Gaussian, or other powered expoential function fit the data 
#'   well, or if another style of semivariogram works better. Examining this also 
#'   allows the user to develop priors such as the approximate split in variance 
#'   between the nugget and partial sill as well as the approximate distance of 
#'   the effective range. Semivariograms are explicitly tied to a corresponding 
#'   spatial correlation function, so determining the former automatically implies 
#'   the latter. See Banerjee, Carlin, and Gelfand for a fuller explanation, as 
#'   well as a guidebook to semivariogram diagnosis (2015, 26-30).
#'   
#'   The function can be applied to a variable, a fitted linear model (\code{lm} 
#'   object) before fitting a spatial model or to a \code{krige} object or \code{semivariance}
#'   object to assess the model fit. When applying to a variable, it will describes 
#'   the raw data; for a \code{lm} object, the default will present empirical 
#'   semivariogram for both the raw data and linear residuals; when applying to a 
#'   \code{krige} object, the default will present empirical semivariogram for the 
#'   raw data and the residuals from linear fit, and the parametric semivariogram 
#'   given the estimates from the geospatial model fitted in \code{metropolis.krige}; 
#'   for a \code{semivariance} object, it will present a plot(s) for whichever the
#'   semivariance is calculated. Users can also specify which semivariogram is needed 
#'   in the \code{terms} argument if there are multiple kinds of semivariogram can 
#'   be plotted. The plots are compatible to the arguments used in base \code{R} 
#'   base graphics. 
#'   
#' @references 
#'   Sudipto Banerjee, Bradley P. Carlin, and Alan E. Gelfand. 2015. \emph{Hierarchical 
#'   Modeling and Analysis for Spatial Data}. 2nd ed. Boca Raton, FL: CRC Press.
#'   
#' @seealso \code{\link{semivariance}}, \code{\link{exponential.semivariance}}   
#' 
#' @examples
#' \dontrun{
#' # Summarize Data
#' summary(ContrivedData)
#' 
#' # Empirical semivariagram for variable y
#' semivariogram(x=ContrivedData$y, coords = cbind(ContrivedData$s.1, ContrivedData$s.2))
#' 
#' # Initial OLS Model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' 
#' # Empirical semivariagram for ols fit
#' semivariogram(contrived.ols, coords = c("s.1","s.2"), bins=13)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, range.tol = 0.05)
#'    
#' # Parametric semivariagram
#' semivariogram(contrived.run, bins=13, terms = c("raw", "residual", "parametric"),
#'   type= c(raw = "p", residual = "p", parametric = "l"), legend = TRUE, col = c("black", 
#'   "blue", "red"), pch = c(1,3,NA), lty = c(NA,NA,1))
#'   
#' # Alternatively, the generic function plot can also be applied to krige object
#' plot(contrived.run)
#' 
#' # Plot semivariance object
#' plot(semivariance(contrived.run, bins=13))
#' }
#' 
#' @importFrom graphics plot axis lines legend points box
#' @export

semivariogram <- function(x, ...){
  UseMethod("semivariogram")
}

#' @rdname semivariogram
#' 
#' @importFrom stats as.formula model.frame
#' @export
semivariogram.krige <- function(x, ..., bins=13, terms = "all", type, pch, lty, 
                                legend, col){
  if (!inherits(x, "krige")) stop("The input x is not a 'krige' x.")
  if (terms[1] == "all") terms <- c("raw", "residual", "parametric")
  dcol <- c(raw="black", residual="blue", parametric="red")
  dtype <- c(raw="p", residual="p", parametric="l")
  dpch <- c(raw=1, residual=3, parametric=2)
  dlty <- c(raw=3, residual=2, parametric=1)
  if (missing(col)) col <- dcol[terms]
  if (is.null(names(col))) names(col) <- terms
  if (missing(type)) type <- dtype[terms]
  if (is.null(names(type))) names(type) <- terms
  type[type == "lines"] <- "l"; type[type == "points"] <- "p"
  if (missing(col)) col <- dcol[terms]
  if (is.null(names(col))) names(col) <- terms
  if (missing(type)) type <- dtype[terms]
  if (is.null(names(type))) names(type) <- terms
  if (missing(pch)) {pch <- ifelse(type == "p", dpch, NA)
  } else { pch <- ifelse(type == "p", pch, NA) }
  if (is.null(names(pch))) names(pch) <- terms
  if (missing(lty)) {lty <- ifelse(type == "l", dlty, NA)
  } else { lty <- ifelse(type == "l", lty, NA) }
  if (is.null(names(lty))) names(lty) <- terms
  if (all(is.na(pch))) pch <- NULL
  if (all(is.na(lty))) lty <- NULL
  if (length(terms) > 1 & missing(legend)) {legend <- TRUE
  } else if (length(terms) == 1 & missing(legend)) {legend <- FALSE}
  
  # Distance Matrix
  distance<-as.numeric(dist(cbind(x$model.data.list$easting,
                                  x$model.data.list$northing)))
  # Raw
  if ("raw" %in% terms) {
    raw <- sv(x = x$model.data.list$y, distance = distance, bins = bins, ...)
  }
  # Residual
  if ("residual" %in% terms) {
    residual <- sv(x = x$init.ols$residuals, distance = distance, bins = bins, ...)
  }
  # Parametric exponential semivariogram from estimation
  if ("parametric" %in% terms) {
    if (!exists("resid")) {residual <- sv(x = x$init.ols$residuals, 
                                          distance = distance, ...)}
    var.terms<-apply(x$mcmc.mat[,1:3],2,quantile,0.5)
    parametric <- exponential.semivariance.default(nugget=var.terms[1], decay=var.terms[2], 
                                            partial.sill = var.terms[3], 
                                            distance = as.numeric(names(residual)), 
                                            power = x$standing.parameter$powered.exp)
  }
  
  for (i in terms) {
    if (i == terms[1]) {
      sv.plot(get(i),distance=distance,bins=bins,add=FALSE,type=type[i],col=col[i],
              pch=pch[i],lty=lty[i])
    } else {
      sv.plot(get(i),distance=distance,bins=bins,type=type[i],col=col[i],pch=pch[i],
              lty=lty[i],add=TRUE)
    }
  }
  if (legend == TRUE) {
    lnames <- c(raw = "Raw data", residual = "OLS residuals", 
                parametric = "Parametric Semivariogram")
    legend("top", legend=lnames[terms], pch = as.vector(pch[terms]), inset=c(0,-.15), xpd=TRUE, 
           bty = "n",col=col[terms], lty=as.vector(lty[terms]), cex=0.8,ncol=length(terms)) }
}


#' @rdname semivariogram
#' @export
plot.krige <- function(...) semivariogram.krige(...)

#' @rdname semivariogram
#' @export
semivariogram.lm <- function(x, ..., coords, bins = 13, terms = c("raw", "residual"),
                             east, north, type, legend, col, pch, lty) {
  if (!inherits(x, "lm")) stop("The input x is not a 'lm' x.")
  if (missing(coords)) {
    if (missing(east) | missing(north)) stop("Coordinates are missing.")
    coords <- cbind(east, north)
    if (is.character(coords) & length(coords) == 2) coords <- as.vector(coords)
  }
  cl <- x$call
  y <- x$model[1]
  if (is.character(coords) & length(coords) == 2) {
    form <- update(as.formula(cl$formula), paste("~ . +",paste(coords, collapse=" + ")))
    na.action <- ifelse("na.action" %in% names(cl), cl$na.action, "na.omit")
    mf <- model.frame(form, data = get(as.character(cl$data)), na.action=na.action)
    if (length(mf[1]) != length(y)) warning("Additional missing data are dropped.")
    y <- mf[1]; coords <- mf[coords]
    x <- lm(cl$formula, data = mf)
  }
  dcol <- c(raw="black", residual="blue")
  dtype <- c(raw="p", residual="p")
  dpch <- c(raw=1, residual=3)
  dlty <- c(raw=1, residual=3)
  if (missing(col)) col <- dcol[terms]
  if (is.null(names(col))) names(col) <- terms
  if (missing(type)) type <- dtype[terms]
  if (is.null(names(type))) names(type) <- terms
  type[type == "lines"] <- "l"; type[type == "points"] <- "p"
  if (missing(pch)) {pch <- ifelse(type == "p", dpch, NA)
  } else { pch <- ifelse(type == "p", pch, NA) }
  if (is.null(names(pch))) names(pch) <- terms
  if (missing(lty)) {lty <- ifelse(type == "l", dlty, NA)
  } else { lty <- ifelse(type == "l", lty, NA) }
  if (is.null(names(lty))) names(lty) <- terms
  if (all(is.na(pch))) pch <- NULL
  if (all(is.na(lty))) lty <- NULL
  if (length(terms) > 1 & missing(legend)) {legend <- TRUE
  } else if (length(terms) == 1 & missing(legend)) {legend <- FALSE}
  
  # Distance Matrix
  distance<-as.numeric(dist(coords))
  # Raw
  if ("raw" %in% terms) {
    raw <- sv(x = y, distance = distance, bins = bins, ...)
  }
  # Residual
  if ("residual" %in% terms) {
    residual <- sv(x = x$residuals, distance = distance, bins = bins, ...)
  }
  
  for (i in terms) {
    if (i == terms[1]) {
      sv.plot(get(i),distance=distance,bins=bins,add=FALSE,type=type[i],col=col[i],
              pch=pch[i],lty=lty[i])
    } else {
      sv.plot(get(i),distance=distance,bins=bins,type=type[i],col=col[i],pch=pch[i],
              lty=lty[i],add=TRUE)
    }
  }
  
  if (legend == TRUE) {
    lnames <- c(raw = "Raw data", residual = "Residuals")
    legend("top", legend=lnames[terms], pch = as.vector(pch[terms]), inset=c(0,-.15), xpd=TRUE, 
           bty = "n",col=col[terms], lty=as.vector(lty[terms]), cex=0.8, ncol=length(terms)) }
}

#' @rdname semivariogram
#' @export
semivariogram.default <- function(x, ..., coords, data, bins=13, east, north, type, 
                                       pch, lty, col){
  if (missing(coords)) {
    coords <- cbind(east, north)
    if (is.character(coords) & length(coords) == 2) coords <- as.vector(coords)
  }
  if (is.character(coords) & length(coords) == 2) coords <- data[coords]
  if (is.character(x) & length(x) == 1) x <- data[x]
  
  if (missing(type)) type <- "p"
  type[type == "points"] <- "p"; type[type == "lines"] <- "l"
  if (type == "p") pch <- ifelse(missing(pch), 1, pch); lty <- NULL
  if (type == "l") lty <- ifelse(missing(lty), 1, lty); pch <- NULL
  if (missing(col)) col <- "black"
  distance <- as.numeric(dist(coords))
  semivariances <- sv(x = x, distance = distance, bins = bins, ...)
  sv.plot(semivariances,distance=distance,bins=bins,add=FALSE,type=type,col=col,
          pch=pch,lty=lty)
}

#' @rdname semivariogram
#' @export
semivariogram.semivariance <- function(x, ..., type, pch, lty, legend, col){
  if (!inherits(x, "semivariance")) stop("The input x is not a 'semivariance' x.")
  if (is.list(x)){
    terms <- names(x)
    dcol <- c(raw="black", residual="blue", parametric="red")
    dtype <- c(raw="p", residual="p", parametric="l")
    dpch <- c(raw=1, residual=3, parametric=2)
    dlty <- c(raw=3, residual=2, parametric=1)
    if (missing(col)) col <- dcol[terms]
    if (is.null(names(col))) names(col) <- terms
    if (missing(type)) type <- dtype[terms]
    if (is.null(names(type))) names(type) <- terms
    type[type == "lines"] <- "l"; type[type == "points"] <- "p"
    if (missing(col)) col <- dcol[terms]
    if (is.null(names(col))) names(col) <- terms
    if (missing(type)) type <- dtype[terms]
    if (is.null(names(type))) names(type) <- terms
    if (missing(pch)) {pch <- ifelse(type == "p", dpch, NA)
    } else { pch <- ifelse(type == "p", pch, NA) }
    if (is.null(names(pch))) names(pch) <- terms
    if (missing(lty)) {lty <- ifelse(type == "l", dlty, NA)
    } else { lty <- ifelse(type == "l", lty, NA) }
    if (is.null(names(lty))) names(lty) <- terms
    if (all(is.na(pch))) pch <- NULL
    if (all(is.na(lty))) lty <- NULL
    if (length(terms) > 1 & missing(legend)) {legend <- TRUE
    } else if (length(terms) == 1 & missing(legend)) {legend <- FALSE}
    
    for (i in terms) {
      if (i == terms[1]) {
        distance <- as.numeric(names(x[[i]]))
        sv.plot(x[[i]],distance=distance,bins=length(distance),add=FALSE,type=type[i],col=col[i],
                pch=pch[i],lty=lty[i])
      } else {
        sv.plot(x[[i]],distance=distance,bins=length(distance),type=type[i],col=col[i],pch=pch[i],
                lty=lty[i],add=TRUE)
      }
      if (legend == TRUE) {
        legend("top", legend=terms, pch = as.vector(pch[terms]), inset=c(0,-.15), xpd=TRUE, 
               bty = "n",col=col[terms], lty=as.vector(lty[terms]), cex=0.8, ncol=length(terms)) }
    }
  } else if (is.matrix(x)) {
    if (length(x) > 100) stop("The length of distance is too large.")
    i <- 1
    if (missing(type)) type <- ifelse(length(x) > 30, "l", "p")
    type[type == "points"] <- "p"; type[type == "lines"] <- "l"
    if (type == "p") pch <- ifelse(missing(pch), 1, pch); lty <- NA
    if (type == "l") lty <- ifelse(missing(lty), 1, lty); pch <- NA
    if (missing(col)) col <- "black"
    distance <- as.numeric(names(x))
    sv.plot(x[[i]],distance=distance,bins=length(distance),add=FALSE,type=type[i],col=col[i],
            pch=pch[i],lty=lty[i])
  }
}

#' @rdname semivariogram
#' @export
plot.semivariance <- semivariogram.semivariance

sv.plot <- function(semivariances,distance,bins,add=FALSE,type="p",col=col,pch=1,lty=1,...){
  breaks<-seq(0,max(distance),length=bins+1)
  digits<-ifelse(max(breaks)<10,2,0)
  labels<-round(breaks[-1],digits)
  cl <- match.call(expand.dots=TRUE)
  if (! "xlab" %in% names(cl)) xlab <- "Distance"
  if (! "ylab" %in% names(cl)) ylab <- "Semivariance"
  if (add==FALSE){
    plot(1, type="n", xlab=xlab, ylab=ylab, xlim=c(0,bins), 
         ylim=c(0,max(semivariances)),axes=F,...)
    axis(1,at=1:bins,labels=labels,cex.axis=.75);axis(2);box()
  }
  if (type=="p"){
    points(semivariances, col = col, pch = pch)
  }
  if (type=="l"){
    lines(semivariances, col = col, lty = lty)
  }
}