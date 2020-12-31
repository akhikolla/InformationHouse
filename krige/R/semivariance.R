############################ SEMIVARIANCE FUNCTION #############################
# LAST UPDATE: 11/15/2020

#' Semivariance for Geospatial Data
#'
#' This function computes the empirical semivariance for a spatially-distributed 
#' variable. Based on the user's chosen level of coarsening, the semivariance is 
#' presented for various distances.
#'
#' @param object An object for which the semivariance is desired. The object can 
#'   be a \code{krige} object, a \code{lm} object, or a vector of variables (or 
#'   variable names in the \code{data}).
#' @param coords A matrix of coordinates for all observations or a vector of variable 
#'   names indicating the coordinates variables in the data. Alternatively, the 
#'   coordinates can also be specified separately using \code{east} and \code{north}.
#' @param data If object is a variable name, a data frame must be provided.
#' @param bins Number of bins into which distances should be divided. The observed 
#'   distances will be split at equal intervals, and the semivariance will be computed 
#'   within each interval. Defaults to 13 intervals.
#' @param terms A vector of strings specifies for which the semivariogram is created. 
#'   Options are "raw" (the semivariogram for raw data), "residual" (the semivariogram
#'   for residuals from linear regression).
#' @param east Alternative specification for the vector of eastings for all observations.
#' @param north Alternative specification for the vector of northing for all observations.
#' @param plot Logical values indicates whether a graph of the empirical semivariogram 
#'   should be presented with a run of the function. Default omits the plot and only 
#'   returns semivariance values. See \code{\link{semivariogram}} for additional 
#'   plotting functions.
#' @param \dots Additional arguments passed to \code{semivariance} methods. 
#' 
#' @return A semivariance object. It will be a numeric vector with each bin's value 
#'   of the semivariance if only one kind of semivariance is computed; a list including 
#'   different kinds of semivariance if both raw and residual semivariance is computed.
#'   
#' @details Semivariance is equal to half of the variance of the difference in a 
#'   variable's values at a given distance. That is, the semivariance is defined 
#'   as: \eqn{\gamma(h)=0.5*E[X(s+h)-X(s)]^2}, where \eqn{X} is the variable of 
#'   interest, s is a location, and h is the distance from s to another location.
#'   
#'   The function can be applied to a variable, a fitted linear model (\code{lm} 
#'   object) before fitting a spatial model or to a \code{krige} object or \code{semivariance}
#'   object to assess the model fit. When applying to a variable, it will describes 
#'   the raw data; for a \code{lm} object, the default will present empirical 
#'   semivariogram for both the raw data and linear residuals. Users can also specify 
#'   which semivariance is needed in the \code{terms} argument if there are multiple 
#'   kinds of semivariogram can be plotted. A \code{semivariance} object can also 
#'   be used to create semivariogram afterwards using generic \code{plot} function 
#'   with more options.
#'   
#' @references 
#'   Sudipto Banerjee, Bradley P. Carlin, and Alan E. Gelfand. 2015. \emph{Hierarchical 
#'   Modeling and Analysis for Spatial Data}. 2nd ed. Boca Raton, FL: CRC Press.
#'   
#' @seealso \code{\link{semivariogram}}, \code{\link{plot.semivariance}}, \code{\link{exponential.semivariance}}   
#' 
#' @examples
#' \dontrun{
#' # Summarize example data
#' summary(ContrivedData)
#' 
#' # Empirical semivariance for variable y
#' semivariance(ContrivedData$y,coords = cbind(ContrivedData$s.1, ContrivedData$s.2))
#' 
#' # Initial OLS Model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData); summary(contrived.ols)
#' 
#' # Empirical semivariance for ols fit
#' (sv.ols <- semivariance(contrived.ols, coords = c("s.1","s.2"), bins=13))
#' plot(sv.ols)
#' 
#' # Estimation using metropolis.krige()
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'   data = ContrivedData, n.iter = M, range.tol = 0.05)
#'   
#' # Parametric semivariance
#' (sv.krige <- semivariance(contrived.run, plot = TRUE))
#' 
#' # Convert to other format for further use
#' as.matrix(sv.krige)
#' as.data.frame(sv.krige)
#' }

#' 
#' @importFrom stats dist as.formula model.frame
#' @importFrom graphics plot axis lines legend points
#' @export

semivariance<-function(object, ...) {
  UseMethod("semivariance")
}

#' @rdname semivariance
#' @export
semivariance.krige <- function(object, bins=13, terms = "all", plot = FALSE, ...){
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object.")
  if (terms == "all") terms <- c("raw", "residual", "parametric")
  
  # Distance Matrix
  distance<-as.numeric(dist(cbind(object$model.data.list$easting,
                                  object$model.data.list$northing)))
  
  # Raw
  if ("raw" %in% terms) {
    raw <- sv(x = object$model.data.list$y, distance = distance, bins = bins, ...)
    class(raw) <- c("semivariance", class(raw))
  }
  
  # Residual
  if ("residual" %in% terms) {
    residual <- sv(x = object$init.ols$residuals, distance = distance, bins = bins, ...)
    class(residual) <- c("semivariance", class(residual))
  }
  
  # Parametric exponential semivariance from estimation
  if ("parametric" %in% terms) {
    if (!exists("resid")) {residual <- sv(x = object$init.ols$residuals, 
                                       distance = distance)}
    var.terms<-apply(object$mcmc.mat[,1:3],2,quantile,0.5)
    parametric <- exponential.semivariance(nugget=var.terms[1], decay=var.terms[2], 
                                            partial.sill = var.terms[3], 
                                            distance = as.numeric(names(residual)), 
                                            power = object$standing.parameter$powered.exp)
    class(parametric) <- c("semivariance", class(residual))
  }
  if (plot == TRUE) {
    semivariogram(x = object, terms = terms, bins=bins)
  }
  if (length(terms) == 1) {
    semivariance <- get(terms)
    class(semivariance) <- c("semivariance", "numeric")
  } else {
    semivariance <- list()
    for (i in terms) {
      semivariance[[i]] <- get(i)
    }
    class(semivariance) <- c("semivariance", "list")
  }
  semivariance
}

#' @rdname semivariance
#' @export
semivariance.lm <- function(object, bins=13, coords, terms = c("raw", "residual"), 
                             east, north, plot = FALSE, ...) {
  if (missing(coords)) {
    if (missing(east) | missing(north)) stop("Coordinates are missing.")
    coords <- cbind(east, north)
    if (is.character(coords) & length(coords) == 2) coords <- as.vector(coords)
  }
  cl <- object$call
  y <- object$model[1]
  if (is.character(coords) & length(coords) == 2) {
    form <- update(as.formula(cl$formula), paste("~ . +",paste(coords, collapse=" + ")))
    na.action <- ifelse("na.action" %in% names(cl), cl$na.action, "na.omit")
    mf <- model.frame(form, data = get(as.character(cl$data)), na.action=na.action)
    y <- mf[1]; coords <- mf[coords]
  }
  # Distance Matrix
  distance<-as.numeric(dist(coords))
  # Raw
  if ("raw" %in% terms) {
    raw <- sv(x = y, distance = distance, bins = bins)
    class(raw) <- c("semivariance",  class(raw))
  }
  # Residual
  if ("residual" %in% terms) {
    residual <- sv(x = object$residuals, distance = distance, bins = bins)
    class(residual) <- c("semivariance", class(residual))
  }
  if (plot == TRUE) semivariogram.lm(object, coords, bins, terms)
  if (length(terms) == 1) {
    semivariance <- get(terms)
    class(semivariance) <- c("semivariance", "numeric")
  } else {
    semivariance <- list()
    for (i in terms) {
      semivariance[[i]] <- get(i)
    }
    class(semivariance) <- c("semivariance", "list")
  }
  semivariance
}

#' @rdname semivariance
#' @export
semivariance.default <- function(object, bins=13, coords, data, east, north, plot=FALSE, ...){
  if (missing(coords)) {
    coords <- cbind(east, north)
    if (is.character(coords) & length(coords) == 2) coords <- as.vector(coords)
  }
  if (is.character(coords) & length(coords) == 2) coords <- data[coords]
  if (is.character(object) & length(object) == 1) object <- data[object]
  semivariance <- sv(x = object, distance = as.numeric(dist(coords)), bins = bins)
  if (plot == TRUE) semivariogram.default(object, coords, bins)
  class(semivariance) <- c("semivariance", class(semivariance))
  semivariance
}

################ PARAMETRIC POWERED EXPONENTIAL SEMIVARIANCE ###################
# LAST UPDATE: 11/15/2020

#' Parametric Exponential Semivariance
#'
#' This function returns the value of a parametric powered exponential semivariogram 
#' given the values of the parameters and the distance between observations.
#'
#' @param object A \code{krige} object of which the values of the estimates are used 
#'   to calculate the exponential semivariance.
#' @param nugget The value of the non-spatial variance, or nugget term.
#' @param decay The value of the decay term that sets the level of correlation given distance.
#' @param partial.sill The value of the spatial variance, or partial sill term.
#' @param distance The distance among observations for which the semivariance value is desired.
#' @param power The exponent specified in the powered exponential semivariogram. 
#'   Defaults to 2, which corresponds to a Gaussian semivariance function.
#' @param \dots Additional arguments
#' 
#' @return A semivariance object. It will be a numeric vector with each bin's value 
#'   of the semivariance.
#'   
#' @details The models estimated by the \code{krige} package assume a powered exponential 
#'   covariance structure. Each parametric covariance function for kriging models 
#'   corresponds to a related semivariance function, given that highly correlated 
#'   values will have a small variance in differences while uncorrelated values 
#'   will vary widely. More specifically, semivariance is equal to half of the 
#'   variance of the difference in a variable's values at a given distance. That is, 
#'   the semivariance is defined as: \eqn{\gamma(h)=0.5*E[X(s+h)-X(s)]^2}, where \eqn{X} 
#'   is the variable of interest, s is a location, and h is the distance from s 
#'   to another location. 
#'   
#'   The powered exponential covariance structure implies that the semivariance 
#'   follows the specific functional form of \eqn{\gamma(d)=\tau^2+\sigma^2(1-\exp(-|\phi d|^p))} 
#'   (Banerjee, Carlin, and Gelfand 2015, 27). A perk of this structure is that 
#'   the special case of \emph{p=1} implies the commonly-used exponential semivariogram, 
#'   and the special case of \emph{p=2} implies the commonly-used Gaussian semivariogram. 
#'   Upon estimating a model, it is advisable to graph the functional form of the 
#'   implied parametric semivariance structure. By substituting estimated values 
#'   of the \code{nugget}, \code{decay}, and \code{partial.sill} terms, as well 
#'   as specifying the correct \code{power} argument, it is possible to compute 
#'   the implied semivariance from the model. The \code{distance} argument easily 
#'   can be a vector of observed distance values.
#'   
#' @references 
#'   Sudipto Banerjee, Bradley P. Carlin, and Alan E. Gelfand. 2015. \emph{Hierarchical 
#'   Modeling and Analysis for Spatial Data}. 2nd ed. Boca Raton, FL: CRC Press.
#'   
#' @seealso \code{\link{semivariogram}}, \code{\link{plot.semivariance}}, \code{\link{exponential.semivariance}}   
#' 
#' @examples
#' \dontrun{
#' # Summarize data
#' summary(ContrivedData)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'   data = ContrivedData, n.iter = M, range.tol = 0.05)
#'   
#' # Parametric powered exponential semivariogram
#' exponential.semivariance(contrived.run)
#' 

#' #OLS Model for Residuals
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' 
#' # Residual semivariance
#' (resid.semivar <- semivariance(contrived.ols, coords = c("s.1", "s.2"), terms = "residual"))
#' 
#' # Parametric exponential semivariance
#' exponential.semivariance(nugget=0.5,decay=2.5,partial.sill=0.5, 
#'                          distance=as.numeric(names(resid.semivar)))
#' }
#' 
#' @importFrom stats dist
#' @export

exponential.semivariance <- function(...){
  UseMethod("exponential.semivariance")
}

#' @rdname exponential.semivariance
#' 
#' @importFrom stats dist
#' @export
exponential.semivariance.krige <- function(object,...){
  var <- apply(object$mcmc.mat[,1:3], 2, mean)
  disMat <- as.numeric(dist(cbind(object$model.data.list$easting,
                                    object$model.data.list$northing)))
  semivariance <- exponential.semivariance.default(nugget=var[1],
                                                   decay=var[2],
                                                   partial.sill=var[3],
                                                   distance=disMat,
                                                   power=object$standing.parameter$powered.exp)
  semivariance
}

#' @rdname exponential.semivariance
#' @export
exponential.semivariance.default <- function(nugget,decay,partial.sill,distance,power=2,...){
  semivariance<-nugget+partial.sill*(1-exp(-abs(decay*distance)^power))
  names(semivariance) <- distance
  class(semivariance) <- "semivariance"
  semivariance
}

#' Convert semivariance to a matrix object
#' 
#' @param x An \code{semivariance} object created by \code{semivariance} or 
#'   \code{exponential.semivariance}.
#' @param \dots Additional arguments passed to \code{as.matrix} and \code{as.data.frame}
#'   methods. Not supported for \code{semivariance} object.
#'   
#' @details The defaults of semivariance methods give a list or numeric vector. These
#'   methods can convert the semivariance output list and vector to \code{matrix} or 
#'   \code{data.frame}.
#'   
#' @return A matrix containing the computed distance and semivariance.
#' 
#' @examples
#' \dontrun{
#' # Summarize Data
#' summary(ContrivedData)
#' 
#' # Empirical semivariance for variable y
#' raw.var <- semivariance(x=ContrivedData$y,coords = cbind(ContrivedData$s.1, 
#'   ContrivedData$s.2))
#' as.matrix(raw.var)
#' 
#' # Estimation using metropolis.krige()
#' #' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, range.tol = 0.05)
#' 
#' # Parametric semivariance
#' parametric.var <- semivariance(contrived.run)
#' as.matrix(parametric.var)
#' as.data.frame(parametric.var)
#' }
#' 
#' @export
as.matrix.semivariance <- function(x, ...) {
  if (!inherits(x, "semivariance")) stop("The input x is not a 'semivariance' x.")
  if (inherits(x, "list")){ 
    mat <- matrix(unlist(x),nrow = length(x[[1]]), ncol = length(x), byrow = FALSE)
    mat <- cbind(as.numeric(names(x[[1]])), mat)
    colnames(mat) <- c("breaks", names(x))
  } else if (inherits(x, "numeric")) {
    mat <- cbind(as.numeric(names(x)), as.vector(x))
    colnames(mat) <- c("break", "var")
  }
  #rownames(mat) <- NULL
  mat
}

#' @rdname as.matrix.semivariance
#' @export
as.data.frame.semivariance <- function(x, ...) {
  if (!inherits(x, "semivariance")) stop("The input object is not a 'semivariance' object.")
  as.data.frame(as.matrix.semivariance(x))
}

#semivariance function from Banerjee, Carlin, and Gelfand (2015, p. 24, Eq. 2.1)
sv<-function(x,distance,east,north,bins){
  if (missing(distance)) {
    distance<-as.numeric(dist(cbind(east,north)))
  } else {
    distance <- distance
  }
  differences<-as.numeric(dist(x))
  breaks<-seq(0,max(distance),length=bins+1)
  categories<-cut(distance,breaks=breaks)
  semivar.calc<-function(x){0.5*mean(x^2)}
  semivariances<-as.numeric(by(differences,INDICES=categories,FUN=semivar.calc))
  names(semivariances)<-breaks[-1]
  return(semivariances)
}
