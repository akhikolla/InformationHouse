#' Dimension reduction
#'
#' This function provides several methods to estimate the central space or central mean space of y on x.
#'  It returns the matrix of central space or central mean space for different dimensions and contains
#'  other information used for dimension selection by \code{\link{mave.dim}}.
#' @param formula the model used in regression
#' @param data the data
#' @param x The n by p design matrix.
#' @param y The n by q respond matrix.
#' @param method This parameter specify which method will be used in dimension reduction. It provides
#' five methods, including "csMAVE","csOPG","meanOPG","meanMAVE","KSIR"
#' by default, method = 'csOPG'
#' \itemize{
#' \item 'meanOPG' and 'meanMAVE' estimate dimension reduction space
#'            for conditional mean
#' \item 'csMAVE' and 'csOPG' estimate the central dimension reduction
#'            space
#' \item 'KSIR' is a kernel version of sliced inverse regression (Li, 1991). It is fast, but
#'            with poor accuracy.
#' }
#' @param max.dim the maximum dimension of dimension reduction space. The default is 10. In practice, max.dim
#' will be equal to min(max.dim,ncol(x),screen).
#' @param screen specify the number of variables retained after screening method. The default is n/log(n).
#' When this number is smaller than max.dim, then max.dim will change to the value of screen
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default
#' is na.action, which wil stop calculations. If na.action is set to be na.omit, the incomplete cases
#' will be removed.
#'
#' @return dr is a list which contains:
#' \itemize{
#' \item dir: dir[[d]] is the central space with d-dimension
#'         d = 1, 2, ..., p reduced direction of different dimensions
#' \item y: the value of response
#' \item idx: the index of variables which survives after screening
#' \item max.dim: the largest dimensions of CS or CMS which have been calculated in mave function
#' \item ky: parameter used for DIM for selection
#' \item x: the original training data
#' }
#'
#' @export
#'
#' @seealso \code{\link{mave.dim}} for dimension selection, \code{\link{predict.mave}} for prediction
#' using the dimension reduction space, \code{\link{coef.mave}} for accessing the basis vectors of
#' dimension reduction space of given dimension, \code{\link{plot.mave}} for plot method for mave class
#'
#' @references Li K C. Sliced inverse regression for dimension reduction[J]. Journal of the American Statistical Association, 1991, 86(414): 316-327.
#' @references Xia Y, Tong H, Li W K, et al. An adaptive estimation of dimension reduction space[J]. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 2002, 64(3): 363-410.
#' @references Xia Y. A constructive approach to the estimation of dimension reduction directions[J]. The Annals of Statistics, 2007: 2654-2690.
#' @references Wang H, Xia Y. Sliced regression for dimension reduction[J]. Journal of the American Statistical Association, 2008, 103(482): 811-821.
#'
#' @examples
#'  x <- matrix(rnorm(400*5),400,5)
#'  b1 <- matrix(c(1,1,0,0,0),5,1)
#'  b2 <- matrix(c(0,0,1,1,0),5,1)
#'  eps <- matrix(rnorm(400),400,1)
#'  y <- x%*%b1 + (x%*%b2)*eps
#'
#'  #finding central space based on OPG method
#'  #dr.csopg <- mave.compute(x,y, method = 'csopg')
#'  #or
#'  dr.csopg <- mave(y ~ x, method = 'csopg')
#'
#   #find central mean space based on MAVE method
#'  #dr.meanopg <- mave.compute(x,y, method = 'meanopg')
#'  #or
#'  dr.meanopg <- mave(y ~ x, method = 'meanopg')
#'
#'  #find central mean space based on ksir method
#'  dr.ksir <- mave(y~x,method='ksir')
#'  #or
#'  #dr.ksir <- mave.compute(x,y,method='ksir')
#'
#'  #See more examples about screening and mutiple responses in the vignette
#'  #Using screening for high dimensional data
#'  #x <- matrix(rnorm(100*50),100,50)
#'  #y1 = as.matrix(x[,1])+rnorm(100)*.2
#'  #y2 = as.matrix(x[,2]+x[,3])*as.matrix(x[,1]+x[,5])+rnorm(100)*.2
#'  #y = cbind(y1,y2)
#'  #dr.sc = mave(y~x,method='CSOPG',max.dim=5,screen=20)
#'  #dr.sc.dim = mave.dim(dr.sc)
#'  #print the directions of central space with the selected variables
#'  #dr.sc.dim$dir[[3]][dr.sc$idx,]
#'
#' @useDynLib MAVE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics pairs plot points
#' @importFrom stats model.extract model.matrix na.fail sd predict
#' @importFrom mda mars


mave<-function(formula, data, method='CSOPG', max.dim=10, screen=NULL, subset, na.action=na.fail){

  call <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- as.matrix(model.extract(mf, "response"))
  X <- model.matrix(mt, mf)
  y.name <- if (is.matrix(Y))
    colnames(Y)
  else as.character(attr(mt, "variables")[2])
  int <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  if (int > 0)
    X <- X[, -int, drop = FALSE]
  if(is.null(screen)){
    screen = ceiling(nrow(X)/log(nrow(X)))
  }
  dr <- mave.compute(X, Y, method, max.dim = max.dim,screen = screen)

  dr$call=call
  return(dr)
}

#' @rdname mave
#' @export
mave.compute<-function(x, y, method='CSOPG', max.dim = 10, screen=nrow(x)/log(nrow(x))){

  method=toupper(method);
  methodvec=c('CSOPG','CSMAVE','MEANOPG','MEANMAVE','KSIR')

  if(!(method %in% methodvec)){
    stop('method should be one of CSMAVE, CSOPG, MEANOPG, MEANMAVE, KSIR')
  }
  y = as.matrix(y)
  if(nrow(x)!=nrow(y)){
    stop('the row of x and the row of y is not compatible')
  }
  if(length(x)==length(y)){
    stop('x is one dimensional, no need do dimension reduction')
  }
  if(screen<ncol(x)){
    cat("screening method is using to select import variables.")
  }

  screen <- min(ncol(x),ceiling(screen))
  max.dim <- min(max.dim,ncol(x),screen)

  x.scaled <- scale(x)  #scale x
  dr <- MAVEfastCpp(x.scaled,y,method,max.dim,screen)

  dr$x <- x
  dr$call <- match.call()
  dr$method <- method
  dr$dir <- M3d2list(dr$dir,colnames(x))
  dr$dir <- dr$dir[1:max.dim]
  len = apply(x,2,sd)
  for(i in 1:max.dim){
    dir <- dr$dir[[i]]
    dir <- dir*matrix(1/len,ncol(x),i) #since x is scaled so we need to transform dir
    lendir = sqrt(apply(dir^2,2,sum))
    dr$dir[[i]] <- dir/matrix(lendir,ncol(x),i,byrow=T)
  }
  dr$max.dim = max.dim
  dr$y=y;
  class(dr) <- 'mave'
  return(dr)
}
