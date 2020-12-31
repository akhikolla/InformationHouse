########################## PREDICTION FUNCTION ################################
## LAST UPDATE: 11/08/2020; Le Bao

#' Predictions by Kriging
#'
#' This function uses the results of a model estimated by \code{metropolis.krige} 
#'   to make kriging-based predictions.
#'
#' @param object An \code{krige} object estimated by \code{metropolis.krige}.
#' @param newdata An optional data frame in which to look for variables with which 
#'   to predict. If omitted, the fitted values are produced. Alternatively, the new 
#'   data can be specified using \code{new.X}, \code{new.east}, and \code{new.north}.
#' @param credible If a credible interval on predictions is desired, a user may 
#'   specify a proportion between 0 and 1 to indicate the interval probability. 
#'   For example, a value of 0.9 would create a 90\% credible interval. If \code{NULL}, 
#'   then no credible interval will be produced.
#' @param new.X The matrix of independent variables for observations to be predicted.
#' @param new.east Vector of eastings for observations to be predicted.
#' @param new.north Vector of northings for observations to be predicted.
#' @param \dots Additional arguments passed to \code{predict} methods. Not supported 
#'   for \code{krige} objects.
#' 
#' @return An object of class \code{matrix} with one prediction per row. By default 
#'  the matrix has one column, as only point predictions are returned. If the \code{credible} 
#'  option is specified, there are three columns respectively indicating a point 
#'  estimate (median prediction from MCMC), lower bound of the credible interval, 
#'  and upper bound of the credible interval.
#'   
#' @details Analysts should use this function if they want to make kriged predictions 
#'    for observations at new locations or predict fitted values at original locations. 
#'    To do this, researchers first must estimate a model using the \code{metropolis.krige} 
#'    function.
#'  
#'  After estimating the model, a \code{krige} object can provide information from 
#'    the model fitted with \code{metropolis.krige}, including the original imput 
#'    data, coordinates, and the model results themselves. The prediction will also 
#'    use the same \code{powered.exp}. \code{new.data} can be specified for predicting 
#'    the observations at new location. Otherwise, the fitted values for the original 
#'    data and locations will be produced.
#'  
#'  By default, the function uses median values of parameters to make a single point 
#'   prediction for every kriged data point. However, if the uses specifies a probability 
#'   with the \code{credible} option, then the function will determine the predictions 
#'   for all iterations of the MCMC sample. The point estimates will then be a median 
#'   of these predictions, and a credible interval will be returned based on percentiles. 
#'   Note that estimating a credible interval is substantially more intensive computationally, 
#'   but has the benefit of reporting uncertainty in predictions.
#'   
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'   \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'   
#' @examples
#' \dontrun{
#' # Summarize Data
#' summary(ContrivedData)
#' 
#' # Initial OLS model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' # summary(contrived.ols)
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
#' # Predict fitted values
#' predict(contrived.run)
#' 
#' # Predict new data
#' euler<-c(0.2,0.7)
#' archimedes<-c(0.3,0.1)
#' pythagoras<-c(0.1,0.4)
#' mathematicians<-rbind(euler,archimedes,pythagoras)
#' basel<-c(0.1,0.8)
#' sicily<-c(0.4,0.1)
#' samos<-c(0.1,0.4)
#' new.locations<-rbind(basel,sicily,samos)
#' newDf <- as.data.frame(cbind(mathematicians, new.locations))
#' colnames(newDf) <- c("x.1", "x.2", "s.1", "s.2")
#' 
#' # Make predictions from median parameter values:
#' (median.pred <- predict(contrived.run, newdata = newDf))
#' 
#' # Make predictions with 90\% credible intervals:
#' (cred.pred <- predict(contrived.run, newdata = newDf, credible=0.9))
#' }
#' 
#' @importFrom stats median quantile 
#' @importFrom Rcpp evalCpp
#' @export
predict.krige <- function(object, newdata, credible=FALSE, new.X, new.east, 
                          new.north, ...){
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object")
  cl <- match.call()
  powered.exp <- object$standing.parameter$powered.exp
    
  # TRAINING DATA
  train.y <- object$model.data.list$y
  train.X <- object$model.data.list$X
  train.east <- object$model.data.list$easting
  train.north <- object$model.data.list$northing
  mcmc.iter<- object$mcmc.mat
  n.iter <- object$n.iter - object$n.burnin
  
  # PREDICT DATA
  if ("newdata" %in% names(cl)) {
    pred.X <- model.matrix(object$formula[-2], data = newdata)
    pred.east <- as.vector(unlist(newdata[object$coords[1]]))
    pred.north <- as.vector(unlist(newdata[object$coords[2]]))
  } else if (sum(c("new.X", "new.east", "new.north") %in% names(cl))==3){
    pred.X <- new.X
    pred.east <- new.east
    pred.north <- new.north
  } else {
    pred.X <- object$model.data.list$X
    pred.east <- object$model.data.list$easting
    pred.north <- object$model.data.list$northing
  }
  
  # DISTANCE MATRIX
  train.coord<-cbind(train.east,train.north)
  test.coord<-cbind(pred.east,pred.north)
  all.dist<-k_distmat(rbind(test.coord,train.coord))
  train.dist<-all.dist[-c(1:nrow(test.coord)),-c(1:nrow(test.coord))]
  pred.dist<-all.dist[1:nrow(test.coord),-c(1:nrow(test.coord))]
  
  if(credible == FALSE) {
    beta<-apply(mcmc.iter[,-c(1:3)],2,median)
    tau2<-median(mcmc.iter[,1])
    phi<-median(mcmc.iter[,2])
    sigma2<-median(mcmc.iter[,3])

    Sigma<-ifelse(train.dist>0, sigma2*exp(-abs(phi*train.dist)^powered.exp), tau2+sigma2)
    inv.Sigma<-k_inv(Sigma)
    smoother<-sigma2*exp(-abs(phi*pred.dist)^powered.exp)
    resid<-train.y-train.X%*%beta
    forecast<-as.vector(pred.X%*%beta+smoother%*%inv.Sigma%*%resid)
    names(forecast)<-rownames(pred.X)
  } else {
    no.pred<-ifelse(is.matrix(pred.X),nrow(pred.X),1)
    forecast.iter<-matrix(NA,nrow=no.pred,ncol=nrow(mcmc.iter))
    
    for(i in 1:n.iter){
      beta<-mcmc.iter[i,-c(1:3)]
      tau2<-mcmc.iter[i,1]
      phi<-mcmc.iter[i,2]
      sigma2<-mcmc.iter[i,3]
      Sigma<-ifelse(train.dist>0, sigma2*exp(-abs(phi*train.dist)^powered.exp), tau2+sigma2)
      smoother<-sigma2*exp(-abs(phi*pred.dist)^powered.exp)
      resid<-train.y-train.X%*%beta
      forecast.iter[,i]<-pred.X%*%beta+smoother%*%k_inv(Sigma)%*%resid
    }
    tail <- (1-credible)/2
    forecast<-t(apply(forecast.iter,1,quantile,c(.5,tail,1-tail)))
    colnames(forecast)<-c("Estimate","Lower C.I.","Upper C.I.")
    rownames(forecast)<-rownames(pred.X)
  }
  return(forecast)
}
