#' @name Fitters
#' @rdname Fitters
#' @aliases Fitters UnivariateFitter GenUnivariateFitter
#' @title Functions used to fit GeDS  objects.
#' @param X a numeric vector containing   \eqn{N} sample values of the covariate chosen to enter the spline
#' regression component of the predictor model.
#' @param Y a vector of size \eqn{N} containing
#' the observed values of the response variable \eqn{y}.
#' @param Z a design matrix with \eqn{N} rows containing
#' other covariates selected to enter the parametric component of the predictor model
#' (see \code{\link[=formula.GeDS]{formula}}). If no such covariates are selected, it is set to \code{NULL} by default.
#' @param family a description of the error distribution and link function to be used in the model. This can be a
#' character string naming a family function (e.g. \code{"gaussian"}),
#' the family function itself (e.g. \code{\link[stats]{gaussian}})
#' or the result of a call to a family function (e.g. \code{gaussian()}).
#'  See \link[stats]{family} for details on family functions.
#' @param weights an optional vector of size \eqn{N} of `prior weights' to be put
#' on the observations in the fitting
#' process in case the user requires weighted GeDS fitting.
#' It is \code{NULL} by default.
#' @param beta numeric parameter in the interval \eqn{[0,1]}
#' tuning the knot placement in stage A of GeDS. See the description of \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the threshold for
#'  the stopping rule  (model selector) in stage A of GeDS. See also \code{stoptype} and
#'  details in the description of \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param min.intknots optional parameter allowing the user to set
#' a minimum number of internal knots required. By default equal to zero.
#' @param max.intknots optional parameter allowing the user to set a maximum number
#'  of internal knots to be added by the GeDS estimation algorithm.
#' By default equal to the number of internal knots \eqn{\kappa} for
#' the saturated GeDS model (i.e. \eqn{\kappa=N-2}).
#' @param q numeric parameter which allows to fine-tune the stopping rule of stage A of GeDS, by default equal to 2.
#' See details in the description of \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param extr numeric vector of 2 elements representing the left-most and right-most limits
#' of the interval embedding the sample values of \code{X}.
#'  By default equal correspondingly to the smallest and largest values of \code{X}.
#' @param show.iters logical variable indicating whether or not to print
#' information at each step. By default equal to \code{FALSE}.
#' @param stoptype a character string indicating the type of GeDS stopping rule to
#' be used. It should be either \code{"SR"}, \code{"RD"} or
#' \code{"LR"}, partial match allowed. See details of \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param offset a vector of size \eqn{N} that can be used to specify a fixed covariate
#' to be included in the predictor model  avoiding the estimation of its corresponding regression coefficient.
#' In case  more than one covariate is fixed, the user should sum the corresponding coordinates of the fixed covariates
#'  to produce one common \eqn{N}-vector of coordinates.
#' The \code{offset} argument is particularly useful when using
#' \code{GenUnivariateFitter} if the link function used is not the identity.
#' @param tol numeric value indicating the tolerance to be used in the knot placement steps in stage A. By default equal to 1E-12. See details below.
#'
#' @description These are computing engines called by \code{\link{NGeDS}} and \code{\link{GGeDS}}, needed for
#' the underlying fitting procedures.
#'
#' @details The functions \code{UnivariateFitter} and \code{GenUnivariateFitter} are in general not intended to be used directly,
#' they should be called through \code{\link{NGeDS}} and \code{\link{GGeDS}}.
#' However, in case there is a need for multiple GeDS fitting (as may be the case e.g. in Monte Carlo simulations)
#' it may be efficient to use the fitters outside the main functions.
#'
#' The argument \code{tol} is used in the knot placement procedure of stage A of the GeDS algorithm in order to check
#' whether the current knot \eqn{\delta^*} is set at an acceptable location or not.
#' If there exists a knot \eqn{\delta_i}  such that \eqn{|\delta^* - \delta_i| < }\code{tol},
#' \eqn{\delta^*}, then the new knot is considered
#' to be coalescent with an existing one, it is discarded and the algorithm seeks alternative knot locations.
#' By default it is equal to 1e-12.
#'
#' See \code{\link{NGeDS}} and \code{\link{GGeDS}}, Kaishev et al. (2016) and Dimitrova et al. (2017) for further details.
#'
#' @seealso  \code{\link{NGeDS}} and \code{\link{GGeDS}}.
#'
#' @return A \code{\link{GeDS-Class}} object, but without the \code{Formula},
#' \code{extcall}, \code{terms} and \code{znames} slots.
#'
#'
#' @examples
#' # Examples similar to the ones
#' # presented in NGeDS and in GGeDS
#'
#' # Generate a data sample for the response variable
#' # Y and the covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' # Specify a model for the mean of Y to include only
#' # a component non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression model using the fitter function
#' (Gmod <- UnivariateFitter(X, Y, beta = 0.6, phi = 0.995,
#'            extr = c(-2,2)))
#'
#' ##############################################################
#' # second: very similar example, but based on Poisson data
#' set.seed(123)
#' X <- sort(runif(N , min = -2, max = 2))
#' means <- exp(f_1(X))
#' Y <- rpois(N,means)
#' (Gmod2 <- GenUnivariateFitter(X, Y, beta = 0.2,
#'             phi = 0.995, family = poisson(), extr = c(-2,2)))
#'
#' # a plot showing quadratic and cubic fits,
#' # in the predictor scale
#' plot(X,log(Y), xlab = "x", ylab = expression(f[1](x)))
#' lines(Gmod2, n = 3, col = "red")
#' lines(Gmod2, n = 4, col = "blue", lty = 2)
#' legend("topleft", c("Quadratic","Cubic"),
#'      col = c("red","blue"), lty = c(1,2))
#'
#'
#' @export
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S., & Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \href{https://doi.org/10.1007/s00180-015-0621-7}{doi.org/10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{http://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
UnivariateFitter <- function(X, Y, Z = NULL, offset = rep(0,NROW(Y)), weights = rep(1,length(X)), beta=.5, phi = 0.5,
                             min.intknots = 0,
                             max.intknots = 300, q = 2,
                             extr=range(X), show.iters=FALSE, tol = as.double(1e-12), stoptype = c("SR","RD","LR")) {
  save <- match.call()
  RSSnew <- numeric()
  stoptype <- match.arg(stoptype)
  phis <- NULL
  phis_star <- NULL
  oldintc <- NULL
  oldslp <- NULL
  Indicator <- table(X)
  indent <- rep(" ",nchar(options()$prompt))
  indent <- paste(indent,collapse="")

  args <- list("X" = X, "Y" = Y, "Z" = Z, "offset"=offset, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots, "q" = q,
               "extr" = extr, "tol" = tol)
  previous <- matrix(nrow=max.intknots+1, ncol=max.intknots+4)
  nz <- if(!is.null(Z)) NCOL(Z) else 0
  oldcoef <- matrix(nrow=max.intknots+1, ncol=max.intknots+2+nz)
  nodi <- NULL
  Xw <- unique(X)
  for(j in 1:min(max.intknots+1,length(Y)-2)){#
    #if(j >1)  nodi<-sort(nodi)
    first.deg<-SplineReg_fast_weighted_zed(X=X, Y=Y, Z=Z, weights=weights, offset=offset,
                                           extr=extr,InterKnots=nodi,n=2) #first regression
    previous[j,1:(j+3)] <- sort(c(nodi,rep(extr,2)))
    oldcoef[j,1:(j+1+nz)] <- first.deg$Theta
    res.tmp <- first.deg$Residuals
    RSS.tmp <- first.deg$RSS
    RSSnew <- c(RSSnew,RSS.tmp)
    prnt <- ""
    if(stoptype=="SR"){
      if(j>q  && length(phis)>=3){
        phis <- c(phis,RSSnew[j]/RSSnew[j-q])
        if(j-q > min.intknots){
          phismod <- log(1-phis)
          ccc<-.lm.fit(cbind(1,(q+1):j),phismod)$coef
          phis_star <- c(phis_star,1-exp(ccc[1])*exp(ccc[2]*j))
          oldintc <- c(oldintc,ccc[1])
          oldslp <- c(oldslp,ccc[2])
          prnt <- paste0(", phi_hat = ",round(1-exp(ccc[1])*exp(ccc[2]*j),3))

          if(1-exp(ccc[1])*exp(ccc[2]*j) >= phi) {
            break
          }
        }
      }
    }
    if(stoptype=="RD" | (stoptype=="SR" & length(phis)<3)){
      phis <- c(phis,RSSnew[j]/RSSnew[j-q])
      if(j>q && (j-q > min.intknots) ){
        prnt <- paste0(", phi = ",round(RSSnew[j]/RSSnew[j-q],3))
        if(RSSnew[j]/RSSnew[j-q] >= phi) {      # stop rule
          break
        }
      }
    }
    if(stoptype=="LR"){
      phis <- c(phis,RSSnew[j-q]-RSSnew[j])
      if(j>q && (j-q > min.intknots) ){
        prnt <- paste0(", p = ",round(pchisq(-(RSSnew[j]-RSSnew[j-q]),df=q),3))
        if(-(RSSnew[j]-RSSnew[j-q]) < qchisq(phi,df=q)) {      # stop rule
          break
        }
      }
    }
    d <- numeric()
    if(any(Indicator>1)){
      res.wkg <- makeNewRes(resold = res.tmp*weights, recurr = as.numeric(Indicator)) #already fast - useless to do in c++
    } else {
      res.wkg <- res.tmp*weights
    }
    segni <- sign(res.wkg)
    for(i in 1:length(Xw)){
      if(all(segni==segni[1])){
        d[i]<-length(segni)
        break
      } else {
        d[i]<-min(which(segni!=segni[1])-1)
        segni <- segni[-(1:d[i])]
      }
    }
    u<-length(d)
    dcum <- cumsum(d)
    wc.range<-numeric(u)
    means <- numeric(u)
    means[1] <- abs(mean(res.wkg[1:dcum[1]]))
    wc.range[1] <- Xw[dcum[1]]-X[1]
    for(i in 2:u){ # embed in C++ useless
      means[i] <- abs(mean((res.wkg[(dcum[i-1]+1):dcum[i]])))
      wc.range[i] <- Xw[dcum[i]]-Xw[dcum[i-1]+1]
    }
    means<-means/max(means)
    wc.range<-wc.range/max(wc.range)
    w <- beta*means + (1-beta)*wc.range
    newknot <- Knotnew(wht=w, restmp=res.wkg, x=Xw, dcm=dcum, oldknots=c(rep(extr,2+1),nodi),tol = tol)[1]# w <- salva <- salva +1
    nodi<-c(nodi,newknot)
    if(show.iters) {
      if(j>q) {
        toprint<-paste0(indent, "Iteration ",j,": New Knot = ", round(newknot,3)  ,
                        ", RSS = " , round(RSSnew[j],3), prnt,"\n")
      } else {
        toprint<-paste0(indent,"Iteration ",j,": New Knot = ", round(newknot,3)  ,", RSS = " , round(RSSnew[j],3), "\n")}
      cat(toprint)}
  }
  if (j == max.intknots + 1) warning("Maximum number of iterations exceeded")
  if (j<max.intknots) {
    previous <- previous[-((j+1):(max.intknots+1)),] #delete also last two$
    previous <- previous[,-((j+4):(max.intknots+4))]
    oldcoef <- oldcoef[-((j+1):(max.intknots+1)),]
    oldcoef <- oldcoef[,1:(j+1+nz)]
  }
  if(j-q<2){
    warning("Too few internal knots found: Linear spline will not be computed. Try to set a different value for 'q' or a different treshold")
    ll <- NULL
    lin <- NULL} else {
      ik <- previous[j-q,-c(1,2,(j+2-q):(j+3))]
      ll <- makenewknots(ik,2)#lin #l'approssimazione
      lin <- SplineReg_LM(X=X,Y=Y,Z=Z,offset=offset,weights=weights,extr=extr,InterKnots=ll,n=2)
    }
  if(j-q<3){
    warning("Too few internal knots found: Quadratic spline will not be computed. Try to set a different value for 'q' or a different treshold")
    qq <- NULL
    squ <- NULL} else {
      qq <- makenewknots(ik,3)#quad
      squ <- SplineReg_LM(X=X,Y=Y,Z=Z,offset=offset,weights=weights,extr=extr,InterKnots=qq,n=3)
    }
  if(j-q < 4){
    warning("Too few internal knots found: Cubic spline will not be computed. Try to set a different value for 'q' or a different treshold")
    cc <- NULL
    cub <- NULL} else {
      cc <- makenewknots(ik,4)#cub
      cub <- SplineReg_LM(X=X,Y=Y,Z=Z,offset=offset,weights=weights,extr=extr,InterKnots=cc,n=4)
    }
  out <- list("Type" = "LM - Univ","Linear.Knots"=ll,"Quadratic.Knots"=qq,"Cubic.Knots"=cc,"Dev.Linear" = lin$RSS,
              "Dev.Quadratic" = squ$RSS,"Dev.Cubic" = cub$RSS,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = previous,
              "Args"= args,"Call"= save, "Nintknots"= j-q-1,"iters" = j, "Guesses" = NULL,
              "Coefficients" = oldcoef, stopinfo = list("phis"=phis,"phis_star"=phis_star,"oldintc"=oldintc,"oldslp"=oldslp))

  class(out) <- "GeDS"
  return(out)
}
#' @rdname Fitters
#' @export
GenUnivariateFitter <- function(X, Y, Z = NULL, offset = rep(0,NROW(Y)),
                                weights=rep(1,length(X)), family=gaussian(), beta=.5, phi = 0.5,
                                min.intknots = 0, max.intknots = 300, q = 2, extr=range(X), show.iters=F, tol = as.double(1e-12),stoptype = c("SR","RD","LR")) {
  save <- match.call()
  stoptype <- match.arg(stoptype)
  RSSnew <- numeric()
  oldintc <- NULL
  oldslp<- NULL
  indent <- rep(" ",nchar(options()$prompt))
  indent <- paste(indent,collapse="")

  args <- list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "weights"=weights,
               "beta" = beta, "phi" = phi, "family"=family, "min.intknots" = min.intknots,
               "max.intknots" = max.intknots, "q" = q, "extr" = extr)
  previous <- matrix(nrow=max.intknots+1, ncol=max.intknots+4)
  nz <- if(!is.null(Z)) NCOL(Z) else 0
  oldcoef <- matrix(nrow=max.intknots+1, ncol=max.intknots+2+nz)
  oldguess <- matrix(nrow=max.intknots+1, ncol=max.intknots+2)
  nodi <- NULL
  Indicator <- table(X)
  Xw <- unique(X)
  iter <- NULL
  flag <- FALSE
  phis <- NULL
  phis_star <- NULL
  firstder <- NULL
  devi <- NULL
  for(j in 1:min(max.intknots+1,length(Y)-2)){
    if(flag) j <- j-2
    if(j >1)  {
      nodi<-sort(nodi)
      oldguess[j,1:(j+1)] <- guess
    }
    if(j == 1) guess <- NULL
    if(j>1) guess <- c(guess,guess_z)
    first.deg<-SplineReg_GLM(X=X,Y=Y,Z=Z,offset=offset,
                             weights=weights, extr=extr, InterKnots=nodi, n=2,
                             family = family, inits = guess) #first regression
    if(anyNA(first.deg$Theta)){
      rango <- Matrix::rankMatrix(first.deg$Basis)
      cols <- NCOL(first.deg$Basis)
      if(rango < cols){
        if(flag) {
          warning("Matrix singular for the second time. Breaking the loop.")
          break
        }
        check <- which(is.na(first.deg$Theta))
        nodi <- nodi[-(check+2)]
        guess <- guess[1:length(first.deg$Theta)][-check]
        toprint <- paste0("Basis Matrix singular, deleting one knot")
        print(toprint)
        flag  <- T
        if(cols == length(Xw)) {
          warning("Number of knots equal to number of unique Xs. Breaking the loop.")
          break
        } else {
          next}
      } else {
        stop("NA(s) in the coefficients")
      }
    } else {
      guess <- first.deg$Theta[1:(j+1)]
    }
    devi <- c(devi,first.deg$deviance)
    iter <- c(iter,length(devi))
    previous[j,1:(j+3)] <- sort(c(nodi,rep(extr,2)))
    guess_z <- if(nz>0) first.deg$Theta[(j+2):(j+1+nz)] else NULL
    oldcoef[j,1:(j+1+nz)] <- first.deg$Theta
    res.tmp <- first.deg$Residuals
    RSS.tmp <- first.deg$temporary$lastdev
    wi <- first.deg$temporary$weights  #rep(1,length(Y))#
    RSSnew <- c(RSSnew,RSS.tmp)
    prnt <- ""
    if(stoptype=="SR"){
      if(j>q  && length(phis)>=3){
        phis <- c(phis,RSSnew[j]/RSSnew[j-q])
        if(j-q > min.intknots){
          phismod <- log(1-phis)
          ccc<-.lm.fit(cbind(1,(q+1):j),phismod)$coef
          phis_star <- c(phis_star,1-exp(ccc[1])*exp(ccc[2]*j))
          oldintc <- c(oldintc,ccc[1])
          oldslp <- c(oldslp,ccc[2])
          prnt <- paste0(", phi_hat = ",round(1-exp(ccc[1])*exp(ccc[2]*j),3))

          if(1-exp(ccc[1])*exp(ccc[2]*j) >= phi) {
            break
          }
        }

      }
    }
    if(stoptype=="RD" | (stoptype=="SR" & length(phis)<3)){
      phis <- c(phis,RSSnew[j]/RSSnew[j-q])
      if(j>q && (j-q > min.intknots) ){
        prnt <- paste0(", phi = ",round(RSSnew[j]/RSSnew[j-q],3))
        if(RSSnew[j]/RSSnew[j-q] >= phi) {      # stop rule
          break
        }
      }
    }
    if(stoptype=="LR"){
      if(j>q && (j-q > min.intknots) ){
        phis <- c(phis,RSSnew[j-q]-RSSnew[j])
        prnt <- paste0(", p = ",round(pchisq(-(RSSnew[j]-RSSnew[j-q]),df=q),3))
        if(-(RSSnew[j]-RSSnew[j-q]) < qchisq(phi,df=q)) {      # stop rule
          break
        }
      }
    }
    d <- numeric()
    if(any(Indicator>1)){
      res.wkg <- makeNewRes2(resold = res.tmp,weights=weights*wi, recurr = as.numeric(Indicator)) #already fast - useless to do in c++
    } else {
      res.wkg <- res.tmp*weights*wi
    }
    segni <- sign(res.wkg)
    for(i in 1:length(Xw)){
      if(all(segni==segni[1])){
        d[i]<-length(segni)
        break
      } else {
        d[i]<-min(which(segni!=segni[1])-1)
        segni <- segni[-(1:d[i])]
      }
    }
    u<-length(d)
    dcum <- cumsum(d)
    wc.range<-numeric(u)
    means <- numeric(u)
    means[1] <- abs(mean(res.wkg[1:dcum[1]]))
    wc.range[1] <- Xw[dcum[1]]-X[1]
    for(i in 2:u){ # embed in C++ useless
      means[i] <- abs(mean((res.wkg[(dcum[i-1]+1):dcum[i]])))
      wc.range[i] <- Xw[dcum[i]]-Xw[dcum[i-1]+1]
    }
    means<-means/max(means)
    wc.range<-wc.range/max(wc.range)
    w <- beta*means + (1-beta)*wc.range+1e-8
    ris <- Knotnew(wht=w, restmp=res.wkg, x=Xw, dcm=dcum, oldknots=c(rep(extr,2+1),nodi), tol = tol)
    newknot <- ris[1]
    indice <- ris[2]
    quale <- sum(nodi < as.numeric(newknot))
    nk.design <- splineDesign(knots=sort(c(nodi,rep(extr,2))),derivs=0,x=newknot,ord=2,outer.ok = T)
    pr.value <- sum(nk.design*guess)
    newguess <- pr.value
    if(quale == 0){guess <- c(guess[1],newguess,guess[-1]) } else {
      if(quale == length(nodi)) {guess <- c(guess[1:(quale+1)],newguess,guess[quale+2]) } else {
        guess <- c(guess[1:(quale+1)],newguess,guess[-(1:(quale+1))])
      }
    }
    nodi<-c(nodi,newknot)

    if(show.iters) {
      if(j>q) {
        toprint<-paste0(indent, "Iteration ",j,": New Knot = ", round(newknot,3)  ,
                        ", DEV = " , round(RSSnew[j],3), prnt,"\n")
      } else {
        toprint<-paste0(indent,"Iteration ",j,": New Knot = ", round(newknot,3)  ,", DEV = " , round(RSSnew[j],3), "\n")}
      cat(toprint)}

  }
  if (j == max.intknots) warning("Maximum number of iterations exceeded")
  if (j<max.intknots) {
    previous <- previous[-((j+1):(max.intknots+1)),] #delete also last two$
    previous <- previous[,-((j+4):(max.intknots+4))]
    oldcoef <- oldcoef[-((j+1):(max.intknots+1)),]
    oldcoef <- oldcoef[,1:(j+1+nz)]
  }
  if(j-q<2){
    warning("Too few internal knots found: Linear spline will not be computed. Try to set a different value for 'q' or a different treshold")
    ll <- NULL
    lin <- NULL} else {
      ik <- previous[j-q,-c(1,2,(j+2-q):(j+3))]
      ll <- makenewknots(ik,2)#lin #l'approssimazione
      lin <- SplineReg_GLM(X=X,Y=Y,Z=Z,offset=offset,
                           weights=weights,extr=extr,InterKnots=ll,n=2,family=family,inits=c(oldguess[j-q,1:(j-q+1)],guess_z))
    }
  if(j-q<3){
    warning("Too few internal knots found: Quadratic spline will not be computed. Try to set a different value for 'q' or a different treshold")
    qq <- NULL
    squ <- NULL} else {
      qq <- makenewknots(ik,3)#quad
      squ <- SplineReg_GLM(X=X,Y=Y,Z=Z,offset=offset,
                           weights=weights,extr=extr,InterKnots=qq,n=3,family=family,mustart=family$linkinv(lin$Predicted))
    }
  if(j-q < 4){
    warning("Too few internal knots found: Cubic spline will not be computed. Try to set a different value for 'q' or a different treshold")
    cc <- NULL
    cub <- NULL} else {
      cc <- makenewknots(ik,4)#cub
      cub <- SplineReg_GLM(X=X,Y=Y,Z=Z,offset=offset,
                           weights=weights,extr=extr,
                           InterKnots=cc,n=4,family=family,mustart=family$linkinv(squ$Predicted))
    }
  out <- list("Type" = "GLM - Univ","Linear.Knots"=ll,"Quadratic.Knots"=qq,"Cubic.Knots"=cc,"Dev.Linear" = lin$RSS,
              "Dev.Quadratic" = squ$RSS,"Dev.Cubic" = cub$RSS,"Knots" = nodi,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = previous,
              "Args"= args,"Call"= save, "Nintknots"= j-q-1,"iters" = j,"Guesses" = oldguess,
              "Coefficients" = oldcoef,
              "deviance" = devi, "iterIrls"=iter, stopinfo = list("phis"=phis,"phis_star"=phis_star,
                                                                  "oldintc"=oldintc,"oldslp"=oldslp))
  class(out) <- "GeDS"
  return(out)
}







