##  Copyright (C) 2015-2016 Leo Belzile
##
##  This file is part of the "lcopula" package for R.  This program
##  is free software; you can redistribute it and/or modify it under the
##  terms of the GNU General Public License as published by the Free
##  Software Foundation; either version 2 of the License, or (at your
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
##  MA 02111-1307 USA or look up the web page
##  http://www.gnu.org/copyleft/gpl.html.


##########################################################################
##              Pickands functions for CDA of Liouville vectors         ##
##     Spectral density for both CDA of copula and survival copula      ##
##                        Version of August 16th, 2015                  ##
##                      Last modified - August 16th, 2015               ##
##########################################################################

#' Pickands dependence function for the copula domain of attraction of Liouville copulas
#'
#' Pickands dependence function as in \cite{Belzile (2014), Proposition 41}
#' Returns the Pickands dependence function of the copula domain of attraction
#' of the copula. This is only derived and implemented in the bivariate case.
#'
#' @param t pseudo-angle in (0,1)
#' @param rho index of regular variation parameter
#' @param alpha vector of Dirichlet allocations. Currently must be of length 2
#'
#' @return value of Pickands function
.pickands.fun.uni<-function(t, rho=0.5, alpha=c(1,1)) {
  if(t==0 && t==1){return(1)}
  a <- lgamma(alpha-rho)-lgamma(alpha)
  kappa <- exp((log(1-t)+a[2])/rho - log(exp((log(1-t)+a[2])/rho)+exp((log(t)+a[1])/rho)))
  exp(log(1-t)+ pbeta(kappa, alpha[1]-rho, alpha[2],log.p=TRUE))+
    exp(log(t)+ pbeta(1-kappa, alpha[2]-rho,alpha[1],log.p=TRUE))

}

.pickands.fun<-Vectorize(.pickands.fun.uni, vectorize.args=c("t"))

#' Pickands dependence function for the copula domain of attraction of Liouville survival copulas
#'
#' Pickands dependence function as in \cite{Belzile (2014), Proposition 40 and Example 4}
#' Returns the Pickands dependence function of the copula domain of attraction
#' of the survival copula, the scaled Dirichlet extreme value model. Currently only implemented in the bivariate case.
#' Setting \code{rho=1} yields the same output as the function in the \code{evd} package.
#'
#' @param t pseudo-angle in (0,1)
#' @param rho index of regular variation parameter
#' @param alpha vector of Dirichlet allocations. Currently must be of length 2
#'
#' @return value of Pickands function for the scaled extremal Dirichlet  model
.pickands.dir.uni<-function(t,alpha,rho){
  if(t==0 && t==1){return(1)}
  k1=lbeta(alpha[1]+rho,alpha[2])
  k2=lbeta(alpha[1],alpha[2]+rho)
  kappa<-exp((1/rho)*(log(1-t)+k1)-log(exp((1/rho)*(log(t)+k2))+exp((1/rho)*(k1+log(1-t)))))
  exp(log(1-t)+ pbeta(kappa, alpha[1], alpha[2]+rho,log.p=T))+
    exp(log(t)+ pbeta(1-kappa, alpha[2],alpha[1]+rho,log.p=T))
}
.pickands.dir<-Vectorize(.pickands.dir.uni, "t")

#' Pickands dependence function for the copula domain of attraction of Liouville survival copulas
#'
#' Pickands dependence function as in \cite{Belzile (2014), Proposition 40 and Example 4} and
#' \cite{Belzile (2014), Proposition 41}, assuming that the parameter alpha is integer-valued.
#' Returns the Pickands dependence function of the copula domain of attraction (CDA)
#' of the survival copula, the scaled Dirichlet extreme value model, or the CDA of the copula, the Liouville EV model.
#'
#' @param t pseudo-angle in (0,1)
#' @param rho index of regular variation parameter
#' @param alpha vector of Dirichlet allocations. Currently must be of length 2
#' @param CDA select the extremal attractor of the copula (\code{C}) or the survival copula (\code{S})
#'
#' @examples
#' pickands.liouv(seq(0,1,by=0.01),1,c(0.1,0.3),CDA="S")
#' pickands.liouv(t = seq(0,1,by=0.01), rho = 0.5, alpha = c(1,3), CDA="C")
#' @return value of Pickands function for the scaled Dirichlet EV model
#' @export
pickands.liouv<-function(t, rho=0.5, alpha=c(1,1),CDA=c("C","S")){
  #Check conditions
  if(missing(CDA)==TRUE){
    CDA="C"
    warning("Setting default to Liouville CDA. Use CDA=`S` for the Dirichlet model")
  }
  if(any(t<0) || any(t>1)){
    stop("t must be a vector whose elements are between 0 and 1")
  }
  if(any(rho<=0)){
    stop("rho must be positive")
  }
  if(length(rho)>1){
    stop("Argument `rho' must be of length 1")
  }
  if(CDA=="C" && any(c(rho >= min(alpha), rho <= 0))){
    stop("`rho must be between (0,min(alpha))")
  }
  if(any(alpha<=0)){
    stop("alpha must be positive")
  }
  if(length(alpha)!=2){
    stop("Not implemented except in the bivariate case")
  }
  if(CDA=="C"){
    return(.pickands.fun(t,rho=rho,alpha=alpha))
  } else if(CDA=="S"){
    return(.pickands.dir(t,rho=rho,alpha=alpha))
  }  else{
    return(NA)
  }
}
#' Plot Pickands dependence function for CDA of Liouville copulas
#'
#' The function will draw the Pickands dependence function for output in \code{tikz} if
#' the corresponding function is selected.
#'
#' @param rho index of regular variation parameter
#' @param alpha vector of Dirichlet allocations. Currently must be of length 2
#' @param plot.new boolean indicating whether a new plotting device should be called
#' @param CDA whether to plot Pickands function for the extremal model of the
#' copula (\code{C}) or the survival copula (\code{S}), which is the scaled Dirichlet
#' @param tikz boolean specifying whether to prepare plot for \code{tikz} output. Defaults to \code{F}
#' @param ... additional arguments passed to \code{lines}
#'
#' @return a plot of the Pickands dependence function
#' @export
#' @examples
#' pickands.plot(rho=0.9, alpha=c(1,1), col="slateblue1", CDA="C")
#' pickands.plot(rho=0.9, alpha=c(2,3), col="slateblue2", CDA="C", plot.new=FALSE)
#' pickands.plot(rho=0.5, alpha=c(2,3), col="slateblue3", CDA="C", plot.new=FALSE)
#' #Parameters for the Pickands function of the scaled Dirichlet need not be integer
#' pickands.plot(rho=0.9, alpha=c(1,1), CDA="S")
#' pickands.plot(rho=0.9, alpha=c(0.2,0.5), col="darkred", CDA="S", plot.new=FALSE)
#' pickands.plot(rho=0.8, alpha=c(1.2,0.1), col="red", CDA="S", plot.new=FALSE)
pickands.plot<-function(rho, alpha, plot.new=T,CDA=c("C","S"), tikz=F, ...){
  #Check conditions
  if(missing(CDA)==TRUE){
    CDA="C"
    warning("Setting default to Liouville CDA. Use CDA=`S` for the Dirichlet model")
  }
  if(length(rho)!=1){stop("rho must be 1-dimensional")}
  if(length(alpha)!=2){
    stop("Not implemented beyond bivariate case")
  }
  if(plot.new==TRUE){
    plot.new()
    plot.window(c(0,1), c(0.5,1))
    axis(side=2, at=seq(0.5,1,by=0.1), pos=0,las=2,tck=0.01)
    axis(side=1, at=seq(0,1,by=0.1), pos=0.5,las=0,tck=0.01)
    #title("Pickands dependence function")
    #Empirical bounds
    lines(c(0,0.5),c(1,0.5),lty=3,col="gray")
    lines(c(0.5,1),c(0.5,1),lty=3,col="gray")
    lines(c(0,1),c(1,1),lty=3,col="gray")
    if(tikz==TRUE){
      mtext("$t$", side=1, line=2)
      mtext("$\\mathrm{A}(t)$", side=2, line=2)
    } else{
      mtext(expression(t), side=1, line=2)
      mtext(expression(A(t)), side=2, line=2)
    }

  }
  x = seq(0,1,by=0.001)
  if(CDA=="C"){
    lines(x=x,y=c(1,.pickands.fun(x[-c(1,length(x))],alpha=alpha, rho=rho),1),type="l",...)
  }
  if(CDA=="S"){
    lines(x=x,y=.pickands.dir(x,alpha=alpha,rho=rho),type="l",...)
  }
}
#' Kendall plot
#'
#' This function plots the expectation of the order statistics under the null
#' hypothesis of independence against the ordered empirical copula values. The
#' data is transformed to ranks.
#'
#' The function uses \code{integrate} and may fail for large \code{d} or large \code{n}. If \eqn{n>200}, the fallback is to generate
#' a corresponding sample of uniform variates and to compare the empirical copula of the sample generated under the null hypothesis with the one
#' obtained from the sample.
#'
#' @author Pr. Christian Genest (the code was adapted for the multivariate case)
#'
#' @param data a \code{n} by \code{d} matrix of observations
#' @param add whether to superimpose lines to an existing graph. Default to \code{F}
#' @param ... additional arguments passed to \code{points}
#'
#' @return The Kendall plot corresponding to the data at hand
#' @export
#' @references Genest & Boies (2003). Detecting Dependence with Kendall Plots, \emph{The American Statistician}, \bold{57}(4), 275--284.
#' @examples
#' #Independence
#' K.plot(matrix(runif(2000),ncol=2))
#' #Negative dependence
#' K.plot(rCopula(n=1000,claytonCopula(param=-0.5,dim=2)),add=TRUE,col=2)
#' #Perfect negative dependence
#' K.plot(rCopula(n=1000,claytonCopula(param=-1,dim=2)),add=TRUE,col=6)
#' #Positive dependence
#' K.plot(rCopula(n=1000,claytonCopula(param=iTau(claytonCopula(0.3),0.5),dim=2)),add=TRUE,col=3)
#' #Perfect positive dependence
#' K.plot(rCopula(n=1000,claytonCopula(param=iTau(claytonCopula(0.3),1),dim=2)),add=TRUE,col=4)
K.plot <- function(data, add=F, ...){
  if(is.null(dim(data)) || ncol(data)<2){
    stop("Invalid input matrix")
  }
  #Empirical copula function
  H <- function(data){
    n <- dim(data)[1]
    d <- dim(data)[2]
    out <- rep(NA,n)
    tmp <- paste("(data[,",1:d,"] <= data[i,",1:d,"])")
    cmnd <- parse(text=paste(tmp,collapse=" & "))
    for (i in 1:n){
      subb <- sum(eval(cmnd))-1
      out[i] <- subb/(n-1)
    }
    out
  }
  # Expectation of order statistics for sample of dimension n x d
  W <- function(n, d){
    fun <- function(w,i,d){
      K0 <- w+w*rowSums(sapply(1:(d-1), function(k){
        exp(k*log(-log(w))-lgamma(k+1))
      }))
      dK0 <- exp((d-1)*log(-log(w))-lgamma(d))
      exp(log(w)+log(dK0)+(i-1)*log(K0)+(n-i)*log(1-K0))
    }
    sapply(1:n, function(i){
      n*choose(n-1,i-1)*integrate(fun,lower=0,upper=1,i=i,d=d,subdivisions=10000,
                                  rel.tol=.Machine$double.eps^0.001)$value
    })
  }
  # Plotting function
  n <- nrow(data)
  d <- ncol(data)
  W <- W(n, d)
  Hs <- sort(H(data))
  seq <- seq(from=0.01,to=1,by=0.01)
  # If n<200, plot the points at location of ranks in the sample
  if(n < 200) {
    if(add==F){
      plot(W,Hs,xlim=c(0,1),ylim=c(0,1),xlab="Independence",ylab="Data",pch=4,cex=0.8)
      lines(c(0,1),c(0,1),lty=2)
      K0 <- seq-seq*log(seq)
      lines(c(0,seq),c(0,K0),lty=2)
      lines(c(1,0),c(0,0),lty=2)
    } else{
      points(W,Hs,xlim=c(0,1),ylim=c(0,1),pch=4,cex=0.8, ...)
    }
  }  else {
    #Generate random uniforms and compute
    ind <- matrix(runif(n*d),nrow=n,ncol=d)
    Wind <- sort(H(ind))
    if(add==F){
      plot(Wind,Hs,xlim=c(0,1),ylim=c(0,1),xlab="Independence",ylab="Data",pch=4,cex=0.5)
      lines(c(0,1),c(0,1),lty=2)
      #This bound is valid for d=2, but the lines are always drawn
      K0 <- seq-seq*log(seq)
      lines(c(0,seq),c(0,K0),lty=2)
      lines(c(1,0),c(0,0),lty=2)
    } else{
      points(Wind, Hs,xlim=c(0,1),ylim=c(0,1),pch=4,cex=0.5, ...)
    }
  }
}

#' Spectral density of the CDA of survival copula and copula of Liouville vectors
#'
#' Computes the Liouville EV model or the scaled Dirichlet EV model spectral density
#'
#' @param w matrix of points in the unit simplex at which to evaluate the density
#' @param alpha vector of Dirichlet allocations (strictly positive).
#' @param rho parameter of limiting model corresponding to index of regular variation
#' @param CDA copula domain of attraction of either the Liouville copula, \code{C}, or its survival copula \code{S}
#' @param logdensity logical; whether to return the log density or not
#' @return a vector with the same number of rows as \code{w}.
#' @export
#' @examples
#' hmvevdliouv(seq(0.01,0.99,by=0.01), alpha=c(1,2), rho=0.2, CDA="C")
#' hmvevdliouv(seq(0.01,0.99,by=0.01), alpha=c(0.1,2), rho=0.2, CDA="S")
#' hmvevdliouv(seq(0.01,0.99,by=0.01), alpha=c(1,2), rho=0.2, CDA="S")
hmvevdliouv <- function(w, alpha, rho, CDA=c("C","S"), logdensity=FALSE){
  if(any(w<0) || any(w>1)){
    stop("w must contain elements between 0 and 1")
  }
  w <- as.matrix(w)
  if(ncol(w)!=length(alpha)){
    w <- cbind(w, 1-rowSums(w))
  } else{
    if(rowSums(w)!=rep(1, nrow(w))){
    stop("Components must sum to 1")
    }
  }
  if(length(rho)>1){
    rho <- rho[1]
    warning("Using only first coordinate of rho")
  }
  if(rho<=0){
    CDA <- "C"
    rho <- -rho
  }
  if(CDA=="C" && rho> min(alpha)){
    stop("`abs(rho)` must be smaller than `min(alpha)`")
  }
  if(any(alpha<=0)){
    stop("alpha must be positive")
  }
  if(missing(CDA)){
    warning("Invalid input for argument `CDA`. Defaulting to `C`, the `negdir` family")
    CDA = "C"
  }
  if(!(CDA %in% c("C","S"))){
    stop("Invalid input for argument `CDA`")
  }
  
  h <- apply(w, 1, function(wrow){ switch(CDA, 
    C=negdirspecdens(dat=t(as.matrix(wrow)), param=c(alpha, rho), d=length(alpha), transform=FALSE),
    S=dirspecdens(dat=t(as.matrix(wrow)),param=c(alpha, rho), d=length(alpha), transform=FALSE)
   )})
  if(logdensity){
    return(h)
  } else{
    return(exp(h))
  }
}
