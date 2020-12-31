#' Compute cubic B-spline basis functions for the given wavenumbers.
#' 
#' This function computes penalised cubic B-splines using the method proposed by
#' Eilers & Marx (1996). The spline coefficients can be computed efficiently
#' using sparse matrix algebra, as described in Sect. 2.3.3 of Green & Silverman
#' (1994) and Appendix B of Ruppert, Wand & Carroll (2003).
#' 
#' @param V a \code{vector} of wavenumbers, \eqn{\Delta \tilde{\nu}}.
#' @param n.b the number of basis functions to use.
#' @param pen the smoothing penalty hyperparameter.
#' @param prec a constant scale factor.
#'   
#' @return a \code{list} containing:
#' \describe{
#'  \item{\code{basis}}{A dense nwl by n.b matrix containing the values of the basis functions.} 
#'  \item{\code{precision}}{A sparse n.b by n.b \code{dsCMatrix}, the inverse of the prior covariance.}
#'  \item{\code{distance}}{The distance between each knot \eqn{(cm^{-1})}.}
#'  \item{\code{knots}}{The knot locations.}
#' }
#' @seealso \code{\link[Matrix:sparseMatrix-class]{sparseMatrix}}
#' @importFrom Matrix Matrix
#' @references
#' Eilers, PHC & Marx, BD (1996) "Flexible smoothing with B-splines and
#' penalties," Statist. Sci.  11(2): 89--121,
#' DOI: \href{http://dx.doi.org/10.1214/ss/1038425655}{10.1214/ss/1038425655}
#' 
#' Green, PJ & Silverman, BW (1994) "Nonparametric Regression and
#' Generalized Linear Models: a roughness penalty approach" Chapman & Hall, Boca
#' Raton, FL, pp. 11--21.
#' 
#' Ruppert, D; Wand, MP & Carroll, RJ (2003) "Semiparametric Regression" CUP,
#' Cambridge, UK, pp. 336--340.
getBsplineBasis <- function(V, n.b, pen, prec=1e-8) {
  NK <- n.b - 3 + 7
  h<-(max(V) - min(V))/(NK-7) # distance between equally-spaced knots (cm^-1)
  V_Ex<-c(min(V)-(3*h):1,V,max(V)+1:(3*h))
  NS_Ex<-length(V_Ex)
  Knots<-seq(min(V_Ex),max(V_Ex),length=NK)
  
  B1<-matrix(rep(0,NK*NS_Ex),ncol=NS_Ex)
  for(i in 1:(NK-1)){
    for(j in 1:NS_Ex){
      if(V_Ex[j]>=Knots[i] && V_Ex[j]<Knots[i+1]){
        B1[i,j]<-1
      }
    }
  }
  B1[NK,NS_Ex]<-1
  
  B4<-matrix(rep(0,(NK-4)*NS_Ex),ncol=NS_Ex)
  for(i in 1:(NK-4)){
    for(j in 1:NS_Ex){
      B4[i,j]<-(1/(6*h^3))*((V_Ex[j]-Knots[i])^3*B1[i,j]+(V_Ex[j]-Knots[i])^2*(Knots[i+2]-V_Ex[j])*B1[i+1,j]+(V_Ex[j]-Knots[i])*(V_Ex[j]-Knots[i+1])*(Knots[i+3]-V_Ex[j])*B1[i+1,j]+(V_Ex[j]-Knots[i])*(Knots[i+3]-V_Ex[j])^2*B1[i+2,j]+(Knots[i+4]-V_Ex[j])*(V_Ex[j]-Knots[i+1])^2*B1[i+1,j]+(Knots[i+4]-V_Ex[j])*(V_Ex[j]-Knots[i+1])*(Knots[i+3]-V_Ex[j])*B1[i+2,j]+(Knots[i+4]-V_Ex[j])^2*(V_Ex[j]-Knots[i+2])*B1[i+2,j]+(Knots[i+4]-V_Ex[j])^3*B1[i+3,j])
    }
  }
  
  X<-t(B4[,1:length(V)+3*h])
  NB<-NK-4
  
  pre<-matrix(rep(0,NB^2),ncol=NB)
  pre[1,1]<-pre[NB,NB]<-1
  pre[2,2]<-pre[NB-1,NB-1]<-5
  for(i in 3:(NB-2)){
    pre[i,i]<-6
  }
  pre[1,2]<-pre[2,1]<--2
  pre[NB-1,NB]<-pre[NB,NB-1]<--2
  for(i in 2:(NB-2)){
    pre[i,i+1]<-pre[i+1,i]<--4
  }
  for(i in 1:(NB-2)){
    pre[i,i+2]<-pre[i+2,i]<-1
  }
  
  g0<-prec*Matrix(diag(rep(1,NB))+pen*pre,sparse=TRUE)
  list(basis=X, precision=g0, distance=h, knots=Knots)
}
