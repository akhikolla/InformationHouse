#ifndef GUARD_detectvectors_cpp
#define GUARD_detectvectors_cpp

#include "detectVectors.h"
//#include<Rcpp.h>


//---------------------------------------------------------------------//
//FFFcd
//---------------------------------------------------------------------//

//' Search for multiple changepoints in the mean using FFF
//'
//' Given a vector \code{x}, a threshold \code{alpha}, a value \code{lambda},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param lambda Value for the fixed forgetting factor in \eqn{(0,1)}.
//'
//' @param alpha Value for the significance threshold in \eqn{(0,1)}.
//' 
//' @param BL Value for the burn-in length.
//' 
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectFFFMeanMultiple(Rcpp::NumericVector x, double lambda, double alpha, int BL){
    FFFChangeDetector fffcd(lambda, alpha, BL);
    return ( fffcd.detectMultiple(x) );
}



//' Find the first changepoint in the mean using FFF
//'
//' Given a vector \code{x}, a threshold \code{alpha}, a value \code{lambda},
//' and a burn-in length \code{BL}, returns a list containing the single 
//' changepoint. Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param lambda Value for the fixed forgetting factor in \eqn{(0,1)}.
//'
//' @param alpha Value for the significance threshold in \eqn{(0,1)}.
//' 
//' @param BL Value for the burn-in length.
//' 
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{The first changepoint found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectFFFMeanSingle(Rcpp::NumericVector x, double lambda, double alpha, int BL){
    FFFChangeDetector fffcd(lambda, alpha, BL);
    return ( fffcd.detectSingle(x) );
}



//' Find the first changepoint in the mean using FFF, assuming prechange known
//'
//' Given a vector \code{x}, a value \code{lambda}, a threshold \code{alpha},
//' and values for known prechange mean and variance, returns a list containing 
//' the single changepoint. Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param lambda Value for the fixed forgetting factor in \eqn{(0,1)}.
//'
//' @param alpha Value for the significance threshold in \eqn{(0,1)}.
//'
//' @param prechangeMean Value of known prechange mean.
//' 
//' @param prechangeSigma Value of known prechange standard deviation.
//' 
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{The index of the first changepoint found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectFFFMeanSinglePrechange(Rcpp::NumericVector x, double lambda, double alpha, double prechangeMean, double prechangeSigma){
    //no need to initialize BL, will be set to zero later
    FFFChangeDetector fffcd(lambda, alpha);
    return ( fffcd.detectSinglePrechangeKnown(x, prechangeMean, prechangeSigma) );
}









//---------------------------------------------------------------------//
//AFFcd
//---------------------------------------------------------------------//





//' Search for multiple changepoints in the mean using AFF
//'
//' Given a vector \code{x}, a threshold \code{alpha}, a step size \code{eta}, 
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param alpha Value for the significance threshold in \eqn{(0,1)}.
//' 
//' @param eta Value for the step size in \eqn{(0,1)}.
//'
//' @param BL Value for the burn-in length.
//' 
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectAFFMeanMultiple(Rcpp::NumericVector x, double alpha, double eta, int BL){
    AFFChangeDetector affcd(alpha, eta, BL);
    return ( affcd.detectMultiple(x) );
}



//' Find the first changepoint in the mean using AFF
//'
//' Given a vector \code{x}, a threshold \code{alpha}, a step size \code{eta}, 
//' and a burn-in length \code{BL}, returns a list containing the single 
//' changepoint. Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//' 
//' @param alpha Value for the significance threshold in \eqn{(0,1)}.
//' 
//' @param eta Value for the step size in \eqn{(0,1)}.
//'
//' @param BL Value for the burn-in length.
//' 
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{The index of the first changepoint found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectAFFMeanSingle(Rcpp::NumericVector x, double alpha, double eta, int BL){
    AFFChangeDetector fffcd(alpha, eta, BL);
    return ( fffcd.detectSingle(x) );
}



//' Find the first changepoint in the mean using AFF, assuming prechange known
//'
//' Given a vector \code{x}, a threshold \code{alpha}, a step size \code{eta}, 
//' and a burn-in length BL, returns a list containing the single changepoint.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param alpha Value for the significance threshold in \eqn{(0,1)}.
//' 
//' @param eta Value for the step size in \eqn{(0,1)}.
//' 
//' @param prechangeMean Value of known prechange mean.
//' 
//' @param prechangeSigma Value of known prechange standard deviation.
//' 
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{The index of the first changepoint found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectAFFMeanSinglePrechange(Rcpp::NumericVector x, double alpha, double eta, double prechangeMean, double prechangeSigma){
    //no need to initialize BL, will be set to zero later
    AFFChangeDetector affcd(alpha, eta, 0);
    return ( affcd.detectSinglePrechangeKnown(x, prechangeMean, prechangeSigma) );
}



//---------------------------------------------------------------------//
//CUSUMcd
//---------------------------------------------------------------------//


//' Search for multiple changepoints in the mean using CUSUM
//'
//' Given a vector \code{x}, control parameters \code{k} and \code{h},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param k control parameter for CUSUM
//'
//' @param h control parameter for CUSUM
//'
//' @param BL Value for the burn-in length.
//' 
//' @details CUSUM updates via: 
//'          \deqn{S_{j} = \max{0, S_{j-1} + (x_{j} - \mu)/ \sigma - k}}
//'          where \eqn{\mu} and \eqn{\sigma} are, respectively, the mean 
//'          and variance of the in-control stream, 
//'          \eqn{x_j} is the observation at time \eqn{j}
//'          and \eqn{k} 
//'          is a control parameter for CUSUM. Then, a change is signalled
//'          if \eqn{S_j > h}, where \eqn{h} is the other control parameter.
//'          This is the formulation for using CUSUM to detect an increase
//'          in the mean; there is a similar formulation for detecting a 
//'          decrease, and usually CUSUM is two-sided (monitors for an 
//'          increase and a decrease in the mean).
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectCUSUMMeanMultiple(Rcpp::NumericVector x, double k, double h, int BL){
    CusumChangeDetector cusumcd(k, h, BL);
    return ( cusumcd.detectMultiple(x) );
}






//' Find the first changepoint in the mean using CUSUM
//'
//' Given a vector \code{x}, control parameters \code{k} and \code{h},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param k control parameter for CUSUM
//'
//' @param h control parameter for CUSUM
//'
//' @param BL Value for the burn-in length.
//' 
//' @details CUSUM updates via: 
//'          \deqn{S_{j} = \max{0, S_{j-1} + (x_{j} - \mu)/ \sigma - k}}
//'          where \eqn{\mu} and \eqn{\sigma} are, respectively, the mean 
//'          and variance of the in-control stream, 
//'          \eqn{x_j} is the observation at time \eqn{j}
//'          and \eqn{k} 
//'          is a control parameter for CUSUM. Then, a change is signalled
//'          if \eqn{S_j > h}, where \eqn{h} is the other control parameter.
//'          This is the formulation for using CUSUM to detect an increase
//'          in the mean; there is a similar formulation for detecting a 
//'          decrease, and usually CUSUM is two-sided (monitors for an 
//'          increase and a decrease in the mean).
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectCUSUMMeanSingle(Rcpp::NumericVector x, double k, double h, int BL){
    CusumChangeDetector cusumcd(k, h, BL);
    return ( cusumcd.detectSingle(x) );
}



//' Find the first changepoint in the mean using CUSUM, assuming prechange 
//' known
//'
//' Given a vector \code{x}, control parameters \code{k} and \code{h},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param k control parameter for CUSUM
//'
//' @param h control parameter for CUSUM
//'
//' @param BL Value for the burn-in length.
//' 
//' @details CUSUM updates via: 
//'          \deqn{S_{j} = \max{0, S_{j-1} + (x_{j} - \mu)/ \sigma - k}}
//'          where \eqn{\mu} and \eqn{\sigma} are, respectively, the mean 
//'          and variance of the in-control stream, 
//'          \eqn{x_j} is the observation at time \eqn{j}
//'          and \eqn{k} 
//'          is a control parameter for CUSUM. Then, a change is signalled
//'          if \eqn{S_j > h}, where \eqn{h} is the other control parameter.
//'          This is the formulation for using CUSUM to detect an increase
//'          in the mean; there is a similar formulation for detecting a 
//'          decrease, and usually CUSUM is two-sided (monitors for an 
//'          increase and a decrease in the mean).
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectCUSUMMeanSinglePrechange(Rcpp::NumericVector x, double k, double h, double prechangeMean, double prechangeSigma){
    CusumChangeDetector cusumcd(k, h, 0);
    return ( cusumcd.detectSinglePrechangeKnown(x, prechangeMean, prechangeSigma) );
}





//---------------------------------------------------------------------//
//EWMAcd
//---------------------------------------------------------------------//


//' Search for multiple changepoints in the mean using EWMA
//'
//' Given a vector \code{x}, control parameters \code{r} and \code{L},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param r control parameter for EWMA. Must be in range \eqn{[0,1]}.
//'
//' @param L control parameter for EWMA.
//'
//' @param BL Value for the burn-in length.
//' 
//' @details EWMA updates via: 
//'          \deqn{Z_{j} = (1-r) Z_{j-1} + r x_{j}}
//'          where \eqn{\mu} is the mean of the in-control stream, 
//'          \eqn{x_j} is the observation at time \eqn{j} and \eqn{r} 
//'          is a control parameter for EWMA. Then, a change is signalled
//'          if \deqn{|Z_j - \mu|  > L \sigma_{Z_j}}, 
//'          where \eqn{L} is the other control parameter, and 
//'          \eqn{\sigma_{Z_j}} is a scaled version of the in-control
//'          variance \eqn{\sigma}.
//'          This is the formulation for using EWMA to detect an increase or
//'          decrease in the mean.
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectEWMAMeanMultiple(Rcpp::NumericVector x, double r, double L, int BL){
    EwmaChangeDetector ewmacd(r, L, BL);
    return ( ewmacd.detectMultiple(x) );
}




//' Find the first changepoint in the mean using EWMA
//'
//' Given a vector \code{x}, control parameters \code{r} and \code{L},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param r control parameter for EWMA. Must be in range \eqn{[0,1]}.
//'
//' @param L control parameter for EWMA.
//'
//' @param BL Value for the burn-in length.
//' 
//' @details EWMA updates via: 
//'          \deqn{Z_{j} = (1-r) Z_{j-1} + r x_{j}}
//'          where \eqn{\mu} is the mean of the in-control stream, 
//'          \eqn{x_j} is the observation at time \eqn{j} and \eqn{r} 
//'          is a control parameter for EWMA. Then, a change is signalled
//'          if \deqn{|Z_j - \mu|  > L \sigma_{Z_j}}, 
//'          where \eqn{L} is the other control parameter, and 
//'          \eqn{\sigma_{Z_j}} is a scaled version of the in-control
//'          variance \eqn{\sigma}.
//'          This is the formulation for using EWMA to detect an increase or
//'          decrease in the mean.
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectEWMAMeanSingle(Rcpp::NumericVector x, double r, double L, int BL){
    EwmaChangeDetector ewmacd(r, L, BL);
    return ( ewmacd.detectSingle(x) );
}



//' Find the first changepoint in the mean using EWMA, assuming prechange 
//' known
//'
//' Given a vector \code{x}, control parameters \code{r} and \code{L},
//' and a burn-in length \code{BL}, returns a list containing the changepoints.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param r control parameter for EWMA. Must be in range \eqn{[0,1]}.
//'
//' @param L control parameter for EWMA.
//'
//' @param BL Value for the burn-in length.
//' 
//' @details EWMA updates via: 
//'          \deqn{Z_{j} = (1-r) Z_{j-1} + r x_{j}}
//'          where \eqn{\mu} is the mean of the in-control stream, 
//'          \eqn{x_j} is the observation at time \eqn{j} and \eqn{r} 
//'          is a control parameter for EWMA. Then, a change is signalled
//'          if \deqn{|Z_j - \mu|  > L \sigma_{Z_j}}, 
//'          where \eqn{L} is the other control parameter, and 
//'          \eqn{\sigma_{Z_j}} is a scaled version of the in-control
//'          variance \eqn{\sigma}.
//'          This is the formulation for using EWMA to detect an increase or
//'          decrease in the mean.
//' 
//' @return A list with 
//' \describe{
//'             \item{\code{tauhat}}{A vector of the changepoints found.}
//'          }
//' 
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_detectEWMAMeanSinglePrechange(Rcpp::NumericVector x, double r, double L, double prechangeMean, double prechangeSigma){
    EwmaChangeDetector ewmacd(r, L, 0);
    return ( ewmacd.detectSinglePrechangeKnown(x, prechangeMean, prechangeSigma) );
}



#endif
