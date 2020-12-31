
#include <Rcpp.h>
using namespace Rcpp;

//################################################################################
//' Mean Squared Prediction Errors, for a single \eqn{h}
//'
//' This function computes the estimated mean squared prediction errors from a
//' given time series and prediction coefficients
//'
//' The array of prediction coefficients \code{coef} is expected to be of
//' dimension \code{P x P x H x length(N) x length(t)} and in the format as
//' it is returned by the function \code{\link{predCoef}}. More precisely, for
//' \eqn{p=1,\ldots,P} and the \code{j.N}th element of \code{N} element of
//' \code{N} the coefficient of the
//' \code{h}-step ahead predictor for \eqn{X_{i+h}} which is computed from
//' the observations \eqn{X_i, \ldots, X_{i-p+1}} has to be available via
//' \code{coef[p, 1:p, h, j.N, t==i]}.
//'
//' Note that \code{t} have to be the indices corresponding to the coefficients.
//'
//' The resulting mean squared prediction error
//' \deqn{\frac{1}{|\mathcal{T}|} \sum_{t \in \mathcal{T}} (X_{t+h} - (X_t, \ldots, X_{t-p+1}) \hat v_{N[j.N],T}^{(p,h)}(t))^2}
//' is then stored in the resulting matrix at position \code{(p, j.N)}. 
//'
//' @name computeMSPEcpp
//'
//' @param X the data
//' @param coef the array of coefficients.
//' @param h which lead time to compute the MSPE for
//' @param t a vector of times from which backward the forecasts are computed
//' @param type indicating what type of measure of accuracy is to be computed;
//'            1: mspe, 2: msae
//' @param trimLo percentage of lower observations to be trimmed away
//' @param trimUp percentage of upper observations to be trimmed away
//'
//' @return Returns a \code{P x length(N)} matrix with the results.
//################################################################################

// [[Rcpp::export]]
NumericVector computeMSPEcpp(NumericVector X, NumericVector coef, int h, IntegerVector t, int type, double trimLo, double trimUp) {
  
   // get the dimensions of the coefficient array 
   IntegerVector d = coef.attr("dim");

   // allocate the array we will return
   NumericMatrix mspe( d[0], d[4]); // p x N  
   std::fill( mspe.begin(), mspe.end(), 0 ) ;

   for (int p = 0; p < d[0]; p++) {
      for (int N = 0; N < d[4]; N++) {
                 
         IntegerVector::iterator it_t;
         
         NumericVector mspeAddends( t.size() );
         std::fill( mspeAddends.begin(), mspeAddends.end(), 0 );
         
         int i_t = 0;
         
         for(it_t = t.begin(); it_t != t.end(); ++it_t) {
            
            // compute one forecast
            double Xhat = 0; 
            for (int i = 0; i <= p; i++) {
               int ii = p+i*d[0]+(h-1)*d[0]*d[1] + i_t*d[0]*d[1]*d[2] + N*d[0]*d[1]*d[2]*d[3];
               Xhat += X[*it_t - i - 1] * coef(ii);
            }
            if( type == 1 ) {
               mspeAddends[i_t] = (X[*it_t + h - 1] - Xhat) * (X[*it_t + h - 1] - Xhat);
            } else if( type == 2 ) {
               mspeAddends[i_t] = std::abs(X[*it_t + h - 1] - Xhat);
            }
            i_t++;
         }
         
         
         std::sort(mspeAddends.begin(), mspeAddends.end());
         int lo = floor(mspeAddends.size() * trimLo) + 1;
         int hi = mspeAddends.size() - floor(mspeAddends.size() * trimUp);

         mspe(p, N) += mean(mspeAddends[seq(lo-1, hi-1)]);
         
         // mspe(p, N) = mspe(p, N) / t.size();
      }
   }
   
   return mspe; // List::create(m, gamma, sigma);
}
