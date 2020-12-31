
#include <Rcpp.h>
using namespace Rcpp;

//################################################################################
//' \eqn{h}-step Prediction coefficients
//'
//' This function computes the localised and iterated Yule-Walker coefficients
//' for h-step ahead forecasting of \eqn{X_{t+h}} from \eqn{X_{t}, ..., X_{t-p+1}},
//' where \eqn{h = 1, \ldots,} \code{H} and \eqn{p = 1, \ldots,} \code{P}.
//' 
//' For every \eqn{t \in} \code{t} and every \eqn{N \in} \code{N} the (iterated) Yule-Walker
//' estimates \eqn{\hat v_{N,T}^{(p,h)}(t)} are computed. They are defined as
//' \deqn{\hat v_{N,T}^{(p,h)}(t) :=  e'_1 \big( e_1 \big( \hat a_{N,T}^{(p)}(t) \big)' + H \big)^h, \quad N \geq 1,}
//' and
//' \deqn{\hat v_{0,T}^{(p,h)}(t) :=  \hat v_{t,T}^{(p,h)}(t),}
//' with
//' \deqn{ e_1 := \left(\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array} \right), \quad H := \left( \begin{array}{ccccc} 0 & 0 & \cdots & 0 & 0 \\ 1 & 0 & \cdots & 0 & 0 \\ 0 & 1 & \cdots & 0 & 0 \\ \vdots & \ddots & \cdots & 0 & 0 \\ 0 & 0 & \cdots & 1 & 0 \end{array} \right)}
//' and
//' \deqn{ \hat a_{N,T}^{(p)}(t) := \big( \hat\Gamma_{N,T}^{(p)}(t) \big)^{-1} \hat\gamma_{N,T}^{(p)}(t),}
//' where
//' \deqn{\hat\Gamma_{N,T}^{(p)}(t) := \big[ \hat \gamma_{i-j;N,T}(t) \big]_{i,j = 1, \ldots, p}, \quad \hat \gamma_{N,T}^{(p)}(t) := \big( \hat \gamma_{1;N,T}(t), \ldots, \hat \gamma_{p;N,T}(t) \big)'}
//' and
//' \deqn{\hat \gamma_{k;N,T}(t) := \frac{1}{N} \sum_{\ell=t-N+|k|+1}^{t} X_{\ell-|k|,T} X_{\ell,T}}
//' is the usual lag-\eqn{k} autocovariance estimator (without mean adjustment),
//' computed from the observations \eqn{X_{t-N+1}, \ldots, X_{t}}.
//'
//' The Durbin-Levinson Algorithm is used to successively compute the solutions to the
//' Yule-Walker equations (cf. Brockwell/Davis (1991), Proposition 5.2.1).
//' To compute the \eqn{h}-step ahead coefficients we use the recursive relationship
//' \deqn{\hat v_{i,N,T}^{(p)}(t,h) = \hat a_{i,N,T}^{(p)}(t) \hat v_{1,N,T}^{(p,h-1)}(t) + \hat v_{i+1,N,T}^{(p,h-1)}(t) I\{i \leq p-1\},}
//' (cf. Section 3.2, Step 3, in Kley et al. (2019)).
//'
//' @name predCoef
//' @export
//'
//' @param X the data \eqn{X_1, \ldots, X_T}
//' @param P the maximum order of coefficients to be computed; has to be a positive integer
//' @param H the maximum lead time; has to be a positive integer
//' @param t a vector of values \eqn{t}; the elements have to satisfy
//'            \code{max(t) <= length(X)} and  \code{min(t) >= min(max(N[N != 0]),p)}.
//' @param N a vector of values \eqn{N}; the elements have to satisfy
//'            \code{max(N[N != 0]) <= min(t)} and \code{min(N[N != 0]) >= 1 + P}.
//'            \eqn{N = 0} corresponds to the case where all data is taken into account.
//'
//' @return Returns a named list with elements \code{coef}, \code{t}, and \code{N},
//'         where \code{coef} is an array of dimension
//'         \code{P} \eqn{\times} \code{P} \eqn{\times} \code{H} \eqn{\times}
//'         \code{length(t)} \eqn{\times} \code{length(N)}, and
//'         \code{t}, and \code{N} are the parameters provided on the call of the
//'         function. See the example on how to access the vector
//'         \eqn{\hat v_{N,T}^{(p,h)}(t)}.
//'
//' @references
//' Brockwell, P. J. & Davis, R. A. (1991).
//' Time Series: Theory and Methods. Springer, New York.
//'
//' @examples
//' T <- 100
//' X <- rnorm(T)
//'
//' P <- 5
//' H <- 1
//' m <- 20
//'
//' Nmin <- 25
//' pcoef <- predCoef(X, P, H, (T - m - H + 1):T, c(0, seq(Nmin, T - m - H, 1)))
//' 
//' ## Access the prediction vector for p = 2, h = 1, t = 95, N = 25
//' p <- 2
//' h <- 1
//' t <- 95
//' N <- 35
//' res <- pcoef$coef[p, 1:p, h, pcoef$t == t, pcoef$N == N]
//################################################################################

int vecmin(IntegerVector x) {
  // Rcpp supports STL-style iterators
  IntegerVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}

// [[Rcpp::export]]
List predCoef(NumericVector X, int P, int H, IntegerVector t, IntegerVector N) {

   // verify validity of parameters
   if (P < 1) {
      stop("P has to be positive");
   }
   if (H < 1) {
      stop("H has to be positive");
   }
   
   IntegerVector::iterator it_t;
   IntegerVector::iterator it_N;
   
   int T = X.size();
   
   for(it_t = t.begin(); it_t != t.end(); ++it_t) {
      for(it_N = N.begin(); it_N != N.end(); ++it_N) {
         if ( *it_t > T ) {
            stop("the elements of t have to satisfy: max(t) <= length(X)");
         }
         if ( *it_N != 0 && *it_N < 1+P ) {
            stop("the elements of N have to satisfy: min(N[N != 0]) >= 1+P");
         }
         if ( *it_N > *it_t ) {
            stop("the elements of N and t have to satisfy: max(N) <= min(t)");
         }
      }
   }
   
   
/*    P >= 1
      H >= 1
      min(t) >= min(max(N),p) 
      max(t) <= lenght(X)-H}
      min(t) >= min(max(N),p)
      min(N) >= 1+P */

   // allocate the array we will return
   NumericVector coef( t.size() * N.size() * H * P * P );
   int tsize = t.size();
   int Nsize = N.size();
   IntegerVector d = IntegerVector::create( P, P, H, tsize, Nsize);
   coef.attr("dim") = d;

   std::fill( coef.begin(), coef.end(), NumericVector::get_na() ) ;

   int i_t = 0;
   for(it_t = t.begin(); it_t != t.end(); ++it_t) {
      int i_N = 0;
      for(it_N = N.begin(); it_N != N.end(); ++it_N) {
         
         IntegerVector idx;
         if (*it_N > 0) {
            idx = *it_t - *it_N - 1 + seq_len(*it_N);        
         } else {
            idx = seq_len(*it_t) - 1;
         }
         NumericVector x = X[idx];
         // Rcout << "The value of x: " << x << std::endl;
         
         int n = x.size();
   
         // BEGIN Durbin-Levinson algorithm 
         std::vector< std::vector<double> > rmat(P, std::vector<double>(P,0));

         NumericVector gamma(P+1);
         NumericVector sigma(P);
         
         for (int k = 0; k < P+1; k++) {
            double res = 0;
            for (int j = 0; j < n-k; j++) {
               res += x[j]*x[j+k];
            }
            //gamma[k] = res / (n-k);
            gamma[k] = res / n;
         }
         
         // Iteration 0
         rmat[0][0] = gamma[1] / gamma[0];
         sigma[0] = gamma[0] * (1 - rmat[0][0]*rmat[0][0]);
        
         for (int j = 2; j <= P; j++) {
            
            rmat[j-1][j-1] = gamma[j];
            for (int i = 1; i < j; i++) { // i here corresponds to i+1 in the notes
               rmat[j-1][j-1] -= rmat[j-1-1][i-1] * gamma[j-i];
            }
            rmat[j-1][j-1] /= sigma[j-1-1];
                  
            for (int k = 1; k < j; k++) {
               rmat[j-1][k-1] = rmat[j-1-1][k-1] - rmat[j-1][j-1] * rmat[j-1-1][j-k-1];
               // Rcout << "The value is " << j << " " << k << std::endl;
               // Rcout << "The value is " << rmat[j][k] << std::endl;
               // rmat[j][k] = gamma[k];
            }
            
            sigma[j-1] = sigma[j-2] * (1 - rmat[j-1][j-1] * rmat[j-1][j-1]);
         }
         
         // END Durbin-Levinson algorithm
   
         std::vector< std::vector<double> > vmat(P, std::vector<double>(P,0));
         for (int i = 0; i < P; i++) {
            for (int j = 0; j <= i; j++) {
               vmat[i][j] = rmat[i][j];
               coef(i+j*P+0*P*P + i_t*P*P*H + i_N*P*P*H*t.size()) = rmat[i][j];
            }
         }
         for (int h = 1; h < H; h++) {
            for (int i = 0; i < P; i++) {
               double v1 = vmat[i][0];
               for (int j = 0; j <= i; j++) {  
                  vmat[i][j] = rmat[i][j] * v1;
                  if (j < i) { vmat[i][j] += vmat[i][j+1]; }
                  coef(i+j*P+h*P*P + i_t*P*P*H + i_N*P*P*H*t.size()) = vmat[i][j];               
               }
            }
         }

         i_N++;
      }
      i_t++;
   }
   
   // return coef; // List::create(coef, gamma, sigma);
   return List::create(Rcpp::Named("coef") = coef,
                       Rcpp::Named("t") = t,
                       Rcpp::Named("N") = N);
}
