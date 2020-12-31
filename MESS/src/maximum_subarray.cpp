// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast computation of maximum sum subarray
//'
//' @description Fast computation of the maximum subarray sum of a vector using Kadane's algorithm. The implementation handles purely negative numbers.
//' @param x A vector
//' @return A list with three elements: sum (the maximum subarray sum), start (the starting index of the subarray) and end (the ending index of the subarray)
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' maximum_subarray(1:4)
//' 
//' maximum_subarray(c(-2, 1, -3, 4, -1, 2, 1, -5, 4))
//'  
//' maximum_subarray(rnorm(100000))
//'
//' @export
// [[Rcpp::export]]
List maximum_subarray(const arma::vec & x) {

  //  arma::colvec X(x.begin(), x.size(), false);
  // arma::uword n = X.size();
  arma::uword n = x.size();

  double max_so_far = x[0];
  double max_ending_here=0;
  arma::uword start=0;
  arma::uword end=0;
  arma::uword s=0;

  for (arma::uword i=0; i<n; i++) {

    max_ending_here += x[i] ;

    if (max_so_far < max_ending_here) {
      max_so_far = max_ending_here;
      start = s;
      end = i;
    }			

    if (max_ending_here < 0) {
      max_ending_here = 0;
      s = i+1;
    }
    
  }



  Rcpp::List RVAL =  Rcpp::List::create(Rcpp::Named("sum") = max_so_far,
                                        Rcpp::Named("start") = start+1,
                                        Rcpp::Named("end") = end+1);


   return(RVAL);


}

/*


   int max_so_far = INT_MIN, max_ending_here = 0, 
       start =0, end = 0, s=0; 
  
    for (int i=0; i< size; i++ ) 
    { 
        max_ending_here += a[i]; 
  
        if (max_so_far < max_ending_here) 
        { 
            max_so_far = max_ending_here; 
            start = s; 
            end = i; 
        } 
  
        if (max_ending_here < 0) 
        { 
            max_ending_here = 0; 
            s = i + 1; 
        } 
    } 
    cout << "Maximum contiguous sum is "
        << max_so_far << endl; 
    cout << "Starting index "<< start 
        << endl << "Ending index "<< end << endl; 



*/
