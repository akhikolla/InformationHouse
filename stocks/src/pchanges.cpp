#include <Rcpp.h>
using namespace Rcpp;

//' Lagged Proportion Changes
//'
//' Calculates proportion changes between subsequent (or lagged) elements of a 
//' vector.
//'
//' @param x Numeric vector.
//' @param lag Numeric value (e.g. 2 for differences between 1st and 3rd
//' element, 2nd and 4th, ...).
//'
//'
//' @return Numeric vector.
//'
//'
//' @examples
//' # Generate 10 values from N(0, 1)
//' x <- rnorm(10)
//' 
//' # Calculate vector of proportion changes between subsequent values
//' (y <- pchanges(x))
//' 
//' # Equivalent base R computation
//' len <- length(x)
//' p1 <- x[2: len] 
//' p2 <- x[1: (len - 1)] 
//' y2 <- p1 / p2 - 1
//' all.equal(y, y2)
//' 
//'
//'@export
// [[Rcpp::export]]
NumericVector pchanges(NumericVector x, int lag = 1) {
  int n = x.size();
  NumericVector y(n - lag);
  if (lag == 1) {
    double current = 0;
    double previous = x(0);
    for (int a = 1; a < n; ++a) {
      current = x(a);
      y(a - 1) = current / previous - 1;
      previous = current;
    }
  }
  else {
    for (int a = lag; a < n; ++a) {
      y(a - lag) = x(a) / x(a - lag) - 1;
    }
  }
  return(y);
}
