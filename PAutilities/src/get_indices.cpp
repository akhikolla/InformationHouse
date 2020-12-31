#include <Rcpp.h>
using namespace Rcpp;

//' Retrieve indices for a rolling window analysis
//'
//' @param y_var NumericVector. Input on which to define the indices for each roll of
//'   the window
//' @param window_size int. The size of the window
//' @return a list in which each element contains \code{window_size} consecutive
//'   integers that indicate which elements of \code{y_var} would be extracted
//'   for that roll of the window
//' @export
//' @note For this function, the output elements contain positions (i.e., indices) from
//'   \code{y_var}, whereas for \code{\link{rolling_groups}} the output elements
//'   contain the raw values found at each index
//' @seealso \code{\link{rolling_groups}}
//' @examples
//' result <- get_indices(1:100, 10)
//' head(result)
//' tail(result)
// [[Rcpp::export]]
List get_indices(NumericVector y_var, int window_size = 15) {
  //--window_size;
  int n = y_var.size() - window_size + 1;
  List indices(n);
  for (int i = 0; i < n; ++i) {
    indices[i] = seq(i+1, i+window_size);
  }
  return indices;
}
