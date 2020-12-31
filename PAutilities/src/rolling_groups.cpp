#include <Rcpp.h>
using namespace Rcpp;

//' Loop along a vector, returning n elements at a time in a list
//'
//' @param values IntegerVector. The vector to loop along
//' @param n int. The number of elements to return in each element of the
//'   resulting list
//' @return a list in which each element contains \code{n} elements from
//'   \code{values}
//' @export
//' @note For this function, the output elements contain raw values from
//'   \code{values}, whereas for \code{\link{get_indices}} the output elements
//'   contain the positions (i.e., indices) rather than the raw values
//' @seealso \code{\link{get_indices}}
//' @examples
//' groups <- rolling_groups(0:50, 3)
//' head(groups)
//' tail(groups)
// [[Rcpp::export]]
List rolling_groups(IntegerVector values, int n = 2) {

  List result(values.size() - (n - 1));

  for (int i = 0; i < result.size(); ++i) {

    IntegerVector group(n);
    for (int j = 0; j < group.size(); ++j) {
      group[j] = values[i + j];
    }

    result[i] = group;

  }

  return result;

}
