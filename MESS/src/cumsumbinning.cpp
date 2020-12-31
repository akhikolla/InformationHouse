#include <Rcpp.h>
using namespace Rcpp ;

//' Binning based on cumulative sum with reset above threshold
//' 
//' Fast binning of cumulative vector sum with new groups when the sum passes a threshold or the group size becomes too large
//'
//' Missing values (NA, Inf, NaN) are completely disregarded and pairwise complete cases are used f
//' 
//' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param threshold The value of the threshold that the cumulative group sum must not cross OR the threshold that each group sum must pass (when the argument cuwhatpassed is set to TRUE). 
//' @param cutwhenpassed A boolean. Should the threshold be the upper limit of the group sum (the default) or the value that each group sum needs to pass (when set to TRUE).
//' @param maxgroupsize An integer that defines the maximum number of elements in each group. NAs count as part of each group but do not add to the group sum. NULL (the default) corresponds to no group size limits.
//' @return An integer vector giving the group indices
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' set.seed(1)
//' x <- sample(10, 20, replace = TRUE)
//' cumsumbinning(x, 15)
//' cumsumbinning(x, 15, 3)
//' 
//' x <- c(3, 4, 5, 12, 1, 5, 3)
//' cumsumbinning(x, 10)
//' cumsumbinning(x, 10, cutwhenpassed=TRUE)
//'
//' @export
// [[Rcpp::export]]
IntegerVector cumsumbinning(NumericVector x, double threshold, bool cutwhenpassed=false, Rcpp::Nullable<int> maxgroupsize=R_NilValue) {
  IntegerVector groupVec (x.size());
  int group = 1;
  int groupsize=0;
  double runSum = 0;

  int maxgroup;

  // Sanity checks
  if ((maxgroupsize.isNull())) {
    maxgroup=0;
  }
  else {     
    NumericVector xxx(maxgroupsize.get());
    maxgroup = xxx[0];
    if (maxgroup<1) {
      stop("maxgroupsize should be larger than 0");
    }
  }
  
  
  for (int i = 0; i < x.size(); i++) {
    if (!NumericVector::is_na(x[i]))
      runSum += x[i];

    groupsize++;

    // Set the group for the current 
    // observation
    groupVec[i] = group;

    if (runSum > threshold) {

      // We will start a new group

      group++;
      groupsize=0; // New group size
      runSum = 0; // THe new running sum starts at zero

      // Special case
      // If cutwhenpassed = FALSE 
      // then we are filling up from below

      if (!cutwhenpassed) {
        groupVec[i] = group; // Update to new group
        // New starting value
        if (NumericVector::is_na(x[i])) {
	  runSum = 0;
        } else {	
  	  runSum = x[i] ; //* cutwhenpassed;
        }
        groupsize=1;
      }
    }

    if ((maxgroup>0) && (groupsize==maxgroup)) {
      group++;
      runSum = 0;
      groupsize=0;      
    }
    
  }
  return groupVec;
}
