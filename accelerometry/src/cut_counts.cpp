#include <Rcpp.h>
using namespace Rcpp;

//' Cut Count Values into Intensity Ranges
//' 
//' Given a vector of accelerometer count values, classifies each count value 
//' into intensity level 1, 2, 3, 4, or 5 (typically representing sedentary, 
//' light, lifestyle, moderate, and vigorous).
//' 
//' 
//' @inheritParams artifacts
//' 
//' @param int_cuts Numeric vector with four cutpoints from which five intensity 
//' ranges are derived. For example, \code{int_cuts = c(100, 760, 2020, 5999)} 
//' creates: 0-99 = intensity 1; 100-759 = intensity level 2; 760-2019 = 
//' intensity 3; 2020-5998 = intensity 4; >= 5999 = intensity 5.
//' 
//' 
//' @return Integer vector.
//' 
//' 
//' @examples
//' # Load accelerometer data for first 5 participants in NHANES 2003-2004
//' data(unidata)
//' 
//' # Get data from ID number 21005
//' counts.part1 <- unidata[unidata[, "seqn"] == 21005, "paxinten"]
//' 
//' # Cut into 5 intensity levels and plot
//' intensity.part1 <- cut_counts(counts = counts.part1)
//' plot(intensity.part1)
//' 
//' 
//' @export
// [[Rcpp::export]]
IntegerVector cut_counts(IntegerVector counts,
                         IntegerVector int_cuts = IntegerVector::create(100, 760, 2020, 5999)) {
  
  // Get length(counts) and initialize output vector
  int n = counts.size();
  IntegerVector out(n, 5);
  
  // Extract cutpoints
  int cut1 = int_cuts(0);
  int cut2 = int_cuts(1);
  int cut3 = int_cuts(2);
  int cut4 = int_cuts(3);

  // Loop through and classify intensities
  for (int a = 0; a < n; ++a) {
    int counts_a = counts(a);
    if (counts_a < cut1) out(a) = 1;
    else if (counts_a < cut2) out(a) = 2;
    else if (counts_a < cut3) out(a) = 3;
    else if (counts_a < cut4) out(a) = 4;
  }
  return out;
  
}
