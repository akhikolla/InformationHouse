#include <Rcpp.h>
using namespace Rcpp;

//' Accelerometer Artifact Correction
//'
//' Corrects abnormally high count values in accelerometer data by replacing 
//' such values with the average of neighboring count values. Returns integer 
//' vector despite the average calculation often producing a decimal; this 
//' follows the convention used in the NCI's SAS programs 
//' (\url{http://riskfactor.cancer.gov/tools/nhanes_pam}). 
//' 
//' 
//' @param counts Integer vector with accelerometer count values.
//' @param thresh Integer value specifying the smallest count value that should 
//' be considered an artifact.
//' @param counts_classify Integer vector with accelerometer count values to 
//' base artifact classification on, but not to adjust. Mainly included for 
//' triaxial data, where you might want to define artifacts based on 
//' vertical-axis counts but then actually adjust the triaxial sum or vector 
//' magnitude counts.
//'
//'
//' @return Integer vector equivalent to \code{counts} except where artifacts 
//' were adjusted.
//' 
//' 
//' @references 
//' National Cancer Institute. Risk factor monitoring and methods: SAS programs 
//' for analyzing NHANES 2003-2004 accelerometer data. Available at: 
//' \url{http://riskfactor.cancer.gov/tools/nhanes_pam}. Accessed Aug. 19, 2018.
//'
//'
//' @examples
//' # Load accelerometer data for first 5 participants in NHANES 2003-2004
//' data(unidata)
//' 
//' # Get data from ID number 21007
//' counts.part3 <- unidata[unidata[, "seqn"] == 21007, "paxinten"]
//' 
//' # Replace counts > 10,000 with average of neighboring values
//' counts.part3.corrected <- artifacts(counts = counts.part3, thresh = 10000)
//'
//'@export
// [[Rcpp::export]]
IntegerVector artifacts(IntegerVector counts,
                        int thresh, 
                        Rcpp::Nullable<Rcpp::IntegerVector> counts_classify = R_NilValue) {
  
  IntegerVector counts_class = counts;
  if (counts_classify.isNotNull()) counts_class = counts_classify;
  
  int n = counts.size();
  IntegerVector out(n);

  int before = -1;
  int after = -1;
  if (counts_class[0] >= thresh) {
    for (int a = 1; a < n; ++a) {
      if (counts_class[a] < thresh) {
        out[0] = counts[a];
        break;
      }
    }
  }
  else out[0] = counts[0];
  for (int b = 1; b < n; ++b) {
    before = -1;
    after = -1;
    if (counts_class[b] >= thresh) {
      for (int c = b - 1; c >= 0; --c) {
        if (counts_class[c] < thresh) {
          before = counts[c];
          break;
        }
      }
      for (int d = b + 1; d < n; ++d) {
        if (counts_class[d] < thresh) {
          after = counts[d];
          break;
        }
      }
      if (before > -1 && after > -1) {
        if (std::fmod(before + after, 2) == 0) out[b] = (before + after) / 2;
        else out[b] = (before + after + 1) / 2;
      }
      else if (before == -1) out[b] = after;
      else if (after == -1) out[b] = before;
    }
    else out[b] = counts[b];
  }
  return(out);
}
