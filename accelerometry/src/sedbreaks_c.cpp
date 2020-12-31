#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int sedbreaks(IntegerVector counts, IntegerVector weartime, int thresh) {
  
  // Get length(counts) and initialize counter for sedentary breaks
  int n = counts.size();
  int sedbreaks = 0;
  
  // Loop through and count sedentary breaks
  for (int a = 0; a < n - 1; ++a) {
    int counts_a = counts[a];
    int counts_ap1 = counts[a + 1];
    int weartime_a = weartime[a];
    int weartime_ap1 = weartime[a + 1];
    if (weartime_a == 1 && weartime_ap1 == 1 && counts_a < thresh && 
        counts_ap1 >= thresh) sedbreaks += 1;
  }
  return(sedbreaks);
  
}

// [[Rcpp::export]]
IntegerVector sedbreaks_flags(IntegerVector counts, IntegerVector weartime, 
                              int thresh) {
  
  // Get length(counts) and initialize output vector
  int n = counts.size();
  IntegerVector out(n);
  
  // Loop through and flag sedentary breaks
  for (int a = 0; a < n - 1; ++a) {
    int weartime_a = weartime[a];
    int weartime_ap1 = weartime[a + 1];
    int counts_a = counts[a];
    int counts_ap1 = counts[a + 1];
    if (weartime_a == 1 && weartime_ap1 == 1 && counts_a < thresh && 
        counts_ap1 >= thresh) out[a + 1] = 1;
  }
  return(out);
  
}
