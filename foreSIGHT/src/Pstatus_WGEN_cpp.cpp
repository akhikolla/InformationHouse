#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Pstatus_WGEN_cpp(NumericVector parPwd, NumericVector parPdd, NumericVector RN, int ndays) {
  
  NumericVector drywet_TS(ndays);
  double pD = parPwd[0]/(1+parPwd[0]-parPdd[0]); // long-run probability of day 0 being dry
    
  if (RN[0] <= pD) {   // starting value based on long-run probability for day 0
    drywet_TS[0] = 0;
  } else {
    drywet_TS[0] = 1;
  }
  
  for(int i=1; i < ndays; ++i) {
    if(drywet_TS[i-1]==0){    // if the day before is dry
      if (RN[i] <= parPdd[i]) {  
        drywet_TS[i] = 0;
      } else {
        drywet_TS[i] = 1;
      }
    } else {                  // else, the day before was wet 
      if (RN[i] <= parPwd[i]) {  
        drywet_TS[i] = 0;
      } else {
        drywet_TS[i] = 1;
      }
    }
  }
  
  return drywet_TS;
}
