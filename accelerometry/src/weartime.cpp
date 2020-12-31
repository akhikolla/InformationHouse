#include <Rcpp.h>
using namespace Rcpp;

//' Wear Time Classification
//' 
//' Classifies wear time vs. non-wear time based on a vector of accelerometer 
//' count values.
//' 
//' If \code{nci = FALSE}, the algorithm uses a moving window to go through 
//' every possible interval of length \code{window} in \code{counts}. Any 
//' interval in which no more than \code{tol} counts are non-zero, and those 
//' are still < \code{tol.upper}, is classified as non-wear time. 
//' 
//' If \code{nci = TRUE}, non-wear time is classified according to the algorithm 
//' used in the NCI's SAS programs. Briefly, this algorithm defines a non-wear 
//' period as an interval of length \code{window} that starts with a count value 
//' of 0, does not contain any periods with \code{(tol + 1)} consecutive 
//' non-zero count values, and does not contain any counts > \code{tol.upper}. 
//' If these criteria are met, the non-wear period continues until there are 
//' \code{(tol + 1)} consecutive non-zero count values or a single count value > 
//' \code{tol.upper}.
//' 
//' 
//' @param counts Integer vector with accelerometer count values.
//'
//' @param window Integer value specifying minimum length of a non-wear 
//' period.
//' 
//' @param tol Integer value specifying tolerance for non-wear algorithm, i.e. 
//' number of seconds/minutes with non-zero counts allowed during a non-wear 
//' interval.
//' 
//' @param tol_upper Integer value specifying maximum count value for a 
//' second/minute with non-zero counts during a non-wear interval.
//' 
//' @param nci Logical value for whether to use algorithm from NCI's SAS 
//' programs. See \bold{Details}.
//' 
//' @param days_distinct Logical value for whether to treat each day of data as 
//' distinct, as opposed to analyzing the entire monitoring period as one 
//' continuous segment. For minute-to-minute counts, strongly recommend setting 
//' to \code{FALSE} to correctly classify time near midnight.
//' 
//' @param units_day Integer value specifying how many data point are in a day. 
//' Typically either 1440 or 86400 depending on whether count values are 
//' minute-to-minute or second-to-second.
//' 
//' 
//' @return Integer vector with 1's for valid wear time and 0's for non-wear 
//' time.
//' 
//' 
//' @references 
//' National Cancer Institute. Risk factor monitoring and methods: SAS programs 
//' for analyzing NHANES 2003-2004 accelerometer data. Available at: 
//' \url{http://riskfactor.cancer.gov/tools/nhanes_pam}. Accessed Aug. 19, 2018.
//' 
//' Acknowledgment: This material is based upon work supported by the National 
//' Science Foundation Graduate Research Fellowship under Grant No. DGE-0940903.
//' 
//' 
//' @examples
//' # Load accelerometer data for first 5 participants in NHANES 2003-2004
//' data(unidata)
//' 
//' # Get data from ID number 21005
//' counts.part1 <- unidata[unidata[, "seqn"] == 21005, "paxinten"]
//' 
//' # Identify periods of valid wear time
//' weartime.flag <- weartime(counts = counts.part1)
//' 
//' 
//' @export
// [[Rcpp::export]]
IntegerVector weartime(IntegerVector counts, int window = 60, int tol = 0, 
                       int tol_upper = 99, bool nci = false, 
                       bool days_distinct = false, int units_day = 1440) {
  
  // Get length(counts) and initialize output vector starting with all 1's
  int n = counts.size();
  IntegerVector out(n, 1);
  
  // Use appropriate version of algorithm given days_distinct, tol, and nci
  if (! days_distinct) {
    
    if (tol == 0) {
      
      int zeros = 0;
      for (int b = 0; b < n; ++b) {
        if (counts[b] == 0) zeros +=1;
        else {
          if (zeros >= window)
            for (int c = b - zeros; c < b; ++c) out[c] = 0;
          zeros = 0;
        }
        if (b == n - 1 && zeros >= window)
          for (int d = b - zeros + 1; d < b + 1; ++d) out[d] = 0;
      }
      
    }
    else if (tol > 0) {
      
      if (! nci) {
        
        IntegerVector status(n);
        for (int b = 0; b < n; ++b) {
          int counts_b = counts[b];
          if (counts_b == 0) status[b] = 0;
          else if (counts_b <= tol_upper) status[b] = 1;
          else if (counts_b > tol_upper) status[b] = tol + 1;
        }
        int sum = 0;
        for (int c = 0; c < window; ++c) 
          sum += status[c];
        if (sum <= tol)
          for (int d = 0; d < window; ++d) out[d] = 0;
        for (int e = window; e < n; ++e) {
          sum = sum - status[e - window] + status[e];
          if (sum <= tol)
            for (int f = e - window + 1; f <= e; ++f) out[f] = 0;
        }
        
      }            
      else if (nci) {
        
        int zeros = 0;
        int tolcount = 0;
        int flag = 0;
        for (int b = 0; b < n; ++b) {
          int counts_b = counts[b];
          if (zeros == 0 && counts_b != 0) continue;
          if (counts_b == 0) {
            zeros += 1;
            tolcount = 0;
          }
          else if (counts_b > 0 && counts_b <= tol_upper) {
            zeros += 1;
            tolcount += 1;
          }
          else if (counts[b]>tol_upper) {
            zeros += 1;
            tolcount += 1;
            flag = 1;
          }
          if (tolcount > tol || flag == 1 || b == n - 1) {
            if (zeros - tolcount >= window)
              for (int c = b - zeros + 1; c < b - tolcount + 1; ++c) out[c] = 0;
            zeros = 0;
            tolcount = 0;
            flag = 0;
          }
        }
        
      }
    }
  }
  else {
    
    if (tol == 0) {
      
      int zeros = 0;
      for (int b = 0; b < n; ++b) {
        if (counts[b] == 0) zeros +=1;
        else {
          if (zeros >= window)
            for (int c = b - zeros; c < b; ++c) out[c] = 0;
          zeros = 0;
        }
        if ((b == n-1 || (b + 1) % units_day == 0) && zeros >= window)
          for (int d = b - zeros + 1; d < b + 1; ++d) out[d] = 0;
        if ((b + 1) % units_day == 0) zeros = 0;
      }
      
    }                
    else if (tol > 0) {
      
      if (! nci) {
        
        IntegerVector status(n);
        for (int b = 0; b < n; ++b) {
          int counts_b = counts[b];
          if (counts_b == 0) status[b] = 0;
          else if (counts_b <= tol_upper) status[b] = 1;
          else if (counts_b > tol_upper) status[b] = tol + 1;
        }
        int sum = 0;
        for (int c = 0; c < window; ++c)
          sum += status[c];
        if (sum <= tol)
          for (int d = 0; d < window; ++d) out[d] = 0;
        for (int e = window; e < n; ++e) {
          sum = sum - status[e - window] + status[e];
          if (sum <= tol && e % units_day > window - 2)
            for (int f = e - window + 1; f <= e; ++f) out[f] = 0;
        }
        
      }                
      else if (nci) {
        
        int zeros = 0;
        int tolcount = 0;
        int flag = 0;
        for (int b = 0; b < n; ++b) {
          int counts_b = counts[b];
          if (zeros == 0 && counts_b != 0) continue;
          if (counts_b == 0) {
            zeros += 1;
            tolcount = 0;
          }
          else if (counts_b > 0 && counts_b <= tol_upper) {
            zeros += 1;
            tolcount += 1;
          }
          else if (counts_b > tol_upper) {
            zeros += 1;
            tolcount += 1;
            flag = 1;
          }
          if (tolcount > tol || flag == 1 || b == n - 1 || (b + 1) % units_day == 0) {
            if (zeros-tolcount>=window)
              for (int c = b-zeros+1; c < b-tolcount+1; ++c) out[c] = 0;
            zeros = 0;
            tolcount = 0;
            flag = 0;
          }
        }
        
      }
    }
  }
  
  // Return output vector
  return(out);
}
