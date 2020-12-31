#include <Rcpp.h>
using namespace Rcpp;

//' Physical Activity Bout Detection
//' 
//' Identify bouts of physical activity based on a vector of accelerometer count 
//' values.
//' 
//' If \code{nci = FALSE}, the algorithm uses a moving window to go through 
//' every possible interval of length \code{bout_length} in \code{counts}. Any 
//' interval in which all counts are >= \code{tol_lower} and <= 
//' \code{tol_upper}, and no more than \code{tol} counts are less than 
//' \code{thresh_lower} or greater than \code{thresh_upper}, is classified as an 
//' activity bout.
//' 
//' If \code{nci = TRUE}, activity bouts are classified according to the 
//' algorithm used in the NCI's SAS programs. Briefly, this algorithm defines an 
//' activity bout as an interval of length \code{bout_length} that starts with a 
//' count value in \code{[thresh_lower, thresh_upper]} and has no more than 
//' \code{tol} counts outside of that range. If these criteria are met, the bout 
//' continues until there are \code{(tol + 1)} consecutive minutes outside of 
//' \code{[thresh_lower, thresh_upper]}. The parameters \code{tol_lower} and 
//' \code{tol_upper} are not used.
//' 
//' If the user allows for a tolerance (e.g. \code{tol = 2}) and does not use 
//' the NCI algorithm (i.e. \code{nci = FALSE}), specifying a non-zero value for 
//' \code{tol_lower} is highly recommended. Otherwise the algorithm will tend to 
//' classify minutes immediately before and after an activity bout as being part 
//' of the bout.
//' 
//' Specifying \code{thresh_lower} while using an arbitrarily large value for 
//' \code{thresh_upper} is generally recommended. Specifying both of these 
//' parameters can be overly restrictive in that the algorithm may miss bouts of 
//' activity in which counts are consistently high, but not exclusively in one 
//' intensity range.
//' 
//' 
//' @inheritParams weartime
//' 
//' @param weartime Integer vector with 1's for wear time minutes and 0's for 
//' non-wear time minutes.
//' 
//' @param bout_length Integer value specifying minimum length of an activity 
//' bout.
//' 
//' @param thresh_lower Integer value specifying lower bound for count values to 
//' be included for the intensity level. 
//' 
//' @param thresh_upper Integer value specifying upper bound for count values to 
//' be included for the intensity level.
//' 
//' @param tol Integer value specifying number of minutes with count values 
//' outside of [\code{thresh_lower}, \code{thresh_upper}] to allow during an 
//' activity bout.
//' 
//' @param tol_lower Integer value specifying lower cut-off for count values 
//' outside of intensity range during an activity bout.
//' 
//' @param tol_upper Integer value specifying upper cut-off for count values 
//' outside of intensity range during an activity bout.
//' 
//' @param days_distinct Logical value for whether to treat each day of data as 
//' distinct, i.e. identify non-wear time and activity bouts for day 1, then day 
//' 2, etc. If \code{FALSE}, algorithm is applied to full monitoring period 
//' continuously. If protocol has participants remove accelerometer for sleep, 
//' strongly recommend setting to \code{FALSE} to capture non-wear periods that 
//' start between 11 pm and midnight. Function assumes that first 1440 data 
//' points are day 1, next 1440 are day 2, and so on.
//' 
//' 
//' @return Integer vector with 1's for minutes that are part of an activity 
//' bout and 0's for minutes that are not.
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
//' wear.part1 <- weartime(counts = counts.part1)
//' 
//' # Identify moderate-to-vigorous activity bouts
//' mvpa.bouts <- bouts(counts = counts.part1, weartime = wear.part1, 
//'                     thresh_lower = 2020)
//' 
//' 
//' @export
// [[Rcpp::export]]
IntegerVector bouts(IntegerVector counts, 
                    Rcpp::Nullable<Rcpp::IntegerVector> weartime = R_NilValue, 
                    int bout_length = 10, 
                    int thresh_lower = 0, 
                    int thresh_upper = 100000, 
                    int tol = 0, 
                    int tol_lower = 0, 
                    int tol_upper = 100000, 
                    bool nci = false, 
                    bool days_distinct = false) {
  
  // Get length(counts) and initialize output vector
  int n = counts.size();
  IntegerVector out(n);
  
  // If weartime unspecified, make it a vector of 1's
  IntegerVector wear(n, 1);
  if (weartime.isNotNull()) wear = weartime;
  
  // Use appropriate version of algorithm given days_distinct, tol, and nci
  if (! days_distinct) {
    
    if (tol == 0) {
      
      int counter = 0;
      for (int a = 0; a < n; ++a) {
        if (wear[a] == 1 && counts[a] >= thresh_lower && 
            counts[a] <= thresh_upper) counter += 1;
        else {
          if (counter >= bout_length)
            for (int b = a - counter; b < a; ++b) out[b] = 1;
          counter = 0;
        }
        if (a == n - 1 && counter >= bout_length)
          for (int c = a - counter + 1; c < a + 1; ++c) out[c] = 1;
      }
      
    }            
    else if (tol > 0) {
      
      if (! nci) {
        
        IntegerVector status(n);
        for (int b = 0; b < n; ++b) {
          if (wear[b] == 1 && counts[b] >= thresh_lower && 
              counts[b] <= thresh_upper) status[b] = 0;
          else if (wear[b] == 1 && counts[b] >= tol_lower 
                     && counts[b] <= tol_upper) status[b] = 1;
          else status[b] = tol + 1;
        }
        int sum = 0;
        for (int c = 0; c < bout_length; ++c)
          sum += status[c];
        if (sum <= tol)
          for (int d = 0; d < bout_length; ++d) out[d] = 1;
        for (int e = bout_length; e < n; ++e) {
          sum -= status[e - bout_length] + status[e];
          if (sum <= tol)
            for (int f = e - bout_length + 1; f <= e ; ++f) out[f] = 1;
        }
        
      }            
      else {
        
        for (int a = 0; a < n - bout_length + 1; ++a) {
          if (counts[a] < thresh_lower || counts[a] > thresh_upper) continue;
          int tolcount = 0;
          int tolcount2 = 0;
          int last = 0;
          for (int b = a; b < a + bout_length; ++b) {
            if (counts[b] >= thresh_lower && counts[b] <= thresh_upper) {
              tolcount2 = 0;
              last = b;
            }
            else {
              tolcount += 1;
              tolcount2 += 1;
            }
            if (tolcount == tol + 1) break;
          }
          if (tolcount < tol + 1) {
            if (a + bout_length == n) {
              for (int c = a; c <= last; ++c) out[c] = 1;
              break;
            }
            for (int c = a + bout_length; c < n; ++c) {
              if (counts[c] >= thresh_lower && counts[c] <= thresh_upper) {
                tolcount2 = 0;
                last = c;
              }
              else tolcount2 += 1;
              if (tolcount2 == tol + 1 || c == n - 1) {
                for (int d = a; d <= last; ++d) out[d] = 1;
                break;
              }
            }
          }
        }
        
      }
    }
  }
  else if (days_distinct) {
    
    if (tol == 0) {
      
      int counter = 0;
      for (int a = 0; a < n; ++a) {
        if (wear[a] == 1 && counts[a] >= thresh_lower && 
            counts[a] <= thresh_upper) counter +=1;
        else {
          if (counter >= bout_length)
            for (int b = a - counter; b < a; ++b) out[b] = 1;
          counter = 0;
        }
        if ((a == n - 1 || (a + 1) % 1440 == 0) && counter >= bout_length)
          for (int c = a - counter + 1; c < a + 1; ++c) out[c] = 1;
        if ((a + 1) % 1440 == 0) counter = 0;
      }
      
    }              
    else if (tol > 0) {
      
      if (! nci) {
        
        IntegerVector status(n);
        for (int b = 0; b < n; ++b) {
          if (wear[b] == 1 && counts[b] >= thresh_lower && 
              counts[b] <= thresh_upper) status[b] = 0;
          else if (wear[b] == 1 && counts[b] >= tol_lower && 
                   counts[b] <= tol_upper) status[b] = 1;
          else status[b] = tol + 1;
        }
        int sum = 0;
        for (int c = 0; c < bout_length; ++c)
          sum += status[c];
        if (sum <= tol)
          for (int d = 0; d < bout_length; ++d) out[d] = 1;
        for (int e = bout_length; e < n; ++e) {
          sum -= status[e - bout_length] + status[e];
          if (sum <= tol && e % 1440 > bout_length - 2)
            for (int f = e - bout_length + 1; f <= e; ++f) out[f] = 1;
        }
        
      }
      else {
        
        for (int a = 0; a < n - bout_length + 1; ++a) {
          if (counts[a] < thresh_lower || counts[a] > thresh_upper) continue;
          int counter = 0;
          int tolcount = 0;
          int tolcount2 = 0;
          int last = 0;
          for (int b = a; b < a + bout_length; ++b) {
            if (counts[b] >= thresh_lower && counts[b] <= thresh_upper) {
              counter += 1;
              tolcount2 = 0;
              last = b;
            }
            else {
              counter += 1;
              tolcount += 1;
              tolcount2 += 1;
            }
            if (tolcount == tol + 1 || (b + 1) % 1440 == 0) break;
          }
          if (counter < bout_length || tolcount == tol + 1) continue;
          else {
            if ((a + bout_length) % 1440 == 0)
              for (int c = a; c <= last; ++c) out[c] = 1;
            else {
              for (int d = a + bout_length; d < n; ++d) {
                if (counts[d] >= thresh_lower && counts[d] <= thresh_upper) {
                  tolcount2 = 0;
                  last = d;
                }
                else tolcount2 += 1;
                if (tolcount2 == tol + 1 || d == n - 1 || (d + 1) % 1440 == 0) {
                  for (int e = a; e <= last; ++e) out[e] = 1;
                  break;
                }
              }
            }
          }
        }
        
      }
    }
  }
  
  // Return vector of 1's and 0's flagging activity bouts
  return(out);
}
