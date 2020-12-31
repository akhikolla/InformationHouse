#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

typedef NumericVector::iterator iter;
typedef NumericVector NV;
typedef NumericMatrix NM;


// ----- roll_mean ---------------------------------------------------

/*
 * This is a simple rolling mean with an alignment parameter.
 */
 
 inline
double roll_mean(NV& x, const int& n, const int& ind, const int& alignCode) {

  // Calculate the average (including ind)
  double total = 0;
  for (int i=0; i<n; i++) {
    if (alignCode == -1) { // index at left edge of window
      total += x[ind + i];
    } else if (alignCode == 0) { // index at center of window
      int k = n / 2;
      total += x[ind -k + i];
    } else { // index at right edge of window
      total += x[ind - i];
    }
  }
  double out_ = total / n;

  return out_;
}

// NOTE:  We will return a vector of the same length as the incoming data 
// NOTE:  with NA's where there is not a full window's worth of data.

// [[Rcpp::export]]
NumericVector roll_mean_numeric_vector( NumericVector x, int n, int increment, int alignCode ) {

  int len = x.size();

  // Initialize the output vector with NA's
  NumericVector out(len, NA_REAL);

  // Get the start and endpoints
  int start = 0;
  int end = len;
  if (alignCode == -1) {
    // "left" aligned means the index is at the left edge of the window
    start = 0;
    end = len - (n - 1);
  } else if (alignCode == 0) {
    // "center" aligned
    start = n / 2;
    end = len - n / 2;
  } else {
    // "right" aligned means the index is at the right edge of the window
    start = n - 1;
    end = len;
  }

  // For the valid region, calculate the result
  for( int ind=start; ind<end; ind+=increment ) {
    out[ind] = roll_mean(x, n, ind, alignCode );
  }

  return out;

}


// ----- roll_median ---------------------------------------------------

/*
 * This is a simple rolling median that can be used to replace outliers
 * detected with the Hampel filter.
 */
 
 inline
double roll_median(NV& x, const int& n, const int& ind) {
  // Get the half-window width
  int k = n / 2;

  // Calculate the median value
  NumericVector tmp(n, NA_REAL);
  for( int i=0; i<n; i++ ) {
    tmp[i] = x[ind - k + i];
  }
  std::sort(tmp.begin(), tmp.end());
  double out_ = n % 2 ? tmp[n / 2] : (tmp[n / 2 - 1 ] + tmp[n / 2]) / 2;

  return out_;
}

// NOTE:  We will return a vector of the same length as the incoming data 
// NOTE:  with NA's in the half-window at either end.

// [[Rcpp::export]]
NumericVector roll_median_numeric_vector( NumericVector x, int n, int increment ) {

  int len = x.size();

  // Initialize the output vector with NA's
  NumericVector out(len, NA_REAL);

  // Get the half-window width
  int k = n / 2;

  // For the valid region, calculate the result
  for( int ind=k; ind<len-k; ind+=increment ) {
    out[ind] = roll_median(x, n, ind );
  }

  return out;

}


// ----- roll_hampel ---------------------------------------------------

/*
 * Here is the R code for the 'hampel' function from the 'pracma' package:
 *
 * NOTE:  R indices start with '1' whereas cpp start with '0'.
 * NOTE:  R's a:b syntax means (i=1; i<b+1; i++)
 *
function (x, k, t0 = 3) 
{
    n <- length(x)
    y <- x
    ind <- c()
    L <- 1.4826
    for (i in (k + 1):(n - k)) {
        x0 <- median(x[(i - k):(i + k)])
        S0 <- L * median(abs(x[(i - k):(i + k)] - x0))
        if (abs(x[i] - x0) > t0 * S0) {
            y[i] <- x0
            ind <- c(ind, i)
        }
    }
    list(y = y, ind = ind)
}
*/

 inline
double roll_hampel(NV& x, const int& n, const int& ind) {
  const double L = 1.4826;

  // Get the half-window width
  int k = n / 2;

  // Calculate x0 = median(x[window])
  NumericVector tmp(n, NA_REAL);
  for( int i=0; i<n; i++ ) {
    tmp[i] = x[ind - k + i];
  }
  std::sort(tmp.begin(), tmp.end());
  double x0 = n % 2 ? tmp[n / 2] : (tmp[n / 2 - 1 ] + tmp[n / 2]) / 2;

  // Calculate the absolute difference between each point and the median
  NumericVector absMinusMedian(n, NA_REAL);
  for ( int i=0; i<n; i++) {
    absMinusMedian[i] = fabs(x[ind - k + i] - x0);
  }

  // Calculate the median of these values
  std::sort(absMinusMedian.begin(), absMinusMedian.end());
  double medianAbsMinusMedian =  n % 2 ? absMinusMedian[n / 2] : (absMinusMedian[n / 2 - 1] + absMinusMedian[n / 2]) / 2;

  // Calculate return value
  double S0 = L * medianAbsMinusMedian;
  double out_ = fabs(x[ind] - x0) / S0;

  return out_;
}

// NOTE:  We will return a vector of the same length as the incoming data 
// NOTE:  with NA's in the half-window at either end.

// [[Rcpp::export]]
NumericVector roll_hampel_numeric_vector( NumericVector x, int n, int increment ) {

  int len = x.size();

  // Initialize the output vector with NA's
  NumericVector out(len, NA_REAL);

  // Get the half-window width
  int k = n / 2;

  // For the valid region, calculate the result
  for( int ind=k; ind<len-k; ind+=increment ) {
    out[ind] = roll_hampel(x, n, ind );
  }

  return out;

}


// ----- roll_sd ---------------------------------------------------------------

/*
 * This is a simple rolling standard deviation with an alignment parameter.
 */
 
 inline
double roll_sd(NV& x, const int& n, const int& ind, const int& alignCode) {

  // Calculate the mean (including ind)
  double total = 0;
  for (int i=0; i<n; i++) {
    if (alignCode == -1) { // index at left edge of window
      total += x[ind + i];
    } else if (alignCode == 0) { // index at center of window
      int k = n / 2;
      total += x[ind -k + i];
    } else { // index at right edge of window
      total += x[ind - i];
    }
  }
  double mean = total / n;
  
  double out_ = 0;
  for( int i=0; i<n; i++ ) {
    if (alignCode == -1) { // index at left edge of window
      out_ += (x[ind + i]-mean)*(x[ind + i]-mean);
    } else if (alignCode == 0) { // index at center of window
      int k = n / 2;
      out_ += (x[ind - k + i]-mean)*(x[ind - k + i]-mean);
    } else { // index at right edge of window
      out_ += (x[ind-i]-mean)*(x[ind-i]-mean);
    }
  }
  out_ = sqrt( out_/(n-1) );

  return out_;
}

// NOTE:  We will return a vector of the same length as the incoming data 
// NOTE:  with NA's where there is not a full window's worth of data.

// [[Rcpp::export]]
NumericVector roll_sd_numeric_vector( NumericVector x, int n, int increment, int alignCode ) {

  int len = x.size();

  // Initialize the output vector with NA's
  NumericVector out(len, NA_REAL);

  // Get the start and endpoints
  int start = 0;
  int end = len;
  if (alignCode == -1) {
    // "left" aligned means the index is at the left edge of the window
    start = 0;
    end = len - (n - 1);
  } else if (alignCode == 0) {
    // "center" aligned
    start = n / 2;
    end = len - n / 2;
  } else {
    // "right" aligned means the index is at the right edge of the window
    start = n - 1;
    end = len;
  }

  // For the valid region, calculate the result
  for( int ind=start; ind<end; ind+=increment ) {
    out[ind] = roll_sd(x, n, ind, alignCode );
  }

  return out;

}


// ----- roll_stalta ---------------------------------------------------

/*
 * This is a simple ratio of two rolling means that is used to detect
 * first breaks.
 */
 
 inline
double roll_stalta(NV& x, const int& n_sta, const int& n_lta, const int& ind) {

  // Calculate the STA, aligned so that the window is right of ind (including ind)
  // NOTE:  Pay attention to for loop components!
  double total = 0;
  for (int i=0; i<n_sta; i++) {
    total += x[ind + i];
  }
  double sta = total / n_sta;

  // Calculate the LTA, aligned so that the window is left of ind (including ind)
  // NOTE:  Pay attention to for loop components!
  total = 0;
  for (int i=0; i<n_lta; i++) {
    total += x[ind - i];
  }
  double lta = total / n_lta;

  // Calcualte return value
  double out_ = sta/lta;
  return out_;
}

// NOTE:  We will return a vector of the same length as the incoming data 
// NOTE:  with NA's in the LTA window length at the left end and in the 
// NOTE:  STA window length at the right end.

// [[Rcpp::export]]
NumericVector roll_stalta_numeric_vector( NumericVector x, int n_sta, int n_lta, int increment ) {

  int len = x.size();

  // Initialize the output vector with NA's
  NumericVector out(len, NA_REAL);

  // For the valid region, calculate the result
  for( int ind=n_lta; ind<len-n_sta; ind+=increment ) {
    out[ind] = roll_stalta(x, n_sta, n_lta, ind );
  }

  return out;

}

// ----- roll_range --------------------------------------------------

 inline
double roll_range(NV& x, const int& n, const int& ind, const int& alignCode) {
   double xMin = x[ind];
   double xMax = x[ind];

   for (int i=0; i<n; i++) {
     if (alignCode == -1) { // index at left edge of window
        if (std::isnan(xMin)) {
            xMin = x[ind +i];
        }
        if (x[ind +i] < xMin ) {
             xMin = x[ind +i];
        }
        if (std::isnan(xMax)) {
             xMax = x[ind +i];
        }
        if (x[ind +i] > xMax ) {
             xMax = x[ind +i];
        }
     } else if (alignCode == 0) { // index at center of window
        int k = n / 2;
        if (std::isnan(xMin)) {
            xMin = x[ind -k + i];
        }
        if (x[ind -k + i] < xMin) {
            xMin = x[ind -k + i];
        }
        if (std::isnan(xMax)) {
             xMax = x[ind -k + i];
        }
        if (x[ind -k + i] > xMax) {
            xMax = x[ind -k + i];
        }
     } else { // index at right edge of window
        if (std::isnan(xMin)) {
            xMin = x[ind - i];
        }
        if (x[ind - i] < xMin) {
             xMin = x[ind - i];
        }
        if (std::isnan(xMax)) {
             xMax = x[ind -i];
        }
        if (x[ind - i] > xMax) {
             xMax = x[ind - i];
        }
     }
   }

   double out_ = xMax - xMin;
   return out_;
}

// [[Rcpp::export]]
NumericVector roll_range_numeric_vector( NumericVector x, int n, int increment, int alignCode ) {

  int len = x.size();

  // Initialize the output vector with NA's
  NumericVector out(len, NA_REAL);

  // Get the start and endpoints
  int start = 0;
  int end = len;
  if (alignCode == -1) {
    // "left" aligned means the index is at the left edge of the window
    start = 0;
    end = len - (n - 1);
  } else if (alignCode == 0) {
    // "center" aligned
    start = n / 2;
    end = len - n / 2;
  } else {
    // "right" aligned means the index is at the right edge of the window
    start = n - 1;
    end = len;
  }

  // For the valid region, calculate the result
  for( int ind=start; ind<end; ind+=increment ) {
    out[ind] = roll_range(x, n, ind, alignCode );
  }

  return out;

}

