#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix rle2_i(IntegerVector x, bool indices) {
  
  // Get length(x)
  int n = x.size();
  
  if (! indices) {
      
    // Loop through and record value, length for runs
    int rows = 1;
    for (int a = 1; a < n; ++a)
      if (x(a) != x(a - 1)) rows += 1;
    IntegerMatrix out(rows, 2);
    int presval = x(0);
    int prespos = 0;
    int index = -1;
    for (int a = 1; a < n; ++a) {
      int x_a = x[a];
      if (x_a != presval) {
        index += 1;
        out(index, 0) = presval;
        out(index, 1) = a - prespos;
        presval = x_a;
        prespos = a;
      }
    }
    index += 1;
    out(index, 0) = presval;
    out(index, 1) = n - prespos;
    return out;
    
  } 
  else {
      
    // Loop through and record value, length, start index, stop index for runs
    int rows = 1;
    for (int a = 1; a < n; ++a)
      if (x(a) != x(a - 1)) rows += 1;
    IntegerMatrix out(rows, 4);
    int presval = x(0);
    int prespos = 0;
    int index = -1;
    for (int a = 1; a < n; ++a) {
      int x_a = x(a);
      if (x_a != presval) {
        index += 1;
        out(index, 0) = presval;
        out(index, 1) = prespos + 1;
        out(index, 2) = a;
        out(index, 3) = a - prespos;
        presval = x_a;
        prespos = a;
      }
    }
    index += 1;
    out(index, 0) = presval;
    out(index, 1) = prespos + 1;
    out(index, 2) = n;
    out(index, 3) = n - prespos;
    return out;
        
  }
}

// [[Rcpp::export]]
NumericMatrix rle2_n(NumericVector x, bool indices) {
  
  // Get length(x)
  int n = x.size();
  
  if (! indices) {
    
    // Loop through and record value, length for runs
    int rows = 1;
    for (int a = 1; a < n; ++a)
      if (x(a) != x(a - 1)) rows += 1;
    NumericMatrix out(rows, 2);
    double presval = x(0);
    double prespos = 0;
    int index = -1;
    for (int a = 1; a < n; ++a) {
      double x_a = x[a];
      if (x_a != presval) {
        index += 1;
        out(index, 0) = presval;
        out(index, 1) = a - prespos;
        presval = x_a;
        prespos = a;
      }
    }
    index += 1;
    out(index, 0) = presval;
    out(index, 1) = n - prespos;
    return out;
      
  } 
  else {
    
    // Loop through and record value, start/stop index, and length for runs
    int rows = 1;
    for (int a = 1; a < n; ++a)
      if (x(a) != x(a - 1)) rows += 1;
    NumericMatrix out(rows, 4);
    double presval = x(0);
    double prespos = 0;
    int index = -1;
    for (int a = 1; a < n; ++a) {
      double x_a = x(a);
      if (x_a != presval) {
        index += 1;
        out(index, 0) = presval;
        out(index, 1) = prespos + 1;
        out(index, 2) = a;
        out(index, 3) = a - prespos;
        presval = x_a;
        prespos = a;
      }
    }
    index += 1;
    out(index, 0) = presval;
    out(index, 1) = prespos + 1;
    out(index, 2) = n;
    out(index, 3) = n - prespos;
    return out;
      
  }
}
