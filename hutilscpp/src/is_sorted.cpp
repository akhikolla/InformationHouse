#include "cpphutils.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
bool is_sorted_ascending_dbl(DoubleVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] < x[i - 1]) {
      return false;
    }
  }
  return true;
}

R_xlen_t do_isntSorted_ascending_dbl(DoubleVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] < x[i - 1]) {
      return i + 1;
    }
  }
  return 0;
}

// [[Rcpp::export]]
bool is_sorted_descending_dbl(DoubleVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] > x[i - 1]) {
      return false;
    }
  }
  return true;
}

R_xlen_t do_isntSorted_descending_dbl(DoubleVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i - 1] < x[i]) {
      return i + 1;
    }
  }
  return 0;
}

// [[Rcpp::export]]
bool is_sorted_ascending_int(IntegerVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] < x[i - 1]) {
      return false;
    }
  }
  return true;
}

R_xlen_t do_isntSorted_ascending_int(IntegerVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] < x[i - 1]) {
      return i + 1;
    }
  }
  return 0;
}

// [[Rcpp::export]]
bool is_sorted_descending_int(IntegerVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] > x[i - 1]) {
      return false;
    }
  }
  return true;
}

R_xlen_t do_isntSorted_descending_int(IntegerVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] > x[i - 1]) {
      return i + 1;
    }
  }
  return 0;
}

// [[Rcpp::export]]
bool is_sorted_int(IntegerVector x) {
  R_xlen_t n = x.length();
  if (n <= 2) {
    return true;
  }
  int xhead = x[0];
  int xtail = x[n - 1];

  if (xhead == xtail) {
    R_xlen_t i = 1;
    while (i < n && x[i] == xhead) {
      ++i;
    }
    return (i == n);
  }
  if (xhead < xtail) {
    return is_sorted_ascending_int(x);
  } else {
    return is_sorted_descending_int(x);
  }
}

// [[Rcpp::export]]
R_xlen_t do_isntSorted_int(IntegerVector x) {
  R_xlen_t n = x.length();
  if (n <= 2) {
    return 0;
  }
  int xhead = x[0];
  int xtail = x[n - 1];

  if (xhead == xtail) {
    R_xlen_t i = 1;
    while (i < n && x[i] == xhead) {
      ++i;
    }
    if (i == n) {
      return 0;
    }
    if (x[i] > xhead) {
      while (i < n && x[i] >= x[i - 1]) {
        ++i;
      }
      return i + 1;
    } else {
      while (i < n && x[i] <= x[i - 1]) {
        ++i;
      }
      return i + 1;
    }
  }
  if (xhead < xtail) {
    return do_isntSorted_ascending_int(x);
  } else {
    return do_isntSorted_descending_int(x);
  }
}

// [[Rcpp::export]]
bool is_sorted_dbl(DoubleVector x) {
  R_xlen_t n = x.length();
  if (n <= 2) {
    return true;
  }
  double xhead = x[0];
  double xtail = x[n - 1];

  if (xhead == xtail) {
    R_xlen_t i = 1;
    while (i < n && x[i] == xhead) {
      ++i;
    }
    return (i == n);
  }
  if (xhead < xtail) {
    return is_sorted_ascending_dbl(x);
  } else {
    return is_sorted_descending_dbl(x);
  }
}

// [[Rcpp::export]]
R_xlen_t do_isntSorted_dbl(DoubleVector x) {
  R_xlen_t n = x.length();
  if (n <= 2) {
    return 0;
  }
  double xhead = x[0];
  double xtail = x[n - 1];

  if (xhead == xtail) {
    R_xlen_t i = 1;
    while (i < n && x[i] == xhead) {
      ++i;
    }
    if (i == n) {
      return 0;
    }
    if (x[i] > xhead) {
      while (i < n && x[i] >= x[i - 1]) {
        ++i;
      }
      return i + 1;
    } else {
      while (i < n && x[i] <= x[i - 1]) {
        ++i;
      }
      return i + 1;
    }
  }
  if (xhead < xtail) {
    return do_isntSorted_ascending_dbl(x);
  } else {
    return do_isntSorted_descending_dbl(x);
  }
}
