#include "cpphutils.h"

// [[Rcpp::export]]
LogicalVector do_or3(LogicalVector x, LogicalVector y, LogicalVector z) {
  R_xlen_t N = x.length();
  if (y.length() != N) {
    stop("y and x have different lengths.");
  }
  LogicalVector out(N);
  if (z.length() != N) {
    if (z.length() > 1) {
      stop("z has the wrong length");
    }
    if (z.length() == 0) {
      for (R_xlen_t i = 0; i < N; ++i) {
        out[i] = x[i] || y[i];
      }
    } else {
     if (z[0]) {
       for (R_xlen_t i = 0; i < N; ++i) {
         out[i] = true;
       }
     } else {
       for (R_xlen_t i = 0; i < N; ++i) {
         out[i] = x[i] || y[i];
       }
     }

    }
  } else {
    for (R_xlen_t i = 0; i < N; ++i) {
      out[i] = x[i] || y[i] || z[i];
    }
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector do_and3(LogicalVector x, LogicalVector y, LogicalVector z) {
  R_xlen_t N = x.length();
  if (y.length() != N) {
    stop("y and x have different lengths.");
  }
  LogicalVector out(N);
  if (z.length() != N) {
    if (z.length() > 1) {
      stop("z has the wrong length");
    }
    // if NULL -> fall back to binary &
    // if TRUE -> equivalent to binary &
    if (z.length() == 0 || z[0]) {
      for (R_xlen_t i = 0; i < N; ++i) {
        out[i] = x[i] && y[i];
      }
    } else {
      // z = false so all are false
      return out;


    }
  } else {
    for (R_xlen_t i = 0; i < N; ++i) {
      out[i] = x[i] && y[i] && z[i];
    }
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector na_and (LogicalVector x) {
  // NA & x
  R_xlen_t n = x.length();
  LogicalVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (x[i] != FALSE) {
      out[i] = NA_LOGICAL;
    }
  }
  return out;
}

// [[Rcpp::export]]
List do_which3(LogicalVector x, LogicalVector y, LogicalVector z,
               bool And = true,
               bool anyNAx = true,
               bool anyNAy = true,
               bool anyNAz = true) {
  R_xlen_t n = (x.length() > 1) ? x.length() : ((y.length() > 1) ? y.length() : z.length());
  const bool nx = x.length() == n;
  const bool ny = y.length() == n;
  const bool nz = z.length() == n;

  IntegerVector out(n);
  int j = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    int xi = nx ? x[i] : x[0];
    int yi = ny ? y[i] : y[0];
    int zi = nz ? z[i] : z[0];

    bool do_i = false;


    if (And) {
      if ((anyNAx ? ((xi != NA_LOGICAL) && xi) : xi) &&
          (anyNAy ? ((yi != NA_LOGICAL) && yi) : yi) &&
          (anyNAz ? ((zi != NA_LOGICAL) && zi) : zi)) {
        do_i = true;
      }
    } else {
      if ((anyNAx ? ((xi != NA_LOGICAL) && xi) : xi) ||
          (anyNAy ? ((yi != NA_LOGICAL) && yi) : yi) ||
          (anyNAz ? ((zi != NA_LOGICAL) && zi) : zi)) {
        do_i = true;
      }
    }
    if (do_i) {
      out[j] = i + 1;
      ++j;
    }
  }
  return List::create(j, out);
}

// [[Rcpp::export]]
IntegerVector do_which3_mem(LogicalVector x, LogicalVector y, LogicalVector z, bool And = true) {
  R_xlen_t n = (x.length() > 1) ? x.length() : ((y.length() > 1) ? y.length() : z.length());
  const bool nx = x.length() == n;
  const bool ny = y.length() == n;
  const bool nz = z.length() == n;

  R_xlen_t Count = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    int xi = nx ? x[i] : x[0];
    int yi = ny ? y[i] : y[0];
    int zi = nz ? z[i] : z[0];
    if (And) {
      if (xi && yi && zi) {
        ++Count;
      }
    } else {
      if (xi || yi || zi) {
        ++Count;
      }
    }
  }
  IntegerVector out(Count);
  R_xlen_t j = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    int xi = nx ? x[i] : x[0];
    int yi = ny ? y[i] : y[0];
    int zi = nz ? z[i] : z[0];
    if (And) {
      if (xi && yi && zi) {
        out[j] = i + 1;
        ++j;
      }
    } else {
      if (xi || yi || zi) {
        out[j] = i + 1;
        ++j;
      }
    }
  }
  return out;
}


