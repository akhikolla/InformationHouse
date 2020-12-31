#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first(LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 0; i < N; ++i) {
    if (x[i] == TRUE) {
      return ++i;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last(LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = N - 1; i >= 0; --i) {
    if (x[i] == TRUE) {
      return ++i;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_false (LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 0; i < N; ++i) {
    if (x[i] == FALSE) {
      return ++i;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_false (LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = N - 1; i >= 0; --i) {
    if (x[i] == FALSE) {
      return ++i;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_notTRUE(LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = 0; i < N; ++i) {
    if (x[i] != TRUE) {
      return ++i;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_notTRUE(LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = N - 1; i >= 0; --i) {
    if (x[i] != TRUE) {
      return ++i;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_notFALSE(LogicalVector x) {
  R_xlen_t N = x.size();
  for (R_xlen_t i = N - 1; i >= 0; --i) {
    if (x[i] != FALSE) {
      return ++i;
    }
  }
  return 0;
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_in_lgl(LogicalVector x, bool anyNA_, bool any_, bool nall_) {
  R_xlen_t N = x.length();

  // anyNA -> has NA
  // any_ -> has TRUE
  // nall -> not FALSE

  // anyNA  any_   nall_
  //     0     0       0  / 0
  //     0     0       1  / 1
  //     0     1       0  / 2
  //     0     1       1  / 3
  //     1     0       0  / 4
  //     1     0       1  / 5
  //     1     1       0  / 6
  //     1     1       1  / 7


  // 0
  if (!anyNA_ & !any_ & !nall_) {
    return 0;
  }

  // 1
  if (!anyNA_ & !any_ &  nall_) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == FALSE) {
        return i + 1;
      }
    }
  }

  // 2
  if (!anyNA_ & any_ & !nall_) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == TRUE) {
        return i + 1;
      }
    }
  }

  // 3
  if (!anyNA_ & any_ & nall_) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] != NA_LOGICAL) {
        return i + 1;
      }
    }
  }

  // 4
  if (anyNA_ & !any_ & !nall_) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == NA_LOGICAL) {
        return i + 1;
      }
    }
  }

  // 5
  if (anyNA_ & !any_ & nall_) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] != TRUE) {
        return i + 1;
      }
    }
  }

  // 6
  if (anyNA_ & any_ & !nall_) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] != FALSE) {
        return i + 1;
      }
    }
  }

  if (anyNA_ & any_ & nall_) {
    return 1;
  }
  return 0;
}

namespace hcwf {
template <int RTYPE>
R_xlen_t do_which_firstNA_(const Vector<RTYPE>& x) {
  R_xlen_t N = x.length();
  for (R_xlen_t i = 0; i < N; ++i) {
    if (Vector<RTYPE>::is_na(x[i])) {
      return i + 1;
    }
  }
  return 0;
}

template <int RTYPE>
R_xlen_t do_which_lastNA_(const Vector<RTYPE>& x) {
  R_xlen_t N = x.length();
  for (R_xlen_t i = N - 1; i >= 0; --i) {
    if (Vector<RTYPE>::is_na(x[i])) {
      return i + 1;
    }
  }
  return 0;
}

}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_firstNA(SEXP x) {
  switch(TYPEOF(x)) {
  case NILSXP:
    return 0;
  case LGLSXP:
    return hcwf::do_which_firstNA_(as<LogicalVector>(x));
  case INTSXP:
    return hcwf::do_which_firstNA_(as<IntegerVector>(x));
  case REALSXP:
    return hcwf::do_which_firstNA_(as<DoubleVector>(x));
  case STRSXP:
    return hcwf::do_which_firstNA_(as<CharacterVector>(x));
  case RAWSXP:
    return 0;
  }
  return 0; // # nocov
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_lastNA(SEXP x) {
  switch(TYPEOF(x)) {
  case NILSXP:
    return 0;
  case LGLSXP:
    return hcwf::do_which_lastNA_(as<LogicalVector>(x));
  case INTSXP:
    return hcwf::do_which_lastNA_(as<IntegerVector>(x));
  case REALSXP:
    return hcwf::do_which_lastNA_(as<DoubleVector>(x));
  case STRSXP:
    return hcwf::do_which_lastNA_(as<CharacterVector>(x));
  case RAWSXP:
    return 0;
  }
  return 0; // # nocov
}






// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_lgl_lgl_op(LogicalVector x, LogicalVector y, int op, bool reverse = false) {
  R_xlen_t N = x.length(), Ny = y.length();
  if (N == 0 || Ny == 0) {
    return 0;
  }
  if (op == OP_BO) {
    // even (F, T) can never occur
    return 0;
  }

  const bool len_eq = Ny == N;
  const bool len1 = Ny == 1;
  if (!len_eq && !len1 && op != OP_IN && op != OP_BW) {
    stop("Lengths differ."); // # nocov
  }
  if (op == OP_IN || op == OP_BW) {
    if (op == OP_BW && Ny != 2) {
      stop("%between% expects RHS to be a vector of length-2.");
    }
    bool hasNA = false;
    bool hasTRUE = false;
    bool hasFALSE = false;

    for (R_xlen_t j = 0; j < Ny; ++j) {
      hasNA = hasNA || y[j] == NA_LOGICAL;
      hasTRUE = hasTRUE || y[j] == TRUE;
      hasFALSE = hasFALSE || y[j] == FALSE;
    }

    if (hasNA && hasTRUE && hasFALSE) {
      if (reverse) {
        return N;
      } else {
        return 1;
      }
    }

    // Two values, for %between%, we need F,F F,T or T,T
    // otherwise will never occur so return 0 now
    if (op == OP_BW) {
      if (y[0] == TRUE && y[1] == FALSE) {
        return 0;
      }
      // c(NA, TRUE) => x <= TRUE
      y[0] = (y[0] == NA_LOGICAL) ? FALSE : y[0];
      y[1] = (y[1] == NA_LOGICAL) ? TRUE  : y[1];

      const bool onlyTRUE  = y[0] == TRUE;
      const bool onlyFALSE = y[1] == FALSE;

      // otherwise we just use the normal:
      for (R_xlen_t k = 0; k < N; ++k) {
        R_xlen_t i = reverse ? (N - k - 1) : k;
        if (onlyTRUE) {
          if (x[i] == TRUE) {
            return i + 1;
          }
          continue;
        }
        if (onlyFALSE) {
          if (x[i] == FALSE) {
            return i + 1;
          }
          continue;
        }
        if (x[i] != NA_LOGICAL) {
          return i + 1;
        }
      }
      return 0;
    }
    //
    for (R_xlen_t k = 0; k < N; ++k) {
      R_xlen_t i = reverse ? (N - k - 1) : k;
      if (hasNA && x[i] == NA_LOGICAL) {
        return i + 1;
      }
      if (hasTRUE && x[i] == TRUE) {
        return i + 1;
      }
      if (hasFALSE && x[i] == FALSE) {
        return i + 1;
      }
    }
    return 0;
  }

  for (R_xlen_t k = 0; k < N; ++k) {
    R_xlen_t i = reverse ? (N - k - 1) : k;
    int xi = x[i];
    int y1 = len1 ? y[0] : y[i];
    int y2 = (N == 2) ? y[1] : y[0];
    if (single_ox_x1_x2(xi, op, y1, y2)) {
      return i + 1;
    }
  }
  return 0;
}


R_xlen_t do_which_first_xd_ad(DoubleVector x,
                              int op,
                              double a) {
  R_xlen_t N = x.length();
  switch(op) {
  case OP_NE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] != a) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == a) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= a) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > a) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= a) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] < a) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

R_xlen_t do_which_first_xd_yd(DoubleVector x, int op, DoubleVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}


R_xlen_t do_which_first_xi_yi(IntegerVector x, int op, IntegerVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}

R_xlen_t do_which_first_xi_yd(IntegerVector x, int op, DoubleVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == NA_INTEGER || ISNAN(y[i])) continue;
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}

R_xlen_t do_which_first_xi_aii(IntegerVector x, int op, int a1, int a2) {
  R_xlen_t N = x.length();

  // a1 = NA_INTEGER is already INT MIN so no need to treat
  if (a2 == NA_INTEGER) {
    a2 = INT_MAX;
  }
  switch(op) {
  case OP_BW:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= a1 && x[i] <= a2) {
        return i + 1;
      }
    }
    break;
  case OP_BO:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > a1 && x[i] < a2) {
        return i + 1;
      }
    }
    break;
  case OP_BC:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= a1 || x[i] >= a2) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

R_xlen_t do_which_first_xi_add(IntegerVector x, int op, double a1, double a2) {
  R_xlen_t N = x.length();
  if (ISNAN(a1) || a1 < -2147483647) {
    a1 = R_NegInf;
  }
  if (ISNAN(a2) || a1 > 2147483647) {
    a2 = R_PosInf;
  }

  switch(op) {
  case OP_BW:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= a1 && x[i] <= a2) {
        return i + 1;
      }
    }
    break;
  case OP_BO:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > a1 && x[i] < a2) {
        return i + 1;
      }
    }
    break;
  case OP_BC:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= a1 || x[i] >= a2) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

R_xlen_t do_which_first_xd_add(DoubleVector x,
                              int op,
                              double a1,
                              double a2) {
  R_xlen_t N = x.length();
  if (ISNAN(a1)) {
    a1 = R_NegInf;
  }
  if (ISNAN(a2)) {
    a2 = R_PosInf;
  }
  switch(op) {
  case OP_BW:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] >= a1 && x[i] <= a2) {
        return i + 1;
      }
    }
    break;
  case OP_BO:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] > a1 && x[i] < a2) {
        return i + 1;
      }
    }
    break;
  case OP_BC:
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] <= a1 || x[i] >= a2) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}


R_xlen_t do_which_first_xi_ai(IntegerVector x,
                              int op,
                              int a) {
  R_xlen_t n = x.length();
  switch(op) {
  case OP_NE:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] != a) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] == a) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] >= a) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] > a) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] <= a) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] < a) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

R_xlen_t do_which_first_xi_ad(IntegerVector x,
                              int op,
                              double ad) {
  R_xlen_t n = x.length();
  if (n == 0) {
    return 0; // # nocov
  }
  int e = 0;
  if (!do_is_safe2int(ad)) {
    if (op == OP_EQ || op == OP_IN || ISNAN(ad)) {
      return 0;
    }
    if (op == OP_NE) {
      return 1;
    }
    if (!R_finite(ad)) {
      if (ad == R_PosInf) {
        switch(op) {
        case OP_GE:
          return 0;
        case OP_LE:
          return 1;
        case OP_GT:
          return 0;
        case OP_LT:
          return 1;
        }
      }
      if (ad == R_NegInf) {
        switch(op) {
        case OP_GE:
          return 1;
        case OP_LE:
          return 0;
        case OP_GT:
          return 1;
        case OP_LT:
          return 0;
        }
      }
      return 0;
    }
    // correct for integer conversion giving
    // incorrect results for >  >=  <  <=
    if (op == OP_GT) {
      // x >  2.5 => x > 2
      // x > -2.5 => x > -3
      e -= (ad < 0);
    }
    if (op == OP_GE) {
      // 3,4,5   >=  4.5  => 3,4,5 >= 5
      e += (ad > 0);
    }
    if (op == OP_LT) {
      // x <  2.5 => x < 3
      // x < -2.5 => x < -2
      e += (ad > 0);
    }
    if (op == OP_LE) {
      // x <= 4.5  => x <= 4
      // x <= -4.5 => x <= -4
      e -= (ad < 0);
    }
  }
  double ade = ((ad < 0) ? ceil(ad) : floor(ad)) + e;
  if (ade >= INT_MAX) {
    ade = INT_MAX;
  } 
  if (ade <= -INT_MAX) {
    ade = -INT_MAX;
  }
  const int a = (int)ade;

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] != a) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] == a) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] >= a) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] > a) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] <= a) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] < a) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first__(SEXP x, int op, SEXP y,
                          int ny,
                          int y1i,
                          int y2i,
                          double y1d,
                          double y2d) {

  if (TYPEOF(x) == INTSXP && TYPEOF(y) == INTSXP) {
    switch(ny) {
    case 1:
      return do_which_first_xi_ai(x, op, y1i);
    case 2:
      return do_which_first_xi_aii(x, op, y1i, y2i);
    default:
      return do_which_first_xi_yi(x, op, y);
    }
  }
  if (TYPEOF(x) == INTSXP && TYPEOF(y) == REALSXP) {
    switch(ny) {
    case 1:
      return do_which_first_xi_ad(x, op, y1d);
    case 2:
      return do_which_first_xi_add(x, op, y1d, y2d);
    default:
      return do_which_first_xi_yd(x, op, y);
    }
  }
  if (TYPEOF(x) == REALSXP && TYPEOF(y) == INTSXP) {
    switch(ny) {
    case 1:
      return do_which_first_xd_ad(x, op, y1d);
    case 2:
      return do_which_first_xd_add(x, op, y1d, y2d);
    default:
      return do_which_first_xd_yd(x, op, y);
    }
  }
  if (TYPEOF(x) == REALSXP && TYPEOF(y) == REALSXP) {
    switch(ny) {
    case 1:
      return do_which_first_xd_ad(x, op, y1d);
    case 2:
      return do_which_first_xd_add(x, op, y1d, y2d);
    default:
      return do_which_first_xd_yd(x, op, y);
    }
  }
  return 0; // # nocov
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_xi_ini(IntegerVector x, IntegerVector y) {
  R_xlen_t Nx = x.length();
  R_xlen_t Ny = y.length();
  if (Ny == 0) {
    return 0; // # nocov
  }

  if (Ny > 100) {
    // Above 100 elements it starts to become faster to use a hash map.
    // 100 fairly approximate but simply a round number
    std::unordered_set<int> H;
    for (R_xlen_t t = 0; t < Ny; ++t) {
      H.insert(y[t]);
    }

    for (R_xlen_t i = 0; i < Nx; ++i) {
      int xi = x[i];
      if (H.count(xi)) {
        return i + 1;
      }
    }
  } else {
    for (R_xlen_t i = 0; i < Nx; ++i) {
      int xi = x[i];
      for (R_xlen_t j = 0; j < Ny; ++j) {
        if (xi == y[j]) {
          return i + 1;
        }
      }
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
IntegerVector dblTable2int(DoubleVector table) {
  R_xlen_t N = table.length();
  std::set<int> o;
  // o.reserve(N);

  for (R_xlen_t i = 0; i < N; ++i) {
    double xd = table[i];
    switch(type_safe2int(xd)) {
    case 2:
      o.insert(NA_INTEGER);
      break;
    case 1:
      o.insert(static_cast<int>(xd));
      break;
    }
  }
  return wrap(o);
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_xi_ind(IntegerVector x, DoubleVector yd) {
  R_xlen_t Nx = x.length();
  if (yd.length() == 0) {
    return 0;
  }

  IntegerVector y = dblTable2int(yd);
  R_xlen_t Ny = y.length();

  if (Ny > 100) {
    // Above 100 elements it starts to become faster to use a hash map.
    // 100 fairly approximate but simply a round number
    std::unordered_set<int> H;
    for (R_xlen_t t = 0; t < Ny; ++t) {
      H.insert(y[t]);
    }

    for (R_xlen_t i = 0; i < Nx; ++i) {
      int xi = x[i];
      if (H.count(xi)) {
        return i + 1;
      }
    }
  } else {
    for (R_xlen_t i = 0; i < Nx; ++i) {
      int xi = x[i];
      for (R_xlen_t j = 0; j < Ny; ++j) {
        if (xi == y[j]) {
          return i + 1;
        }
      }
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_first_xd_ind(DoubleVector x, DoubleVector y, bool y_has_na) {
  R_xlen_t Nx = x.length();
  R_xlen_t Ny = y.length();
  if (Ny == 0) {
    return 0;
  }



  if (Ny > 100) {
    // Above 100 elements it starts to become faster to use a hash map.
    // 100 fairly approximate but simply a round number
    std::unordered_set<double> H;
    for (R_xlen_t t = 0; t < Ny; ++t) {
      H.insert(y[t]);
    }

    for (R_xlen_t i = 0; i < Nx; ++i) {
      double xi = x[i];
      if (y_has_na && ISNAN(xi)) {
        return i + 1;
      }
      if (H.count(xi)) {
        return i + 1;
      }
    }
  } else {
    for (R_xlen_t i = 0; i < Nx; ++i) {
      double xi = x[i];
      for (R_xlen_t j = 0; j < Ny; ++j) {
        if (y_has_na && ISNAN(xi)) {
          return i + 1;
        }
        if (xi == y[j]) {
          return i + 1;
        }
      }
    }
  }
  return 0;
}
