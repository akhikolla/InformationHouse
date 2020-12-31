#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xd_ad(DoubleVector x,
                             int op,
                             double a) {
  R_xlen_t N = x.length();
  switch(op) {
  case OP_NE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != a) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == a) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= a) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > a) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= a) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < a) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xi_ai(IntegerVector x,
                             int op,
                             int a) {
  R_xlen_t N = x.length();
  switch(op) {
  case OP_NE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != a) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == a) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= a) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > a) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= a) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < a) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xi_ad(IntegerVector x,
                             int op,
                             double ad) {
  R_xlen_t N = x.length();
  if (N == 0) {
    return 0; // # nocov
  }
  int e = 0;
  if (!do_is_safe2int(ad)) {
    if (op == OP_EQ || op == OP_IN || ISNAN(ad)) {
      return 0;
    }
    if (op == OP_NE) {
      return N;
    }
    if (!R_finite(ad)) {
      if (ad == R_PosInf) {
        switch(op) {
        case OP_GE:
          return 0;
        case OP_LE:
          return N;
        case OP_GT:
          return 0;
        case OP_LT:
          return N;
        }
      }
      if (ad == R_NegInf) {
        switch(op) {
        case OP_GE:
          return N;
        case OP_LE:
          return 0;
        case OP_GT:
          return N;
        case OP_LT:
          return 0;
        }
      }
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
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != a) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == a) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= a) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > a) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= a) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < a) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}



// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xd_yd(DoubleVector x, int op, DoubleVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xi_yi(IntegerVector x, int op, IntegerVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xi_yd(IntegerVector x, int op, DoubleVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == NA_INTEGER || ISNAN(y[i])) continue;
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xd_yi(DoubleVector x, int op, IntegerVector y) {
  R_xlen_t N = x.length();
  if (N != y.length()) {
    stop("Lengths differ."); // # nocov
  }

  switch(op) {
  case OP_NE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_EQ:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_GT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LE:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= y[i]) {
        return i + 1;
      }
    }
    break;
  case OP_LT:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] < y[i]) {
        return i + 1;
      }
    }
    break;
  }

  return 0;
}



R_xlen_t do_which_last_xi_aii(IntegerVector x, int op, int a1, int a2) {
  R_xlen_t N = x.length();
  if (a2 == NA_INTEGER) {
    a2 = INT_MAX;
  }
  switch(op) {
  case OP_BW:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= a1 && x[i] <= a2) {
        return i + 1;
      }
    }
    break;
  case OP_BO:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > a1 && x[i] < a2) {
        return i + 1;
      }
    }
    break;
  case OP_BC:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= a1 || x[i] >= a2) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

R_xlen_t do_which_last_xi_add(IntegerVector x, int op, double a1, double a2) {
  R_xlen_t N = x.length();
  if (ISNAN(a1) || a1 < -2147483647) {
    a1 = R_NegInf;
  }
  if (ISNAN(a2) || a1 > 2147483647) {
    a2 = R_PosInf;
  }

  switch(op) {
  case OP_BW:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= a1 && x[i] <= a2) {
        return i + 1;
      }
    }
    break;
  case OP_BO:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > a1 && x[i] < a2) {
        return i + 1;
      }
    }
    break;
  case OP_BC:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= a1 || x[i] >= a2) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xd_add(DoubleVector x, int op, double a1, double a2) {
  R_xlen_t N = x.length();
  if (ISNAN(a1)) {
    a1 = R_NegInf;
  }
  if (ISNAN(a2)) {
    a2 = R_PosInf;
  }
  switch(op) {
  case OP_BW:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] >= a1 && x[i] <= a2) {
        return i + 1;
      }
    }
    break;
  case OP_BO:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] > a1 && x[i] < a2) {
        return i + 1;
      }
    }
    break;
  case OP_BC:
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] <= a1 || x[i] >= a2) {
        return i + 1;
      }
    }
    break;
  }
  return 0;
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last__(SEXP x, int op, SEXP y,
                         int ny,
                         int y1i,
                         int y2i,
                         double y1d,
                         double y2d) {
  if (TYPEOF(x) == INTSXP && TYPEOF(y) == INTSXP) {
    switch(ny) {
    case 1:
      return do_which_last_xi_ai(x, op, y1i);
    case 2:
      return do_which_last_xi_aii(x, op, y1i, y2i);
    default:
      return do_which_last_xi_yi(x, op, y);
    }
  }
  if (TYPEOF(x) == INTSXP && TYPEOF(y) == REALSXP) {
    switch(ny) {
    case 1:
      return do_which_last_xi_ad(x, op, y1d);
    case 2:
      return do_which_last_xi_add(x, op, y1d, y2d);
    default:
      return do_which_last_xi_yd(x, op, y);
    }
  }
  if (TYPEOF(x) == REALSXP && TYPEOF(y) == INTSXP) {
    switch(ny) {
    case 1:
      return do_which_last_xd_ad(x, op, y1d);
    case 2:
      return do_which_last_xd_add(x, op, y1d, y2d);
    default:
      return do_which_last_xd_yi(x, op, y);
    }
  }
  if (TYPEOF(x) == REALSXP && TYPEOF(y) == REALSXP) {
    switch(ny) {
    case 1:
      return do_which_last_xd_ad(x, op, y1d);
    case 2:
      return do_which_last_xd_add(x, op, y1d, y2d);
    default:
      return do_which_last_xd_yd(x, op, y);
    }
  }
  return 0; // # nocov
}




// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_in_lgl(LogicalVector x, bool anyNA_, bool any_, bool nall_) {
  R_xlen_t N = x.length();

  // anyNA -> has NA
  // any_ -> has TRUE
  // nall -> has FALSE

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
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == FALSE) {
        return i + 1;
      }
    }
  }

  // 2
  if (!anyNA_ & any_ & !nall_) {
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == TRUE) {
        return i + 1;
      }
    }
  }

  // 3
  if (!anyNA_ & any_ & nall_) {
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != NA_LOGICAL) {
        return i + 1;
      }
    }
  }

  // 4
  if (anyNA_ & !any_ & !nall_) {
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] == NA_LOGICAL) {
        return i + 1;
      }
    }
  }

  // 5
  if (anyNA_ & !any_ & nall_) {
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != TRUE) {
        return i + 1;
      }
    }
  }

  // 6
  if (anyNA_ & any_ & !nall_) {
    for (R_xlen_t i = N - 1; i >= 0; --i) {
      if (x[i] != FALSE) {
        return i + 1;
      }
    }
  }

  if (anyNA_ & any_ & nall_) {
    return N;
  }
  return 0;
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_last_xi_ini(IntegerVector x, IntegerVector y) {
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

    for (R_xlen_t i = Nx - 1; i >= 0; --i) {
      int xi = x[i];
      if (H.count(xi)) {
        return i + 1;
      }
    }
  } else {
    for (R_xlen_t i = Nx - 1; i >= 0; --i) {
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
R_xlen_t do_which_last_xi_ind(IntegerVector x, DoubleVector yd) {
  R_xlen_t Nx = x.length();
  if (yd.length() == 0) {
    return 0; // # nocov
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

    for (R_xlen_t i = Nx - 1; i >= 0; --i) {
      int xi = x[i];
      if (H.count(xi)) {
        return i + 1;
      }
    }
  } else {
    for (R_xlen_t i = Nx - 1; i >= 0; --i) {
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
R_xlen_t do_which_last_xd_ind(DoubleVector x, DoubleVector y, bool y_has_na) {
  R_xlen_t Nx = x.length();
  R_xlen_t Ny = y.length();
  if (Ny == 0) {
    return 0; // # nocov
  }

  if (Ny > 100) {
    // Above 100 elements it starts to become faster to use a hash map.
    // 100 fairly approximate but simply a round number
    std::unordered_set<double> H;
    for (R_xlen_t t = 0; t < Ny; ++t) {
      H.insert(y[t]);
    }

    for (R_xlen_t i = Nx - 1; i >= 0; --i) {
      double xi = x[i];
      if (y_has_na && ISNAN(xi)) {
        return i + 1;
      }
      if (H.count(xi)) {
        return i + 1;
      }
    }
  } else {
    for (R_xlen_t i = Nx - 1; i >= 0; --i) {
      double xi = x[i];
      if (y_has_na && ISNAN(xi)) {
        return i + 1;
      }
      for (R_xlen_t j = 0; j < Ny; ++j) {
        if (xi == y[j]) {
          return i + 1;
        }
      }
    }
  }
  return 0;
}

