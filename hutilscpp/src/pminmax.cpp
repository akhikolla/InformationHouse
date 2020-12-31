#include "cpphutils.h"



void message(std::string txt) {
  Rcpp::Function msg("message");
  msg(std::string("") + txt);
}






// [[Rcpp::export(rng = false)]]
R_xlen_t do_firstNonNegativeRadix_int(IntegerVector x,
                                      R_xlen_t mini = 0,
                                      R_xlen_t maxi = -1,
                                      bool desc = false,
                                      int depth = 0) {

  R_xlen_t xsize = x.length();
  if (maxi < 0 || maxi > xsize) {
    if (xsize < 1) {
      return xsize; // # nocov
    }
    maxi = xsize;
  }
  if (mini < 0) {
    mini = 0;
  }
  int lastx = x[maxi - 1];

  if (desc) {
    if (x[mini] < 0 || lastx > 0) {
      return mini;
    }
  } else {
    if (x[mini] > 0 || lastx < 0) {
      return mini;
    }
  }

  if (mini > maxi - 1024 || depth > 31) {
    for (R_xlen_t i = mini; i < maxi; ++i) {
      if (desc) {
        if (x[i] <= 0) {
          return i;
        }
      } else {
        if (x[i] >= 0) {
          return i;
        }
      }
    }
    return maxi; // # nocov
  }
  R_xlen_t medi = mini + (maxi - mini) / 2;
  bool lhs = (x[medi] < 0) ? desc : !desc;
  R_xlen_t left = lhs ? mini : medi - 1;
  R_xlen_t right = lhs ? medi + 2 : maxi;
  return do_firstNonNegativeRadix_int(x, left, right, desc, depth + 1);
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_firstNonNegativeRadix_dbl(DoubleVector x,
                                      R_xlen_t mini = 0,
                                      R_xlen_t maxi = -1,
                                      bool desc = false,
                                      int depth = 0) {

  R_xlen_t xsize = x.length();
  if (maxi < 0 || maxi > xsize) {
    if (xsize < 1) {
      return 0; // # nocov
    }
    maxi = xsize;
  }
  if (mini < 0) {
    mini = 0;
  }
  double lastx = x[maxi - 1];


  if (desc) {
    if (x[mini] < 0 || lastx > 0) {
      return mini;
    }
  } else {
    if (x[mini] > 0 || lastx < 0) {
      return mini;
    }
  }

  if (mini > maxi - 1024 || depth > 31) {
    // showValuex("depth = ", depth);
    for (R_xlen_t i = mini; i < maxi; ++i) {
      if (desc) {
        if (x[i] <= 0) {
          return i;
        }
      } else {
        if (x[i] >= 0) {
          return i;
        }
      }
    }
    return maxi; // # nocov
  }
  R_xlen_t medi = mini + (maxi - mini) / 2;
  bool lhs = (x[medi] < 0) ? desc : !desc;
  R_xlen_t left = lhs ? mini : medi - 1;
  R_xlen_t right = lhs ? medi + 1 : maxi;
  return do_firstNonNegativeRadix_dbl(x, left, right, desc, depth + 1);
}

// [[Rcpp::export(rng = false)]]
DoubleVector do_pmax0_radix_sorted_dbl(DoubleVector x,
                                       bool in_place = false,
                                       int nThread = 1) {
  R_xlen_t N = x.size();
  if (N == 0) {
    return x;
  }
  if (N == 1) {
    return (x[0] >= 0) ? x : DoubleVector(1);
  }
  bool x0_positive = x[0] > 0;
  bool xn_positive = x[N - 1] > 0;
  if (x0_positive && xn_positive) {
    return x;
  }
  if (!x0_positive && !xn_positive) {
    return DoubleVector(N);
  }
  bool desc = x[0] > 0;
  R_xlen_t root = do_firstNonNegativeRadix_dbl(x, 0, N, desc);
  if (in_place) {
    DoubleVector out(x);
    if (desc) {
      for (R_xlen_t i = root; i < N; ++i) {
        out[i] = 0;
      }
    } else {
      for (R_xlen_t i = 0; i < root; ++i) {
        out[i] = 0;
      }
    }
    return out;
  } else {
    DoubleVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (R_xlen_t i = 0; i < N; ++i) {
      bool zero = desc ? i >= root : i < root;
      out[i] = zero ? 0 : x[i];
    }
    return out;
  }
}

// [[Rcpp::export(rng = false)]]
DoubleVector do_pmin0_radix_sorted_dbl(DoubleVector x,
                                       bool in_place = false,
                                       int nThread = 1) {
  R_xlen_t N = x.size();
  bool x0_positive = x[0] > 0;
  bool xn_positive = x[N - 1] > 0;
  if (!x0_positive && !xn_positive) {
    return x;
  }
  if (x0_positive && xn_positive) {
    return DoubleVector(N);
  }
  const bool desc = x[0] > 0;
  R_xlen_t root = do_firstNonNegativeRadix_dbl(x, 0, N, desc);
  if (in_place) {
    DoubleVector out(x);
    R_xlen_t a = desc ? 0 : root;
    R_xlen_t b = desc ? root : N;
    for (R_xlen_t i = a; i < b; ++i) {
      out[i] = 0;
    }
    return out;
  } else {
    DoubleVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (R_xlen_t i = 0; i < N; ++i) {
      bool zero = desc ? i < root : i >= root;
      out[i] = zero ? 0 : x[i];
    }
    return out;
  }
}


// [[Rcpp::export(rng = false)]]
IntegerVector do_pmax0_radix_sorted_int(IntegerVector x,
                                        bool in_place = false,
                                        int nThread = 1) {
  R_xlen_t N = x.size();
  bool x0_positive = x[0] > 0;
  bool xn_positive = x[N - 1] > 0;
  if (x0_positive && xn_positive) {
    return x;
  }
  if (!x0_positive && !xn_positive) {
    return IntegerVector(N);
  }
  bool desc = x[0] > 0;
  R_xlen_t root = do_firstNonNegativeRadix_int(x, 0, N, desc);
  if (in_place) {
    IntegerVector out(x);
    if (desc) {
      for (R_xlen_t i = root; i < N; ++i) {
        out[i] = 0;
      }
    } else {
      for (R_xlen_t i = 0; i < root; ++i) {
        out[i] = 0;
      }
    }
    return out;
  } else {
    IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (R_xlen_t i = 0; i < N; ++i) {
      bool zero = desc ? i >= root : i < root;
      out[i] = zero ? 0 : x[i];
    }
    return out;
  }
}

// [[Rcpp::export(rng = false)]]
IntegerVector do_pmin0_radix_sorted_int(IntegerVector x,
                                        bool in_place = false,
                                        int nThread = 1) {
  R_xlen_t N = x.size();
  bool x0_positive = x[0] > 0;
  bool xn_positive = x[N - 1] > 0;
  if (!x0_positive && !xn_positive) {
    return x;
  }
  if (x0_positive && xn_positive) {
    return IntegerVector(N);
  }
  const bool desc = x[0] > 0;
  R_xlen_t root = do_firstNonNegativeRadix_int(x, 0, N, desc);
  if (in_place) {
    IntegerVector out(x);
    R_xlen_t a = desc ? 0 : root;
    R_xlen_t b = desc ? root : N;
    for (R_xlen_t i = a; i < b; ++i) {
      out[i] = 0;
    }
    return out;
  } else {
    IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (R_xlen_t i = 0; i < N; ++i) {
      bool zero = desc ? i < root : i >= root;
      out[i] = zero ? 0 : x[i];
    }
    return out;
  }
}

// Significantly slower to combine pmax and pmin!!

// [[Rcpp::export(rng = false)]]
IntegerVector do_pmax0_bitwise(IntegerVector x, int nThread = 1) {
  R_xlen_t N = x.size();
  R_xlen_t j = 0;
  while (j < N && x[j] >= 0) {
    ++j;
  }
  if (j == N) {
    return x;
  }

  IntegerVector out = no_init(N);

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    if (i < j) {
      out[i] = x[i];
      continue;
    }
    // https://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
    int xi = x[i];
    int r = xi - (xi & (xi >> (sizeof(int) * CHAR_BIT - 1)));
    out[i] = r;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
IntegerVector do_pmin0_bitwise(IntegerVector x, int nThread = 1) {
  R_xlen_t N = x.size();
  R_xlen_t j = 0;
  while (j < N && x[j] <= 0) {
    ++j;
  }
  if (j == N) {
    return x;
  }
  IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {

    // https://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
    int xi = x[i];
    if (i < j) {
      out[i] = xi;
      continue;
    }
    int r = (xi & (xi >> (sizeof(int) * CHAR_BIT - 1)));
    out[i] = r;
  }
  return out;
}



// [[Rcpp::export(rng = false)]]
SEXP do_pminpmax(SEXP X, SEXP Y,
                 bool do_min = false,
                 bool in_place = false,
                 bool keep_nas = false,
                 bool dbl_ok = false,
                 bool swap_xy = false,
                 int nThread = 1) {
  if (TYPEOF(X) == INTSXP && TYPEOF(Y) == INTSXP) {
    IntegerVector x = X;
    IntegerVector y = Y;
    R_xlen_t N = x.length();
    R_xlen_t Ny = y.length();
    if (N != Ny && Ny != 1) {
      stop("Lengths differ."); // # nocov
    }
    const bool y_lenN = Ny == N;
    const int y0 = y[0];
    if (in_place) {
      for (R_xlen_t i = 0; i < N; ++i) {
        int xi = x[i];
        int a = (y_lenN) ? y[i] : y0;
        bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
        int res = (choose_xi) ? xi : a;
        x[i] = res;
      }
      return x;
    } else {
      IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
      for (R_xlen_t i = 0; i < N; ++i) {
        int xi = x[i];
        int a = (y_lenN) ? y[i] : y0;
        bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
        int res = (choose_xi) ? xi : a;
        out[i] = res;
      }
      return out;
    }
  }

  if (TYPEOF(X) == REALSXP && TYPEOF(Y) == REALSXP) {
    DoubleVector x = X;
    DoubleVector y = Y;
    R_xlen_t N = x.length();
    R_xlen_t Ny = y.length();
    if (N != Ny && Ny != 1) {
      stop("Lengths differ."); // # nocov
    }
    const double y0 = y[0];
    if (in_place) {
      for (R_xlen_t i = 0; i < N; ++i) {
        double xi = x[i];
        double a = (Ny == N) ? y[i] : y0;
        bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
        double res = (choose_xi) ? xi : a;
        x[i] = res;
      }
      return x;
    } else {
      DoubleVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
      for (R_xlen_t i = 0; i < N; ++i) {
        double xi = x[i];
        double a = (Ny == N) ? y[i] : y0;
        bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
        double res = (choose_xi) ? xi : a;
        out[i] = res;
      }
      return out;
    }
  }

  if (TYPEOF(X) == INTSXP && TYPEOF(Y) == REALSXP) {
    IntegerVector x = X;
    DoubleVector y = Y;
    R_xlen_t N = x.length();
    R_xlen_t Ny = y.length();
    if (N != Ny && Ny != 1) {
      stop("Lengths differ."); // # nocov
    }
    const double y0 = y[0];
    bool int_out = (Ny == 1) && do_is_safe2int(y0);
    if (int_out) {
      if (in_place) {
        for (R_xlen_t i = 0; i < N; ++i) {
          int xi = x[i];
          int a = y0; // Ny == N not possible for now
          bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
          int res = (choose_xi) ? xi : a;
          x[i] = res;
        }
        return x;
      } else {
        IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
        for (R_xlen_t i = 0; i < N; ++i) {
          int xi = x[i];
          int a = y0;
          bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
          int res = (choose_xi) ? xi : a;
          out[i] = res;
        }
        return out;
      }
    } else {
      if (in_place) {
        stop("Internal error(do_pminpmax): INT x DBL y INPLACE"); // # nocov
      }
      DoubleVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
      for (R_xlen_t i = 0; i < N; ++i) {
        double xi = x[i];
        double a = (Ny == N) ? y[i] : y0;
        bool choose_xi = (keep_nas && (xi == NA_INTEGER)) || (do_min xor (xi > a));
        double res = (choose_xi) ? xi : a;
        out[i] = res;
      }
      message("Output is double");
      return out;
    }
  }
  stop("Internal error: unreachable pminmax.cpp"); // # nocov
  return R_NilValue; // # nocov
}


