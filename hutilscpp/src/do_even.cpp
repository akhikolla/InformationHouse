#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
LogicalVector do_are_even(IntegerVector x, DoubleVector y, int wb = 0, int nThread = 1) {
  R_xlen_t N = x.size();
  R_xlen_t M = y.size();
  const bool is_int = N > 0;

  if (is_int) {
    LogicalVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (R_xlen_t i = 0; i < N; ++i) {
      if (x[i] == NA_INTEGER) {
        out[i] = NA_LOGICAL;
        continue;
      }
      out[i] = !(x[i] % 2);
    }
    return out;
  } else {
    LogicalVector out = no_init(M);
    int wc = (wb > 1) ? (wb - 1) : N;
    for (int i = 0; i < M; ++i) {
      double yi = y[i];
      if (i >= wc) {
        if (!R_finite(yi)) {
          out[i] = NA_LOGICAL;
        } else {
          out[i] = std::fmod(yi, 2) == 0;
        }
        continue;
      }
      int yii = yi;
      out[i] = !(yii % 2);
    }
    return out;
  }
}


// [[Rcpp::export]]
IntegerVector do_which_even(IntegerVector x, DoubleVector y, int wb = 0) {
  R_xlen_t Norig = x.size();
  R_xlen_t Morig = y.size();
  if (Morig >= INT_MAX || Norig >= INT_MAX) {
    stop("Internal error: long vectors are not supported."); // # nocov
  }
  int N = (int)Norig;
  int M = (int)Morig;
  int wc = (wb > 1) ? (wb - 1) : N;
  const bool is_int = N > 0;
  std::vector<int> out1;
  out1.reserve(N / 2);
  if (is_int) {
    for (int i = 0; i < N; ++i) {
      // NA integer is even
      if ((x[i] % 2) == 0) {
        int j = i + 1;  // indexing
        out1.push_back(j);
      }
    }
  } else {
    for (int i = 0; i < M; ++i) {
      double yi0 = y[i];
      if (i >= wc && !R_finite(yi0)) {
        continue;
      }
      int yi = yi0;
      if ((yi % 2) == 0) {
        int j = i + 1;  // indexing
        out1.push_back(j);
      }
    }
  }
  return wrap(out1);
}
