#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
bool is_altrep(SEXP x) {
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 5, 0)
  return ALTREP(x);
#else
  return false;
#endif
}

// [[Rcpp::export(rng = false)]]
IntegerVector allocate0_int(R_xlen_t N, int nThread = 1) {
  IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    out[i] = 0;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
DoubleVector allocate0_dbl(R_xlen_t N, int nThread = 1) {
  DoubleVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    out[i] = 0;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
IntegerVector allocate0_except(R_xlen_t N, DoubleVector India, IntegerVector Victor, int nThread = 1) {
  IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    out[i] = 0;
  }
  R_xlen_t tn = India.length();
  R_xlen_t vn = Victor.length();
  if (vn < 1 || tn < 1) {
    warning("Ignoring Victor."); // # nocov
    return out; // # nocov
  }
  for (R_xlen_t j = 0; j < tn; ++j) {
    R_xlen_t i = India[j];
    if (i < 0 || i >= N) {
      continue; // # nocov
    }
    int v = (tn == vn) ? Victor[j] : Victor[0];
    out[i] = v;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
IntegerVector allocate_with_root(R_xlen_t N,
                                 int a,
                                 R_xlen_t r,
                                 bool left,
                                 bool do_pmin,
                                 int nThread = 1) {
  IntegerVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    bool post = i >= r;
    if (post xor left) {
      out[i] = a;
      continue;
    }
    if (do_pmin) {
      out[i] = post ? (r - i) : (i - r);
    } else {
      out[i] = post ? (i - r) : (r - i);
    }
    out[i] += a;
  }
  return out;
}

/*
/./ [[Rcpp::export(rng = false)]]
IntegerVector all_integers(int nThread = 1) {
  IntegerVector out = no_init(4294967296);
#pragma omp parallel for num_threads(nThread)
  for (unsigned int i = 0; i < 4294967295; ++i) {
    int outi = static_cast<int>(i);
    out[i] = outi;
  }
  out[4294967295] = -1;
  return out;
}
*/


