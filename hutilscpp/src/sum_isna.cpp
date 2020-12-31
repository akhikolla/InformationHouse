#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
R_xlen_t sum_isna_int(IntegerVector x, int nThread = 1) {
  R_xlen_t n = x.size();
  R_xlen_t out = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    out += x[i] == NA_INTEGER;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t sum_isna_dbl(DoubleVector x, int nThread = 1) {
  R_xlen_t n = x.size();
  R_xlen_t out = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    out += ISNAN(x[i]);
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t sum_isna_complx(ComplexVector x, int nThread = 1) {
  R_xlen_t n = x.size();
  R_xlen_t out = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    Rcomplex v = COMPLEX_ELT(x, i);
    if (ISNAN(v.r) || ISNAN(v.i)) {
      out += 1;
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t sum_isna_char(CharacterVector x, int nThread = 1) {
  R_xlen_t n = x.size();
  R_xlen_t out = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    out += x[i] == NA_STRING;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t sum_isfalse(LogicalVector x, int nThread = 1) {
  R_xlen_t n = x.size();
  R_xlen_t out = n;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    bool xi = x[i];
    out -= xi;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t sum_isna_logi(LogicalVector x, int nThread = 1) {
  R_xlen_t n = x.size();
  R_xlen_t out = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    out += x[i] == NA_LOGICAL;
  }
  return out;
}
