#include "cpphutils.h"

constexpr bool ON_TWOS_COMPLEMENT = (-1 == ~0);

// [[Rcpp::export]]
LogicalVector do_divisible(IntegerVector x, int d, int nThread = 1) {
  R_xlen_t N = x.length();
  LogicalVector out = no_init(N);

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    out[i] = !(x[i] % d);
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
LogicalVector do_divisible2(IntegerVector x, int nThread = 1) {
  R_xlen_t N = x.length();
  // # nocov start
  if (!ON_TWOS_COMPLEMENT) {
    stop("Unexpected arch.");
  }
  // # nocov end

  LogicalVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) schedule(static)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    uint_fast32_t ui = (uint_fast32_t)x[i];
    out[i] = !(ui & 1U);
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector do_divisible16(IntegerVector x, int nThread = 1) {
  R_xlen_t N = x.length();
  LogicalVector out = no_init(N);

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    unsigned int xi = (unsigned int)x[i];
    out[i] = !(xi % 16);
  }
  return out;
}

