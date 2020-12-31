#include "cpphutils.h"

namespace impl {
template <int RTYPE>
bool is_constant_(const Vector<RTYPE>& x, int nThread)
{
  R_xlen_t N = x.length();
  if (N <= 1) {
    return true;
  }

  if (nThread <= 1) {
    for (R_xlen_t i = 1; i < N; ++i) {
      if (x[i] != x[0]) {
        return false;
      }
    }
    return true;
  }

  R_xlen_t n_neq = 0;

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : n_neq)
#endif
  for (R_xlen_t i = 1; i < N; ++i) {
    n_neq += x[i] != x[0];
  }
  return n_neq == 0;
}

template <int RTYPE>
R_xlen_t isntConstant_(const Vector<RTYPE>& x) {
  R_xlen_t N = x.length();
  for (R_xlen_t i = 1; i < N; ++i) {
    if (x[i] != x[0]) {
      return i + 1;
    }
  }
  return 0;
}

}

// [[Rcpp::export(rng = false)]]
bool all_na_real(DoubleVector x, int nThread = 1) {
  R_xlen_t N = x.length();
  R_xlen_t n_na = ISNAN(x[0]);
  if (n_na) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : n_na)
#endif
    for (R_xlen_t i = 1; i < N; ++i) {
      n_na += ISNAN(x[i]);
    }
  }
  return n_na == N;
}

// [[Rcpp::export(rng = false)]]
R_xlen_t do_isntConstant_dbl(DoubleVector x) {
  R_xlen_t N = x.length();
  R_xlen_t n_na = ISNAN(x[0]);
  double x0 = x[0];
  for (R_xlen_t i = 1; i < N; ++i) {
    if (ISNAN(x[i])) {
        if (!n_na) {
          return i + 1;
        }
    } else {
      if (n_na) {
        return i + 1;
      }
      if (x[i] != x0) {
        return i + 1;
      }
    }
  }
  return 0;
}


// [[Rcpp::export(rng = false)]]
bool do_is_constant(SEXP x, int nThread = 1) {
  switch(TYPEOF(x)) {
  case NILSXP:
    return true;
  case LGLSXP:
    return impl::is_constant_(as<LogicalVector>(x), nThread);
  case INTSXP:
    return impl::is_constant_(as<IntegerVector>(x), nThread);
  case REALSXP:
    return impl::is_constant_(as<DoubleVector>(x), nThread);
  case STRSXP:
    return impl::is_constant_(as<CharacterVector>(x), nThread);
  case RAWSXP:
    return impl::is_constant_(as<RawVector>(x), nThread);
  }
  return false; // # nocov
}


// [[Rcpp::export(rng = false)]]
R_xlen_t do_isntConstant(SEXP x) {
  switch(TYPEOF(x)) {
  case NILSXP:
    return 0;
  case LGLSXP:
    return impl::isntConstant_(as<LogicalVector>(x));
  case INTSXP:
    return impl::isntConstant_(as<IntegerVector>(x));
  case REALSXP:
    return do_isntConstant_dbl(as<DoubleVector>(x));
  case STRSXP:
    return impl::isntConstant_(as<CharacterVector>(x));
  case RAWSXP:
    return impl::isntConstant_(as<RawVector>(x));
  }
  return 0; // # nocov
}
