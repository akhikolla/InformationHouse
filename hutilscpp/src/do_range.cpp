#include "cpphutils.h"


DoubleVector do_range_int(IntegerVector x, int halt_if_min, int halt_if_max) {
  R_xlen_t n = x.size();

  R_xlen_t which_max = n - 1;
  int maxi = x[n - 1];
  while (maxi == NA_INTEGER) {
    --which_max;
    maxi = x[which_max];
  }

  int mini = x[0];
  R_xlen_t which_min = 0;
  while (mini == NA_INTEGER) {
    ++which_min;
    mini = x[which_min];
  }

  const bool do_halt = halt_if_min < halt_if_max;
  if (do_halt) {
    maxi = x[which_min];
    which_max = which_min + 0;
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    int xi = x[i];
    if (xi == NA_INTEGER) {
      continue;
    }
    if (xi < mini) {
      mini = xi;
      which_min = i;
    } else {
      if (xi > maxi) {
        maxi = xi;
        which_max = i;
      }
    }
    if (do_halt) {
      if (mini <= halt_if_min ||
          maxi >= halt_if_max) {
        break;
      }
    }
  }
  DoubleVector out(4); // need to use double in case x is long
  ++which_min;
  ++which_max;

  out[0] = mini;
  out[1] = maxi;
  out[2] = which_min;
  out[3] = which_max;

  return out;
}


// [[Rcpp::export]]
DoubleVector do_range_dbl_simple(DoubleVector x) {
  R_xlen_t n = x.size();
  double mini = x[0];
  double maxi = x[0];
  R_xlen_t which_min = 0;
  R_xlen_t which_max = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i];
    if (xi < mini) {
      mini = xi;
      which_min = i;
    } else {
      if (xi > maxi) {
        maxi = xi;
        which_max = i;
      }
    }
  }
  ++which_min;
  ++which_max;
  DoubleVector out(4);
  out[0] = mini;
  out[1] = maxi;
  out[2] = which_min;
  out[3] = which_max;
  return out;
}



DoubleVector do_range_dbl(DoubleVector x, double halt_if_min, double halt_if_max) {
  R_xlen_t n = x.size();
  R_xlen_t which_max = n - 1;
  double maxi = x[n - 1];

  while ((R_IsNA(maxi) || R_IsNaN(maxi)) && which_max) {
    --which_max;
    maxi = x[which_max];
  }

  double mini = x[0];
  R_xlen_t which_min = 0;
  while ((R_IsNA(mini) || R_IsNaN(mini)) && which_min < n) {
    ++which_min;
    mini = x[which_min];
  }

  const bool do_halt = halt_if_min < halt_if_max;
  if (do_halt) {
    maxi = x[which_min];
    which_max = which_min + 0;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i];
    if (R_IsNA(xi) || R_IsNaN(xi)) {
      continue;
    }
    if (xi < mini) {
      mini = xi;
      which_min = i;
    } else {
      if (xi > maxi) {
        maxi = xi;
        which_max = i;
      }
    }
    if (do_halt) {
      if (mini <= halt_if_min ||
          maxi >= halt_if_max) {
        break;
      }
    }
  }
  ++which_min;
  ++which_max;
  DoubleVector out(4);
  out[0] = mini;
  out[1] = maxi;
  out[2] = which_min;
  out[3] = which_max;
  return out;

}

/*
/// [[Rcpp::export(rng = false)]]
IntegerVector do_prange2_int(IntegerVector x, int nThread = 1) {
  R_xlen_t N = x.length();
  IntegerVector out(2);
  if (N <= 1) {
    out[0] = (N == 0) ? INT_MAX : x[0];
    out[1] = (N == 0) ? INT_MIN : x[0];
    return out;
  }
  int xmax = x[0];
  int xmin = x[0];
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(min : xmin) reduction(max : xmax)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    if (xmin == 1 && xmax == 1000000000) {
      continue;
    }
    int xi = x[i];
    if (xi <= xmin) {
      xmin = xi;
    } else if (xi >= xmax) {
      xmax = xi;
    }
  }
  out[0] = xmin;
  out[1] = xmax;
  return out;
}


*/

