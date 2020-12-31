#include "cpphutils.h"

// returns a list of three elements (intended to be passed to an if statement
// immediately after so no NA_LOGICALs)
// whether to warn/error
// whether to warn
// messages

// # nocov start

// [[Rcpp::export(rng = false)]]
bool has_openmp() {
#ifdef _OPENMP
  return true;
#endif
  return false;
}


// [[Rcpp::export(rng = false)]]
List diagnose_omp(int threads_requested,
                  String msg_no_openmp,
                  String msg_threads_neg,
                  String msg_unknown_issues,
                  String msg_too_many_threads) {
  LogicalVector False(1);
  LogicalVector True(1);
  True[0] = TRUE;
  CharacterVector out(1);


#ifndef _OPENMP
  return List::create(False, False, out);
#endif

  int n_procs = 1;
#ifdef _OPENMP
n_procs = omp_get_num_procs();
#endif




  if (threads_requested > 0 && threads_requested <= n_procs) {
    out[0] = msg_threads_neg;
    return List::create(False, False, out);
  }

  if (threads_requested < 0) {
    out[0] = msg_threads_neg;
    return List::create(True, False, out);
  }
  if (threads_requested > n_procs) {
    out[0] = msg_too_many_threads;
    return List::create(True, False, out);
  }
  out[0] = msg_unknown_issues;

  return List::create(True, True, out);
}

// # nocov end

