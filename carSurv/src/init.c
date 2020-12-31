#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _carSurv_weightedCovarRcpp(SEXP, SEXP, SEXP);
extern SEXP _carSurv_weightedCovarRcppN(SEXP, SEXP, SEXP);
extern SEXP _carSurv_weightedVarRcpp(SEXP, SEXP);
extern SEXP _carSurv_weightedVarRcppN(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_carSurv_weightedCovarRcpp",  (DL_FUNC) &_carSurv_weightedCovarRcpp,  3},
  {"_carSurv_weightedCovarRcppN", (DL_FUNC) &_carSurv_weightedCovarRcppN, 3},
  {"_carSurv_weightedVarRcpp",    (DL_FUNC) &_carSurv_weightedVarRcpp,    2},
  {"_carSurv_weightedVarRcppN",   (DL_FUNC) &_carSurv_weightedVarRcppN,   2},
  {NULL, NULL, 0}
};

void R_init_carSurv(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
