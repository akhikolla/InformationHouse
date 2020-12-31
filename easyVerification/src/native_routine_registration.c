#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP easyVerification_EnsRocaCpp(SEXP, SEXP);
extern SEXP easyVerification_rankCpp(SEXP);
extern SEXP easyVerification_rankEnsCpp(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"easyVerification_EnsRocaCpp", (DL_FUNC) &easyVerification_EnsRocaCpp, 2},
    {"easyVerification_rankCpp",    (DL_FUNC) &easyVerification_rankCpp,    1},
    {"easyVerification_rankEnsCpp", (DL_FUNC) &easyVerification_rankEnsCpp, 1},
    {NULL, NULL, 0}
};

void R_init_easyVerification(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
