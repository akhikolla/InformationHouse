#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SpecsVerification_auc_cpp(SEXP, SEXP);
extern SEXP _SpecsVerification_aucdiff_cpp(SEXP, SEXP, SEXP);
extern SEXP _SpecsVerification_dresscrps_cpp(SEXP, SEXP, SEXP);
extern SEXP _SpecsVerification_enscrps_cpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SpecsVerification_auc_cpp",       (DL_FUNC) &_SpecsVerification_auc_cpp,       2},
    {"_SpecsVerification_aucdiff_cpp",   (DL_FUNC) &_SpecsVerification_aucdiff_cpp,   3},
    {"_SpecsVerification_dresscrps_cpp", (DL_FUNC) &_SpecsVerification_dresscrps_cpp, 3},
    {"_SpecsVerification_enscrps_cpp",   (DL_FUNC) &_SpecsVerification_enscrps_cpp,   3},
    {NULL, NULL, 0}
};

void R_init_SpecsVerification(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
