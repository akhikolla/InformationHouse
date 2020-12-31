#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _GeneralizedUmatrix_addRowWiseC(SEXP, SEXP);
extern SEXP _GeneralizedUmatrix_Delta3DWeightsC(SEXP, SEXP);
extern SEXP _GeneralizedUmatrix_trainstepC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GeneralizedUmatrix_addRowWiseC",     (DL_FUNC) &_GeneralizedUmatrix_addRowWiseC,     2},
    {"_GeneralizedUmatrix_Delta3DWeightsC", (DL_FUNC) &_GeneralizedUmatrix_Delta3DWeightsC, 2},
    {"_GeneralizedUmatrix_trainstepC",      (DL_FUNC) &_GeneralizedUmatrix_trainstepC,      8},
    {NULL, NULL, 0}
};

void R_init_GeneralizedUmatrix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
/* 
output of call
tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
*/
