#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _foreSIGHT_latentX_calc_cpp(SEXP, SEXP, SEXP);
extern SEXP _foreSIGHT_Pstatus_WGEN_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _foreSIGHT_residualGenerator_cpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_foreSIGHT_latentX_calc_cpp",      (DL_FUNC) &_foreSIGHT_latentX_calc_cpp,      3},
    {"_foreSIGHT_Pstatus_WGEN_cpp",      (DL_FUNC) &_foreSIGHT_Pstatus_WGEN_cpp,      4},
    {"_foreSIGHT_residualGenerator_cpp", (DL_FUNC) &_foreSIGHT_residualGenerator_cpp, 2},
    {NULL, NULL, 0}
};

void R_init_foreSIGHT(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
