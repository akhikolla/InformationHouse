#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP lclGWAS_alphaEst(SEXP, SEXP);
extern SEXP lclGWAS_betaEst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lclGWAS_varEst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"lclGWAS_alphaEst", (DL_FUNC) &lclGWAS_alphaEst, 2},
    {"lclGWAS_betaEst",  (DL_FUNC) &lclGWAS_betaEst,  8},
    {"lclGWAS_varEst",   (DL_FUNC) &lclGWAS_varEst,   8},
    {NULL, NULL, 0}
};

void R_init_lclGWAS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

