#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP landscapeR_assignValues_cpp(SEXP, SEXP, SEXP);
extern SEXP landscapeR_contigCells_cpp(SEXP, SEXP, SEXP);
extern SEXP landscapeR_indexTranspose_cpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"landscapeR_assignValues_cpp",   (DL_FUNC) &landscapeR_assignValues_cpp,   3},
  {"landscapeR_contigCells_cpp",    (DL_FUNC) &landscapeR_contigCells_cpp,    3},
  {"landscapeR_indexTranspose_cpp", (DL_FUNC) &landscapeR_indexTranspose_cpp, 3},
  {NULL, NULL, 0}
};

void R_init_landscapeR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
