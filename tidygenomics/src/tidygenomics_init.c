#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _tidygenomics_cluster_interval(SEXP, SEXP, SEXP);
extern SEXP _tidygenomics_sort_indices(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_tidygenomics_cluster_interval", (DL_FUNC) &_tidygenomics_cluster_interval, 3},
    {"_tidygenomics_sort_indices",     (DL_FUNC) &_tidygenomics_sort_indices,     1},
    {NULL, NULL, 0}
};

void R_init_tidygenomics(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
