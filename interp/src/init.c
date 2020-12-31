#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP interp_inHull(SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_interpDeltri(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_interpShull(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_left(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_nearestNeighbours(SEXP, SEXP);
extern SEXP interp_on(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_onHull(SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_partDerivGrid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_partDerivPoints(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP interp_shullDeltri(SEXP, SEXP);
extern SEXP interp_triFind(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"interp_inHull",            (DL_FUNC) &interp_inHull,             4},
    {"interp_interpDeltri",      (DL_FUNC) &interp_interpDeltri,       6},
    {"interp_interpShull",       (DL_FUNC) &interp_interpShull,        8},
    {"interp_left",              (DL_FUNC) &interp_left,               7},
    {"interp_nearestNeighbours", (DL_FUNC) &interp_nearestNeighbours,  2},
    {"interp_on",                (DL_FUNC) &interp_on,                 7},
    {"interp_onHull",            (DL_FUNC) &interp_onHull,             4},
    {"interp_partDerivGrid",     (DL_FUNC) &interp_partDerivGrid,     12},
    {"interp_partDerivPoints",   (DL_FUNC) &interp_partDerivPoints,   12},
    {"interp_shullDeltri",       (DL_FUNC) &interp_shullDeltri,        2},
    {"interp_triFind",           (DL_FUNC) &interp_triFind,            8},
    {NULL, NULL, 0}
};

void R_init_interp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
