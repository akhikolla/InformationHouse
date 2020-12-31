#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BondValuation_CppPrevDate(SEXP);
extern SEXP _BondValuation_CppSuccDate(SEXP);
extern SEXP _BondValuation_Date_LDM(SEXP);
extern SEXP _BondValuation_DayDiff(SEXP);
extern SEXP _BondValuation_DaysInMonth(SEXP);
extern SEXP _BondValuation_DaysInYear(SEXP);
extern SEXP _BondValuation_DIST(SEXP);
extern SEXP _BondValuation_FirstMatch(SEXP);
extern SEXP _BondValuation_LDM(SEXP);
extern SEXP _BondValuation_leap(SEXP);
extern SEXP _BondValuation_LeapDayInside(SEXP);
extern SEXP _BondValuation_NumToDate(SEXP);
extern SEXP _BondValuation_PayCalc(SEXP);
extern SEXP _BondValuation_sumC(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BondValuation_CppPrevDate",   (DL_FUNC) &_BondValuation_CppPrevDate,   1},
    {"_BondValuation_CppSuccDate",   (DL_FUNC) &_BondValuation_CppSuccDate,   1},
    {"_BondValuation_Date_LDM",      (DL_FUNC) &_BondValuation_Date_LDM,      1},
    {"_BondValuation_DayDiff",       (DL_FUNC) &_BondValuation_DayDiff,       1},
    {"_BondValuation_DaysInMonth",   (DL_FUNC) &_BondValuation_DaysInMonth,   1},
    {"_BondValuation_DaysInYear",    (DL_FUNC) &_BondValuation_DaysInYear,    1},
    {"_BondValuation_DIST",          (DL_FUNC) &_BondValuation_DIST,          1},
    {"_BondValuation_FirstMatch",    (DL_FUNC) &_BondValuation_FirstMatch,    1},
    {"_BondValuation_LDM",           (DL_FUNC) &_BondValuation_LDM,           1},
    {"_BondValuation_leap",          (DL_FUNC) &_BondValuation_leap,          1},
    {"_BondValuation_LeapDayInside", (DL_FUNC) &_BondValuation_LeapDayInside, 1},
    {"_BondValuation_NumToDate",     (DL_FUNC) &_BondValuation_NumToDate,     1},
    {"_BondValuation_PayCalc",       (DL_FUNC) &_BondValuation_PayCalc,       1},
    {"_BondValuation_sumC",          (DL_FUNC) &_BondValuation_sumC,          1},
    {NULL, NULL, 0}
};

void R_init_BondValuation(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
