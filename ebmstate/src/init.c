#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "survS.h"
#include <Rversion.h>
#include "survproto.h"

/*
  The following symbols/expressions for .NAME have been omitted

    _ebmstate_convolute_semiMarkov
    _ebmstate_convolute_Markov

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void agmssurv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void Ccoxfit5a(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void Ccoxfit5b(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void Ccoxfit5c( void *, void *, void *, void *, void *);

extern void Cagfit5a(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void Cagfit5b(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void Cagfit5c(void *);


static const R_CMethodDef CEntries[] = {
    {"agmssurv", (DL_FUNC) &agmssurv, 20},
    {"Ccoxfit5a",   (DL_FUNC) &coxfit5_a, 20},
    {"Ccoxfit5b",   (DL_FUNC) &coxfit5_b, 19},
    {"Ccoxfit5c",   (DL_FUNC) &coxfit5_c,  5},
    {"Cagfit5a",    (DL_FUNC) &agfit5a,  20},
    {"Cagfit5b",    (DL_FUNC) &agfit5b,  19},
    {"Cagfit5c",    (DL_FUNC) &agfit5c,   1},
    {"Ccoxscho",    (DL_FUNC) &coxscho,    8},
    {"Ccoxscore",    (DL_FUNC) &coxscore,    10},
    {"Cagscore",    (DL_FUNC) &agscore,   10},
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP _ebmstate_convolute_Markov(SEXP, SEXP, SEXP,SEXP);
extern SEXP _ebmstate_convolute_semiMarkov(SEXP, SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_ebmstate_convolute_Markov",            (DL_FUNC) &_ebmstate_convolute_Markov,    4},
  {"_ebmstate_convolute_semiMarkov",    (DL_FUNC) &_ebmstate_convolute_semiMarkov,    3},
  {"Ccoxcount1",    (DL_FUNC) &coxcount1,    2},
  {"Ccoxcount2",    (DL_FUNC) &coxcount2,    4},
  {"Cagmart3",      (DL_FUNC) &agmart3,      6},
  {NULL, NULL, 0}
};

void R_init_ebmstate(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
