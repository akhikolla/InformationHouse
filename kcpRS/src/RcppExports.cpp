// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getScatterMatrix
List getScatterMatrix(NumericMatrix II_, NumericMatrix X_, NumericMatrix H_);
RcppExport SEXP _kcpRS_getScatterMatrix(SEXP II_SEXP, SEXP X_SEXP, SEXP H_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type II_(II_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H_(H_SEXP);
    rcpp_result_gen = Rcpp::wrap(getScatterMatrix(II_, X_, H_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kcpRS_getScatterMatrix", (DL_FUNC) &_kcpRS_getScatterMatrix, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_kcpRS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}