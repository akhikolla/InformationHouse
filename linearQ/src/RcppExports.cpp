// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// convHullConstrSum
DataFrame convHullConstrSum(const NumericVector& Px1, const NumericVector& Px2, const NumericVector& Qx1, const NumericVector& Qx2, double thd);
RcppExport SEXP _linearQ_convHullConstrSum(SEXP Px1SEXP, SEXP Px2SEXP, SEXP Qx1SEXP, SEXP Qx2SEXP, SEXP thdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type Px1(Px1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Px2(Px2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Qx1(Qx1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Qx2(Qx2SEXP);
    Rcpp::traits::input_parameter< double >::type thd(thdSEXP);
    rcpp_result_gen = Rcpp::wrap(convHullConstrSum(Px1, Px2, Qx1, Qx2, thd));
    return rcpp_result_gen;
END_RCPP
}
// mStat
double mStat(NumericVector data, int type);
RcppExport SEXP _linearQ_mStat(SEXP dataSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(mStat(data, type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_linearQ_convHullConstrSum", (DL_FUNC) &_linearQ_convHullConstrSum, 5},
    {"_linearQ_mStat", (DL_FUNC) &_linearQ_mStat, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_linearQ(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}