// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// parse_names
DataFrame parse_names(CharacterVector names);
RcppExport SEXP _humaniformat_parse_names(SEXP namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type names(namesSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_names(names));
    return rcpp_result_gen;
END_RCPP
}
// format_reverse
CharacterVector format_reverse(CharacterVector names);
RcppExport SEXP _humaniformat_format_reverse(SEXP namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type names(namesSEXP);
    rcpp_result_gen = Rcpp::wrap(format_reverse(names));
    return rcpp_result_gen;
END_RCPP
}
// format_period
CharacterVector format_period(CharacterVector names);
RcppExport SEXP _humaniformat_format_period(SEXP namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type names(namesSEXP);
    rcpp_result_gen = Rcpp::wrap(format_period(names));
    return rcpp_result_gen;
END_RCPP
}
// get_
CharacterVector get_(CharacterVector names, int element);
RcppExport SEXP _humaniformat_get_(SEXP namesSEXP, SEXP elementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type names(namesSEXP);
    Rcpp::traits::input_parameter< int >::type element(elementSEXP);
    rcpp_result_gen = Rcpp::wrap(get_(names, element));
    return rcpp_result_gen;
END_RCPP
}
// set_
CharacterVector set_(CharacterVector names, int element, String replacement);
RcppExport SEXP _humaniformat_set_(SEXP namesSEXP, SEXP elementSEXP, SEXP replacementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type names(namesSEXP);
    Rcpp::traits::input_parameter< int >::type element(elementSEXP);
    Rcpp::traits::input_parameter< String >::type replacement(replacementSEXP);
    rcpp_result_gen = Rcpp::wrap(set_(names, element, replacement));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_humaniformat_parse_names", (DL_FUNC) &_humaniformat_parse_names, 1},
    {"_humaniformat_format_reverse", (DL_FUNC) &_humaniformat_format_reverse, 1},
    {"_humaniformat_format_period", (DL_FUNC) &_humaniformat_format_period, 1},
    {"_humaniformat_get_", (DL_FUNC) &_humaniformat_get_, 2},
    {"_humaniformat_set_", (DL_FUNC) &_humaniformat_set_, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_humaniformat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}