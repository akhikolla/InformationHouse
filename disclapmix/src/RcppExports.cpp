// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// disclapglm_linkfun
NumericVector disclapglm_linkfun(NumericVector mu);
RcppExport SEXP _disclapmix_disclapglm_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(disclapglm_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// disclapglm_linkinv
NumericVector disclapglm_linkinv(NumericVector eta);
RcppExport SEXP _disclapmix_disclapglm_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(disclapglm_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// disclapglm_mu_eta
NumericVector disclapglm_mu_eta(NumericVector eta);
RcppExport SEXP _disclapmix_disclapglm_mu_eta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(disclapglm_mu_eta(eta));
    return rcpp_result_gen;
END_RCPP
}
// disclapglm_varfunc
NumericVector disclapglm_varfunc(NumericVector mu);
RcppExport SEXP _disclapmix_disclapglm_varfunc(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(disclapglm_varfunc(mu));
    return rcpp_result_gen;
END_RCPP
}
// disclapglm_loglikeh
double disclapglm_loglikeh(double mu, double y);
RcppExport SEXP _disclapmix_disclapglm_loglikeh(SEXP muSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(disclapglm_loglikeh(mu, y));
    return rcpp_result_gen;
END_RCPP
}
// disclapglm_deviance
double disclapglm_deviance(NumericVector y, NumericVector mu, NumericVector wt);
RcppExport SEXP _disclapmix_disclapglm_deviance(SEXP ySEXP, SEXP muSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(disclapglm_deviance(y, mu, wt));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_create_design_matrix
IntegerMatrix rcpp_create_design_matrix(IntegerMatrix x, int clusters);
RcppExport SEXP _disclapmix_rcpp_create_design_matrix(SEXP xSEXP, SEXP clustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type clusters(clustersSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_create_design_matrix(x, clusters));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_create_new_weight_vector
NumericVector rcpp_create_new_weight_vector(NumericMatrix vic, int loci);
RcppExport SEXP _disclapmix_rcpp_create_new_weight_vector(SEXP vicSEXP, SEXP lociSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type vic(vicSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_create_new_weight_vector(vic, loci));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_create_response_vector
IntegerVector rcpp_create_response_vector(IntegerMatrix x, IntegerMatrix y);
RcppExport SEXP _disclapmix_rcpp_create_response_vector(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_create_response_vector(x, y));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_calculate_wic
NumericMatrix rcpp_calculate_wic(IntegerMatrix x, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP _disclapmix_rcpp_calculate_wic(SEXP xSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_calculate_wic(x, y, p, tau));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_calculate_vic
NumericMatrix rcpp_calculate_vic(NumericMatrix wic);
RcppExport SEXP _disclapmix_rcpp_calculate_vic(SEXP wicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type wic(wicSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_calculate_vic(wic));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_calculate_haplotype_probabilities
NumericVector rcpp_calculate_haplotype_probabilities(IntegerMatrix new_data, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP _disclapmix_rcpp_calculate_haplotype_probabilities(SEXP new_dataSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type new_data(new_dataSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_calculate_haplotype_probabilities(new_data, y, p, tau));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_calculate_haplotype_probabilities_se
NumericVector rcpp_calculate_haplotype_probabilities_se(IntegerMatrix new_data, IntegerMatrix y, NumericVector theta_clusters, NumericVector theta_loci, NumericMatrix vcov, NumericVector tau);
RcppExport SEXP _disclapmix_rcpp_calculate_haplotype_probabilities_se(SEXP new_dataSEXP, SEXP ySEXP, SEXP theta_clustersSEXP, SEXP theta_lociSEXP, SEXP vcovSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type new_data(new_dataSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_clusters(theta_clustersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_loci(theta_lociSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vcov(vcovSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_calculate_haplotype_probabilities_se(new_data, y, theta_clusters, theta_loci, vcov, tau));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_calculate_haplotype_probabilities_clusterwise
NumericMatrix rcpp_calculate_haplotype_probabilities_clusterwise(IntegerMatrix new_data, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP _disclapmix_rcpp_calculate_haplotype_probabilities_clusterwise(SEXP new_dataSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type new_data(new_dataSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_calculate_haplotype_probabilities_clusterwise(new_data, y, p, tau));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_simulate
IntegerMatrix rcpp_simulate(int nsim, IntegerMatrix y, NumericVector tau_cumsum, NumericMatrix disclap_parameters);
RcppExport SEXP _disclapmix_rcpp_simulate(SEXP nsimSEXP, SEXP ySEXP, SEXP tau_cumsumSEXP, SEXP disclap_parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau_cumsum(tau_cumsumSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type disclap_parameters(disclap_parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_simulate(nsim, y, tau_cumsum, disclap_parameters));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_find_haplotype_in_matrix
int rcpp_find_haplotype_in_matrix(const IntegerMatrix subpop, const IntegerVector h);
RcppExport SEXP _disclapmix_rcpp_find_haplotype_in_matrix(SEXP subpopSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type subpop(subpopSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_find_haplotype_in_matrix(subpop, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_disclapmix_disclapglm_linkfun", (DL_FUNC) &_disclapmix_disclapglm_linkfun, 1},
    {"_disclapmix_disclapglm_linkinv", (DL_FUNC) &_disclapmix_disclapglm_linkinv, 1},
    {"_disclapmix_disclapglm_mu_eta", (DL_FUNC) &_disclapmix_disclapglm_mu_eta, 1},
    {"_disclapmix_disclapglm_varfunc", (DL_FUNC) &_disclapmix_disclapglm_varfunc, 1},
    {"_disclapmix_disclapglm_loglikeh", (DL_FUNC) &_disclapmix_disclapglm_loglikeh, 2},
    {"_disclapmix_disclapglm_deviance", (DL_FUNC) &_disclapmix_disclapglm_deviance, 3},
    {"_disclapmix_rcpp_create_design_matrix", (DL_FUNC) &_disclapmix_rcpp_create_design_matrix, 2},
    {"_disclapmix_rcpp_create_new_weight_vector", (DL_FUNC) &_disclapmix_rcpp_create_new_weight_vector, 2},
    {"_disclapmix_rcpp_create_response_vector", (DL_FUNC) &_disclapmix_rcpp_create_response_vector, 2},
    {"_disclapmix_rcpp_calculate_wic", (DL_FUNC) &_disclapmix_rcpp_calculate_wic, 4},
    {"_disclapmix_rcpp_calculate_vic", (DL_FUNC) &_disclapmix_rcpp_calculate_vic, 1},
    {"_disclapmix_rcpp_calculate_haplotype_probabilities", (DL_FUNC) &_disclapmix_rcpp_calculate_haplotype_probabilities, 4},
    {"_disclapmix_rcpp_calculate_haplotype_probabilities_se", (DL_FUNC) &_disclapmix_rcpp_calculate_haplotype_probabilities_se, 6},
    {"_disclapmix_rcpp_calculate_haplotype_probabilities_clusterwise", (DL_FUNC) &_disclapmix_rcpp_calculate_haplotype_probabilities_clusterwise, 4},
    {"_disclapmix_rcpp_simulate", (DL_FUNC) &_disclapmix_rcpp_simulate, 4},
    {"_disclapmix_rcpp_find_haplotype_in_matrix", (DL_FUNC) &_disclapmix_rcpp_find_haplotype_in_matrix, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_disclapmix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}