#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP xyz_absolute_covariates(SEXP, SEXP);
extern SEXP xyz_absolute_covariates_pairs(SEXP, SEXP, SEXP);
extern SEXP xyz_calculate_residuals(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_calculate_xbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_clean_all_effects(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_clean_pairs(SEXP);
extern SEXP xyz_colsum_index(SEXP, SEXP);
extern SEXP xyz_create_lambda_sequence(SEXP, SEXP, SEXP);
extern SEXP xyz_equalpairs(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_estimate_background_interaction_frequency(SEXP, SEXP, SEXP);
extern SEXP xyz_find_strongest_pairs(SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_gaussiglmnet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_interaction_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_interaction_search_low_level(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_iterate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_naive_interaction_search(SEXP, SEXP, SEXP);
extern SEXP xyz_order_vector(SEXP, SEXP);
extern SEXP xyz_prod_matrix_vector(SEXP, SEXP);
extern SEXP xyz_projected_equal_pairs(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_sample_int_replace(SEXP, SEXP);
extern SEXP xyz_sample_uniform(SEXP, SEXP);
extern SEXP xyz_scale_intr(SEXP, SEXP, SEXP);
extern SEXP xyz_scan_intr_effects(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_scan_main_effects(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_soft_threshold(SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_sort_using_order_intmat(SEXP, SEXP);
extern SEXP xyz_sort_using_order_numvec(SEXP, SEXP);
extern SEXP xyz_translate_to_binary(SEXP, SEXP);
extern SEXP xyz_update_intr_final(SEXP, SEXP);
extern SEXP xyz_update_intr_vars(SEXP, SEXP, SEXP, SEXP);
extern SEXP xyz_warm_start(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"xyz_absolute_covariates",                       (DL_FUNC) &xyz_absolute_covariates,                        2},
  {"xyz_absolute_covariates_pairs",                 (DL_FUNC) &xyz_absolute_covariates_pairs,                  3},
  {"xyz_calculate_residuals",                       (DL_FUNC) &xyz_calculate_residuals,                        9},
  {"xyz_calculate_xbeta",                           (DL_FUNC) &xyz_calculate_xbeta,                            9},
  {"xyz_clean_all_effects",                         (DL_FUNC) &xyz_clean_all_effects,                          5},
  {"xyz_clean_pairs",                               (DL_FUNC) &xyz_clean_pairs,                                1},
  {"xyz_colsum_index",                              (DL_FUNC) &xyz_colsum_index,                               2},
  {"xyz_create_lambda_sequence",                    (DL_FUNC) &xyz_create_lambda_sequence,                     3},
  {"xyz_equalpairs",                                (DL_FUNC) &xyz_equalpairs,                                 5},
  {"xyz_estimate_background_interaction_frequency", (DL_FUNC) &xyz_estimate_background_interaction_frequency,  3},
  {"xyz_find_strongest_pairs",                      (DL_FUNC) &xyz_find_strongest_pairs,                       4},
  {"xyz_gaussiglmnet",                              (DL_FUNC) &xyz_gaussiglmnet,                               9},
  {"xyz_interaction_search",                        (DL_FUNC) &xyz_interaction_search,                         6},
  {"xyz_interaction_search_low_level",              (DL_FUNC) &xyz_interaction_search_low_level,               5},
  {"xyz_iterate",                                   (DL_FUNC) &xyz_iterate,                                   14},
  {"xyz_naive_interaction_search",                  (DL_FUNC) &xyz_naive_interaction_search,                   3},
  {"xyz_order_vector",                              (DL_FUNC) &xyz_order_vector,                               2},
  {"xyz_prod_matrix_vector",                        (DL_FUNC) &xyz_prod_matrix_vector,                         2},
  {"xyz_projected_equal_pairs",                     (DL_FUNC) &xyz_projected_equal_pairs,                      5},
  {"xyz_sample_int_replace",                        (DL_FUNC) &xyz_sample_int_replace,                         2},
  {"xyz_sample_uniform",                            (DL_FUNC) &xyz_sample_uniform,                             2},
  {"xyz_scale_intr",                                (DL_FUNC) &xyz_scale_intr,                                 3},
  {"xyz_scan_intr_effects",                         (DL_FUNC) &xyz_scan_intr_effects,                         12},
  {"xyz_scan_main_effects",                         (DL_FUNC) &xyz_scan_main_effects,                         10},
  {"xyz_soft_threshold",                            (DL_FUNC) &xyz_soft_threshold,                             4},
  {"xyz_sort_using_order_intmat",                   (DL_FUNC) &xyz_sort_using_order_intmat,                    2},
  {"xyz_sort_using_order_numvec",                   (DL_FUNC) &xyz_sort_using_order_numvec,                    2},
  {"xyz_translate_to_binary",                       (DL_FUNC) &xyz_translate_to_binary,                        2},
  {"xyz_update_intr_final",                         (DL_FUNC) &xyz_update_intr_final,                          2},
  {"xyz_update_intr_vars",                          (DL_FUNC) &xyz_update_intr_vars,                           4},
  {"xyz_warm_start",                                (DL_FUNC) &xyz_warm_start,                                 5},
  {NULL, NULL, 0}
};

void R_init_xyz(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
