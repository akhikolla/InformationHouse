#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _cna_C_append_intList(SEXP, SEXP);
extern SEXP _cna_C_charList2string(SEXP, SEXP, SEXP, SEXP);
extern SEXP _cna_C_checkHall_k(SEXP, SEXP);
extern SEXP _cna_C_checkHallsCondition(SEXP);
extern SEXP _cna_C_concat(SEXP, SEXP);
extern SEXP _cna_C_conCov(SEXP, SEXP, SEXP);
extern SEXP _cna_C_conj_conCov(SEXP, SEXP, SEXP, SEXP);
extern SEXP _cna_C_conjScore(SEXP, SEXP);
extern SEXP _cna_C_countUniques(SEXP);
extern SEXP _cna_C_disj_conCov(SEXP, SEXP, SEXP, SEXP);
extern SEXP _cna_C_disj_contained(SEXP, SEXP, SEXP);
extern SEXP _cna_C_duplicatedMat(SEXP);
extern SEXP _cna_C_find_asf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _cna_C_find_first_false(SEXP);
extern SEXP _cna_C_hasSubsetIn(SEXP, SEXP);
extern SEXP _cna_C_hasSubsetInM(SEXP, SEXP);
extern SEXP _cna_C_hasSupersetIn(SEXP, SEXP, SEXP);
extern SEXP _cna_C_increment(SEXP, SEXP, SEXP, SEXP);
extern SEXP _cna_C_init_ii(SEXP, SEXP);
extern SEXP _cna_C_intList_minimal_old(SEXP, SEXP, SEXP);
extern SEXP _cna_C_is_submodel(SEXP, SEXP, SEXP);
extern SEXP _cna_C_isSubsetOf(SEXP, SEXP);
extern SEXP _cna_C_max_which(SEXP);
extern SEXP _cna_C_mconcat(SEXP, SEXP, SEXP);
extern SEXP _cna_C_minimal(SEXP, SEXP, SEXP);
extern SEXP _cna_C_minimal_old(SEXP, SEXP, SEXP);
extern SEXP _cna_C_mredund(SEXP, SEXP);
extern SEXP _cna_C_recCharList2char(SEXP, SEXP, SEXP, SEXP);
extern SEXP _cna_C_redund(SEXP);
extern SEXP _cna_C_relist_Char(SEXP, SEXP);
extern SEXP _cna_C_relist_Int(SEXP, SEXP);
extern SEXP _cna_C_relist_List(SEXP, SEXP);
extern SEXP _cna_C_relist_Log(SEXP, SEXP);
extern SEXP _cna_C_relist_Num(SEXP, SEXP);
extern SEXP _cna_C_rowSubsetColAnys(SEXP, SEXP);
extern SEXP _cna_C_selectCols(SEXP, SEXP);
extern SEXP _cna_C_set_lim(SEXP, SEXP);
extern SEXP _cna_C_subsetMin(SEXP, SEXP);
extern SEXP _cna_C_uniqueCombs(SEXP, SEXP);
extern SEXP _cna_C_uniqueMat(SEXP);
extern SEXP _cna_intList_equal(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_cna_C_append_intList",      (DL_FUNC) &_cna_C_append_intList,      2},
    {"_cna_C_charList2string",     (DL_FUNC) &_cna_C_charList2string,     4},
    {"_cna_C_checkHall_k",         (DL_FUNC) &_cna_C_checkHall_k,         2},
    {"_cna_C_checkHallsCondition", (DL_FUNC) &_cna_C_checkHallsCondition, 1},
    {"_cna_C_concat",              (DL_FUNC) &_cna_C_concat,              2},
    {"_cna_C_conCov",              (DL_FUNC) &_cna_C_conCov,              3},
    {"_cna_C_conj_conCov",         (DL_FUNC) &_cna_C_conj_conCov,         4},
    {"_cna_C_conjScore",           (DL_FUNC) &_cna_C_conjScore,           2},
    {"_cna_C_countUniques",        (DL_FUNC) &_cna_C_countUniques,        1},
    {"_cna_C_disj_conCov",         (DL_FUNC) &_cna_C_disj_conCov,         4},
    {"_cna_C_disj_contained",      (DL_FUNC) &_cna_C_disj_contained,      3},
    {"_cna_C_duplicatedMat",       (DL_FUNC) &_cna_C_duplicatedMat,       1},
    {"_cna_C_find_asf",            (DL_FUNC) &_cna_C_find_asf,            7},
    {"_cna_C_find_first_false",    (DL_FUNC) &_cna_C_find_first_false,    1},
    {"_cna_C_hasSubsetIn",         (DL_FUNC) &_cna_C_hasSubsetIn,         2},
    {"_cna_C_hasSubsetInM",        (DL_FUNC) &_cna_C_hasSubsetInM,        2},
    {"_cna_C_hasSupersetIn",       (DL_FUNC) &_cna_C_hasSupersetIn,       3},
    {"_cna_C_increment",           (DL_FUNC) &_cna_C_increment,           4},
    {"_cna_C_init_ii",             (DL_FUNC) &_cna_C_init_ii,             2},
    {"_cna_C_intList_minimal_old", (DL_FUNC) &_cna_C_intList_minimal_old, 3},
    {"_cna_C_is_submodel",         (DL_FUNC) &_cna_C_is_submodel,         3},
    {"_cna_C_isSubsetOf",          (DL_FUNC) &_cna_C_isSubsetOf,          2},
    {"_cna_C_max_which",           (DL_FUNC) &_cna_C_max_which,           1},
    {"_cna_C_mconcat",             (DL_FUNC) &_cna_C_mconcat,             3},
    {"_cna_C_minimal",             (DL_FUNC) &_cna_C_minimal,             3},
    {"_cna_C_minimal_old",         (DL_FUNC) &_cna_C_minimal_old,         3},
    {"_cna_C_mredund",             (DL_FUNC) &_cna_C_mredund,             2},
    {"_cna_C_recCharList2char",    (DL_FUNC) &_cna_C_recCharList2char,    4},
    {"_cna_C_redund",              (DL_FUNC) &_cna_C_redund,              1},
    {"_cna_C_relist_Char",         (DL_FUNC) &_cna_C_relist_Char,         2},
    {"_cna_C_relist_Int",          (DL_FUNC) &_cna_C_relist_Int,          2},
    {"_cna_C_relist_List",         (DL_FUNC) &_cna_C_relist_List,         2},
    {"_cna_C_relist_Log",          (DL_FUNC) &_cna_C_relist_Log,          2},
    {"_cna_C_relist_Num",          (DL_FUNC) &_cna_C_relist_Num,          2},
    {"_cna_C_rowSubsetColAnys",    (DL_FUNC) &_cna_C_rowSubsetColAnys,    2},
    {"_cna_C_selectCols",          (DL_FUNC) &_cna_C_selectCols,          2},
    {"_cna_C_set_lim",             (DL_FUNC) &_cna_C_set_lim,             2},
    {"_cna_C_subsetMin",           (DL_FUNC) &_cna_C_subsetMin,           2},
    {"_cna_C_uniqueCombs",         (DL_FUNC) &_cna_C_uniqueCombs,         2},
    {"_cna_C_uniqueMat",           (DL_FUNC) &_cna_C_uniqueMat,           1},
    {"_cna_intList_equal",         (DL_FUNC) &_cna_intList_equal,         2},
    {NULL, NULL, 0}
};

void R_init_cna(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
