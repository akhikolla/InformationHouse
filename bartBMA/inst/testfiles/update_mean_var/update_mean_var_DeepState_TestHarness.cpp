// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// update_mean_var_DeepState_TestHarness_generation.cpp and update_mean_var_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector update_mean_var(NumericMatrix tree_table, NumericMatrix tree_matrix, NumericVector resids, double a);

TEST(bartBMA_deepstate_test,update_mean_var_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix tree_table  = RcppDeepState_NumericMatrix();
  qs::c_qsave(tree_table,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/update_mean_var/inputs/tree_table.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "tree_table values: "<< tree_table << std::endl;
  NumericMatrix tree_matrix  = RcppDeepState_NumericMatrix();
  qs::c_qsave(tree_matrix,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/update_mean_var/inputs/tree_matrix.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "tree_matrix values: "<< tree_matrix << std::endl;
  NumericVector resids  = RcppDeepState_NumericVector();
  qs::c_qsave(resids,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/update_mean_var/inputs/resids.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "resids values: "<< resids << std::endl;
  NumericVector a(1);
  a[0]  = RcppDeepState_double();
  qs::c_qsave(a,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/update_mean_var/inputs/a.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "a values: "<< a << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    update_mean_var(tree_table,tree_matrix,resids,a[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}