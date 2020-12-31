// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// ids_to_count_matrix_DeepState_TestHarness_generation.cpp and ids_to_count_matrix_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::IntegerMatrix ids_to_count_matrix(Rcpp::IntegerVector ids);

TEST(RcppCWB_deepstate_test,ids_to_count_matrix_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector ids  = RcppDeepState_IntegerVector();
  qs::c_qsave(ids,"/home/akhila/fuzzer_packages/fuzzedpackages/RcppCWB/inst/testfiles/ids_to_count_matrix/inputs/ids.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ids values: "<< ids << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    ids_to_count_matrix(ids);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}