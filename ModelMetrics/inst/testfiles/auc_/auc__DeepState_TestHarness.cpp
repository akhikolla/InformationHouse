// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// auc__DeepState_TestHarness_generation.cpp and auc__DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double auc_(NumericVector actual, NumericVector predicted);

TEST(ModelMetrics_deepstate_test,auc__test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector actual  = RcppDeepState_NumericVector();
  qs::c_qsave(actual,"/home/akhila/fuzzer_packages/fuzzedpackages/ModelMetrics/inst/testfiles/auc_/inputs/actual.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "actual values: "<< actual << std::endl;
  NumericVector predicted  = RcppDeepState_NumericVector();
  qs::c_qsave(predicted,"/home/akhila/fuzzer_packages/fuzzedpackages/ModelMetrics/inst/testfiles/auc_/inputs/predicted.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "predicted values: "<< predicted << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    auc_(actual,predicted);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}