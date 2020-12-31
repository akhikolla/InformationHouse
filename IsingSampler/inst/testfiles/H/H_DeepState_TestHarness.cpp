// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// H_DeepState_TestHarness_generation.cpp and H_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double H(NumericMatrix J, IntegerVector s, NumericVector h);

TEST(IsingSampler_deepstate_test,H_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix J  = RcppDeepState_NumericMatrix();
  qs::c_qsave(J,"/home/akhila/fuzzer_packages/fuzzedpackages/IsingSampler/inst/testfiles/H/inputs/J.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "J values: "<< J << std::endl;
  IntegerVector s  = RcppDeepState_IntegerVector();
  qs::c_qsave(s,"/home/akhila/fuzzer_packages/fuzzedpackages/IsingSampler/inst/testfiles/H/inputs/s.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "s values: "<< s << std::endl;
  NumericVector h  = RcppDeepState_NumericVector();
  qs::c_qsave(h,"/home/akhila/fuzzer_packages/fuzzedpackages/IsingSampler/inst/testfiles/H/inputs/h.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "h values: "<< h << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    H(J,s,h);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}