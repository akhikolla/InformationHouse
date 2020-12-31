// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// interpolate_C_DeepState_TestHarness_generation.cpp and interpolate_C_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector interpolate_C(NumericVector original_samples, int target_frequency);

TEST(AGread_deepstate_test,interpolate_C_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector original_samples  = RcppDeepState_NumericVector();
  qs::c_qsave(original_samples,"/home/akhila/fuzzer_packages/fuzzedpackages/AGread/inst/testfiles/interpolate_C/inputs/original_samples.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "original_samples values: "<< original_samples << std::endl;
  IntegerVector target_frequency(1);
  target_frequency[0]  = RcppDeepState_int();
  qs::c_qsave(target_frequency,"/home/akhila/fuzzer_packages/fuzzedpackages/AGread/inst/testfiles/interpolate_C/inputs/target_frequency.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "target_frequency values: "<< target_frequency << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    interpolate_C(original_samples,target_frequency[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}