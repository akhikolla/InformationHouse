// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// rcpp_calculate_vic_DeepState_TestHarness_generation.cpp and rcpp_calculate_vic_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix rcpp_calculate_vic(NumericMatrix wic);

TEST(disclapmix_deepstate_test,rcpp_calculate_vic_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix wic  = RcppDeepState_NumericMatrix();
  qs::c_qsave(wic,"/home/akhila/fuzzer_packages/fuzzedpackages/disclapmix/inst/testfiles/rcpp_calculate_vic/inputs/wic.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "wic values: "<< wic << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    rcpp_calculate_vic(wic);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}