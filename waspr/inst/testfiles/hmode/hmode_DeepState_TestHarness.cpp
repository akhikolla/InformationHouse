// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// hmode_DeepState_TestHarness_generation.cpp and hmode_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double hmode(NumericVector x, double cip);

TEST(waspr_deepstate_test,hmode_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/waspr/inst/testfiles/hmode/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector cip(1);
  cip[0]  = RcppDeepState_double();
  qs::c_qsave(cip,"/home/akhila/fuzzer_packages/fuzzedpackages/waspr/inst/testfiles/hmode/inputs/cip.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cip values: "<< cip << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    hmode(x,cip[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}