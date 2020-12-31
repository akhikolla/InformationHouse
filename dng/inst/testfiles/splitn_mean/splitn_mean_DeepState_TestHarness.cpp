// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// splitn_mean_DeepState_TestHarness_generation.cpp and splitn_mean_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector splitn_mean(NumericVector mu, NumericVector sigma, NumericVector lmd);

TEST(dng_deepstate_test,splitn_mean_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector mu  = RcppDeepState_NumericVector();
  qs::c_qsave(mu,"/home/akhila/fuzzer_packages/fuzzedpackages/dng/inst/testfiles/splitn_mean/inputs/mu.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "mu values: "<< mu << std::endl;
  NumericVector sigma  = RcppDeepState_NumericVector();
  qs::c_qsave(sigma,"/home/akhila/fuzzer_packages/fuzzedpackages/dng/inst/testfiles/splitn_mean/inputs/sigma.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "sigma values: "<< sigma << std::endl;
  NumericVector lmd  = RcppDeepState_NumericVector();
  qs::c_qsave(lmd,"/home/akhila/fuzzer_packages/fuzzedpackages/dng/inst/testfiles/splitn_mean/inputs/lmd.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lmd values: "<< lmd << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    splitn_mean(mu,sigma,lmd);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}