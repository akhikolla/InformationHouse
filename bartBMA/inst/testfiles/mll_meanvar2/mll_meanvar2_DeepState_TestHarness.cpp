// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// mll_meanvar2_DeepState_TestHarness_generation.cpp and mll_meanvar2_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double mll_meanvar2(double x, double x2, int n);

TEST(bartBMA_deepstate_test,mll_meanvar2_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x(1);
  x[0]  = RcppDeepState_double();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/mll_meanvar2/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector x2(1);
  x2[0]  = RcppDeepState_double();
  qs::c_qsave(x2,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/mll_meanvar2/inputs/x2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x2 values: "<< x2 << std::endl;
  IntegerVector n(1);
  n[0]  = RcppDeepState_int();
  qs::c_qsave(n,"/home/akhila/fuzzer_packages/fuzzedpackages/bartBMA/inst/testfiles/mll_meanvar2/inputs/n.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "n values: "<< n << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    mll_meanvar2(x[0],x2[0],n[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}