// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// d_gauss_cpp_DeepState_TestHarness_generation.cpp and d_gauss_cpp_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector d_gauss_cpp(NumericVector X, double x, double sigma);

TEST(activegp_deepstate_test,d_gauss_cpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector X  = RcppDeepState_NumericVector();
  qs::c_qsave(X,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/d_gauss_cpp/inputs/X.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "X values: "<< X << std::endl;
  NumericVector x(1);
  x[0]  = RcppDeepState_double();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/d_gauss_cpp/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector sigma(1);
  sigma[0]  = RcppDeepState_double();
  qs::c_qsave(sigma,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/d_gauss_cpp/inputs/sigma.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "sigma values: "<< sigma << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    d_gauss_cpp(X,x[0],sigma[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}