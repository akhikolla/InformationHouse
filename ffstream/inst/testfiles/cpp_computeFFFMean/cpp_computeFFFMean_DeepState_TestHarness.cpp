// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// cpp_computeFFFMean_DeepState_TestHarness_generation.cpp and cpp_computeFFFMean_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double cpp_computeFFFMean(Rcpp::NumericVector x, double lambda);

TEST(ffstream_deepstate_test,cpp_computeFFFMean_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/ffstream/inst/testfiles/cpp_computeFFFMean/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector lambda(1);
  lambda[0]  = RcppDeepState_double();
  qs::c_qsave(lambda,"/home/akhila/fuzzer_packages/fuzzedpackages/ffstream/inst/testfiles/cpp_computeFFFMean/inputs/lambda.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lambda values: "<< lambda << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    cpp_computeFFFMean(x,lambda[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}