// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// eucdist2_DeepState_TestHarness_generation.cpp and eucdist2_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix eucdist2(NumericVector x1, NumericVector y1, NumericVector x2, NumericVector y2, double eps);

TEST(gear_deepstate_test,eucdist2_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x1  = RcppDeepState_NumericVector();
  qs::c_qsave(x1,"/home/akhila/fuzzer_packages/fuzzedpackages/gear/inst/testfiles/eucdist2/inputs/x1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x1 values: "<< x1 << std::endl;
  NumericVector y1  = RcppDeepState_NumericVector();
  qs::c_qsave(y1,"/home/akhila/fuzzer_packages/fuzzedpackages/gear/inst/testfiles/eucdist2/inputs/y1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "y1 values: "<< y1 << std::endl;
  NumericVector x2  = RcppDeepState_NumericVector();
  qs::c_qsave(x2,"/home/akhila/fuzzer_packages/fuzzedpackages/gear/inst/testfiles/eucdist2/inputs/x2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x2 values: "<< x2 << std::endl;
  NumericVector y2  = RcppDeepState_NumericVector();
  qs::c_qsave(y2,"/home/akhila/fuzzer_packages/fuzzedpackages/gear/inst/testfiles/eucdist2/inputs/y2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "y2 values: "<< y2 << std::endl;
  NumericVector eps(1);
  eps[0]  = RcppDeepState_double();
  qs::c_qsave(eps,"/home/akhila/fuzzer_packages/fuzzedpackages/gear/inst/testfiles/eucdist2/inputs/eps.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "eps values: "<< eps << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    eucdist2(x1,y1,x2,y2,eps[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}