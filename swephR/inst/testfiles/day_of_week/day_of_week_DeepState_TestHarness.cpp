// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// day_of_week_DeepState_TestHarness_generation.cpp and day_of_week_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::IntegerVector day_of_week(Rcpp::NumericVector jd);

TEST(swephR_deepstate_test,day_of_week_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector jd  = RcppDeepState_NumericVector();
  qs::c_qsave(jd,"/home/akhila/fuzzer_packages/fuzzedpackages/swephR/inst/testfiles/day_of_week/inputs/jd.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "jd values: "<< jd << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    day_of_week(jd);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}