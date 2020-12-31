// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// make_density_multilevel_DeepState_TestHarness_generation.cpp and make_density_multilevel_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::NumericVector make_density_multilevel(Rcpp::NumericVector ordered, Rcpp::NumericVector interval);

TEST(imagerExtra_deepstate_test,make_density_multilevel_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector ordered  = RcppDeepState_NumericVector();
  qs::c_qsave(ordered,"/home/akhila/fuzzer_packages/fuzzedpackages/imagerExtra/inst/testfiles/make_density_multilevel/inputs/ordered.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ordered values: "<< ordered << std::endl;
  NumericVector interval  = RcppDeepState_NumericVector();
  qs::c_qsave(interval,"/home/akhila/fuzzer_packages/fuzzedpackages/imagerExtra/inst/testfiles/make_density_multilevel/inputs/interval.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "interval values: "<< interval << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    make_density_multilevel(ordered,interval);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}