// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// median_pillars_DeepState_TestHarness_generation.cpp and median_pillars_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix median_pillars(NumericVector arr3d);

TEST(autothresholdr_deepstate_test,median_pillars_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector arr3d  = RcppDeepState_NumericVector();
  qs::c_qsave(arr3d,"/home/akhila/fuzzer_packages/fuzzedpackages/autothresholdr/inst/testfiles/median_pillars/inputs/arr3d.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "arr3d values: "<< arr3d << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    median_pillars(arr3d);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}