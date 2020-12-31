// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// numeric_at_interpolator_DeepState_TestHarness_generation.cpp and numeric_at_interpolator_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector numeric_at_interpolator(NumericVector from, NumericVector to, NumericVector at, CharacterVector ease);

TEST(tweenr_deepstate_test,numeric_at_interpolator_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector from  = RcppDeepState_NumericVector();
  qs::c_qsave(from,"/home/akhila/fuzzer_packages/fuzzedpackages/tweenr/inst/testfiles/numeric_at_interpolator/inputs/from.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "from values: "<< from << std::endl;
  NumericVector to  = RcppDeepState_NumericVector();
  qs::c_qsave(to,"/home/akhila/fuzzer_packages/fuzzedpackages/tweenr/inst/testfiles/numeric_at_interpolator/inputs/to.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "to values: "<< to << std::endl;
  NumericVector at  = RcppDeepState_NumericVector();
  qs::c_qsave(at,"/home/akhila/fuzzer_packages/fuzzedpackages/tweenr/inst/testfiles/numeric_at_interpolator/inputs/at.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "at values: "<< at << std::endl;
  CharacterVector ease  = RcppDeepState_CharacterVector();
  qs::c_qsave(ease,"/home/akhila/fuzzer_packages/fuzzedpackages/tweenr/inst/testfiles/numeric_at_interpolator/inputs/ease.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ease values: "<< ease << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    numeric_at_interpolator(from,to,at,ease);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}