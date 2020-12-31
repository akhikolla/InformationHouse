// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// expHet_DeepState_TestHarness_generation.cpp and expHet_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector expHet(NumericMatrix af);

TEST(diveRsity_deepstate_test,expHet_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix af  = RcppDeepState_NumericMatrix();
  qs::c_qsave(af,"/home/akhila/fuzzer_packages/fuzzedpackages/diveRsity/inst/testfiles/expHet/inputs/af.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "af values: "<< af << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    expHet(af);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}