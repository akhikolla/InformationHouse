// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// int_randomShuffle_DeepState_TestHarness_generation.cpp and int_randomShuffle_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::IntegerVector int_randomShuffle(Rcpp::IntegerVector a);

TEST(BoostMLR_deepstate_test,int_randomShuffle_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector a  = RcppDeepState_IntegerVector();
  qs::c_qsave(a,"/home/akhila/fuzzer_packages/fuzzedpackages/BoostMLR/inst/testfiles/int_randomShuffle/inputs/a.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "a values: "<< a << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    int_randomShuffle(a);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}