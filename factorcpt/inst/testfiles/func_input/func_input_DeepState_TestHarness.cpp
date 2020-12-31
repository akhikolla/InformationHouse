// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// func_input_DeepState_TestHarness_generation.cpp and func_input_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix func_input(NumericMatrix coef, NumericMatrix sgn);

TEST(factorcpt_deepstate_test,func_input_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix coef  = RcppDeepState_NumericMatrix();
  qs::c_qsave(coef,"/home/akhila/fuzzer_packages/fuzzedpackages/factorcpt/inst/testfiles/func_input/inputs/coef.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "coef values: "<< coef << std::endl;
  NumericMatrix sgn  = RcppDeepState_NumericMatrix();
  qs::c_qsave(sgn,"/home/akhila/fuzzer_packages/fuzzedpackages/factorcpt/inst/testfiles/func_input/inputs/sgn.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "sgn values: "<< sgn << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    func_input(coef,sgn);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}