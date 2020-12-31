// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// c_MM_L2_SQ_WASS_D_DeepState_TestHarness_generation.cpp and c_MM_L2_SQ_WASS_D_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double c_MM_L2_SQ_WASS_D(NumericMatrix MM);

TEST(HistDAWass_deepstate_test,c_MM_L2_SQ_WASS_D_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix MM  = RcppDeepState_NumericMatrix();
  qs::c_qsave(MM,"/home/akhila/fuzzer_packages/fuzzedpackages/HistDAWass/inst/testfiles/c_MM_L2_SQ_WASS_D/inputs/MM.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "MM values: "<< MM << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    c_MM_L2_SQ_WASS_D(MM);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}