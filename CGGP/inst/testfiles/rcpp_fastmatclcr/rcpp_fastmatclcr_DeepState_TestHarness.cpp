// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// rcpp_fastmatclcr_DeepState_TestHarness_generation.cpp and rcpp_fastmatclcr_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

void rcpp_fastmatclcr(NumericMatrix I, NumericVector w, NumericMatrix MSEmat, NumericVector S, int maxlevel);

TEST(CGGP_deepstate_test,rcpp_fastmatclcr_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix I  = RcppDeepState_NumericMatrix();
  qs::c_qsave(I,"/home/akhila/fuzzer_packages/fuzzedpackages/CGGP/inst/testfiles/rcpp_fastmatclcr/inputs/I.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "I values: "<< I << std::endl;
  NumericVector w  = RcppDeepState_NumericVector();
  qs::c_qsave(w,"/home/akhila/fuzzer_packages/fuzzedpackages/CGGP/inst/testfiles/rcpp_fastmatclcr/inputs/w.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "w values: "<< w << std::endl;
  NumericMatrix MSEmat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(MSEmat,"/home/akhila/fuzzer_packages/fuzzedpackages/CGGP/inst/testfiles/rcpp_fastmatclcr/inputs/MSEmat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "MSEmat values: "<< MSEmat << std::endl;
  NumericVector S  = RcppDeepState_NumericVector();
  qs::c_qsave(S,"/home/akhila/fuzzer_packages/fuzzedpackages/CGGP/inst/testfiles/rcpp_fastmatclcr/inputs/S.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "S values: "<< S << std::endl;
  IntegerVector maxlevel(1);
  maxlevel[0]  = RcppDeepState_int();
  qs::c_qsave(maxlevel,"/home/akhila/fuzzer_packages/fuzzedpackages/CGGP/inst/testfiles/rcpp_fastmatclcr/inputs/maxlevel.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "maxlevel values: "<< maxlevel << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    rcpp_fastmatclcr(I,w,MSEmat,S,maxlevel[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}