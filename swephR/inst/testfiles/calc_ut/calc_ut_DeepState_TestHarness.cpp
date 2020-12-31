// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// calc_ut_DeepState_TestHarness_generation.cpp and calc_ut_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::List calc_ut(Rcpp::NumericVector jd_ut, Rcpp::IntegerVector ipl, int iflag);

TEST(swephR_deepstate_test,calc_ut_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector jd_ut  = RcppDeepState_NumericVector();
  qs::c_qsave(jd_ut,"/home/akhila/fuzzer_packages/fuzzedpackages/swephR/inst/testfiles/calc_ut/inputs/jd_ut.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "jd_ut values: "<< jd_ut << std::endl;
  IntegerVector ipl  = RcppDeepState_IntegerVector();
  qs::c_qsave(ipl,"/home/akhila/fuzzer_packages/fuzzedpackages/swephR/inst/testfiles/calc_ut/inputs/ipl.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ipl values: "<< ipl << std::endl;
  IntegerVector iflag(1);
  iflag[0]  = RcppDeepState_int();
  qs::c_qsave(iflag,"/home/akhila/fuzzer_packages/fuzzedpackages/swephR/inst/testfiles/calc_ut/inputs/iflag.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "iflag values: "<< iflag << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    calc_ut(jd_ut,ipl,iflag[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}