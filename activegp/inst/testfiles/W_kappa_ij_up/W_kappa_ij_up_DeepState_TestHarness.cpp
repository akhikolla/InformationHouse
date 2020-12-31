// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// W_kappa_ij_up_DeepState_TestHarness_generation.cpp and W_kappa_ij_up_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

void W_kappa_ij_up(NumericMatrix W, NumericMatrix design, NumericVector theta, int i1, int i2, int start, int ct);

TEST(activegp_deepstate_test,W_kappa_ij_up_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix W  = RcppDeepState_NumericMatrix();
  qs::c_qsave(W,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/W.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "W values: "<< W << std::endl;
  NumericMatrix design  = RcppDeepState_NumericMatrix();
  qs::c_qsave(design,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/design.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "design values: "<< design << std::endl;
  NumericVector theta  = RcppDeepState_NumericVector();
  qs::c_qsave(theta,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/theta.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "theta values: "<< theta << std::endl;
  IntegerVector i1(1);
  i1[0]  = RcppDeepState_int();
  qs::c_qsave(i1,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/i1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "i1 values: "<< i1 << std::endl;
  IntegerVector i2(1);
  i2[0]  = RcppDeepState_int();
  qs::c_qsave(i2,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/i2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "i2 values: "<< i2 << std::endl;
  IntegerVector start(1);
  start[0]  = RcppDeepState_int();
  qs::c_qsave(start,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/start.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "start values: "<< start << std::endl;
  IntegerVector ct(1);
  ct[0]  = RcppDeepState_int();
  qs::c_qsave(ct,"/home/akhila/fuzzer_packages/fuzzedpackages/activegp/inst/testfiles/W_kappa_ij_up/inputs/ct.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ct values: "<< ct << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    W_kappa_ij_up(W,design,theta,i1[0],i2[0],start[0],ct[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}