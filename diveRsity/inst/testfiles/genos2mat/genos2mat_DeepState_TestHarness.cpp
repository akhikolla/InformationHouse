// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// genos2mat_DeepState_TestHarness_generation.cpp and genos2mat_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix genos2mat(NumericMatrix mat, IntegerVector ip, NumericVector na);

TEST(diveRsity_deepstate_test,genos2mat_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix mat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(mat,"/home/akhila/fuzzer_packages/fuzzedpackages/diveRsity/inst/testfiles/genos2mat/inputs/mat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "mat values: "<< mat << std::endl;
  IntegerVector ip  = RcppDeepState_IntegerVector();
  qs::c_qsave(ip,"/home/akhila/fuzzer_packages/fuzzedpackages/diveRsity/inst/testfiles/genos2mat/inputs/ip.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ip values: "<< ip << std::endl;
  NumericVector na  = RcppDeepState_NumericVector();
  qs::c_qsave(na,"/home/akhila/fuzzer_packages/fuzzedpackages/diveRsity/inst/testfiles/genos2mat/inputs/na.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "na values: "<< na << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    genos2mat(mat,ip,na);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}