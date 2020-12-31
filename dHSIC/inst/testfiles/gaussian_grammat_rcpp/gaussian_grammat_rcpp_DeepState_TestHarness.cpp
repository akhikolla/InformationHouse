// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// gaussian_grammat_rcpp_DeepState_TestHarness_generation.cpp and gaussian_grammat_rcpp_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix gaussian_grammat_rcpp(NumericMatrix x, double bandwidth, int n, int d);

TEST(dHSIC_deepstate_test,gaussian_grammat_rcpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/dHSIC/inst/testfiles/gaussian_grammat_rcpp/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector bandwidth(1);
  bandwidth[0]  = RcppDeepState_double();
  qs::c_qsave(bandwidth,"/home/akhila/fuzzer_packages/fuzzedpackages/dHSIC/inst/testfiles/gaussian_grammat_rcpp/inputs/bandwidth.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "bandwidth values: "<< bandwidth << std::endl;
  IntegerVector n(1);
  n[0]  = RcppDeepState_int();
  qs::c_qsave(n,"/home/akhila/fuzzer_packages/fuzzedpackages/dHSIC/inst/testfiles/gaussian_grammat_rcpp/inputs/n.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "n values: "<< n << std::endl;
  IntegerVector d(1);
  d[0]  = RcppDeepState_int();
  qs::c_qsave(d,"/home/akhila/fuzzer_packages/fuzzedpackages/dHSIC/inst/testfiles/gaussian_grammat_rcpp/inputs/d.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "d values: "<< d << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    gaussian_grammat_rcpp(x,bandwidth[0],n[0],d[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}