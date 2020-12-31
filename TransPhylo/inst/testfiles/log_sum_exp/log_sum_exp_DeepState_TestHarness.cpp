// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// log_sum_exp_DeepState_TestHarness_generation.cpp and log_sum_exp_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double log_sum_exp(double u, double v);

TEST(TransPhylo_deepstate_test,log_sum_exp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector u(1);
  u[0]  = RcppDeepState_double();
  qs::c_qsave(u,"/home/akhila/fuzzer_packages/fuzzedpackages/TransPhylo/inst/testfiles/log_sum_exp/inputs/u.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "u values: "<< u << std::endl;
  NumericVector v(1);
  v[0]  = RcppDeepState_double();
  qs::c_qsave(v,"/home/akhila/fuzzer_packages/fuzzedpackages/TransPhylo/inst/testfiles/log_sum_exp/inputs/v.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "v values: "<< v << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    log_sum_exp(u[0],v[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}