// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// newey_west_pw_DeepState_TestHarness_generation.cpp and newey_west_pw_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List newey_west_pw(NumericMatrix hhat_mat, NumericMatrix xpxi_mat, NumericMatrix D_mat, int h);

TEST(lpirfs_deepstate_test,newey_west_pw_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix hhat_mat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(hhat_mat,"/home/akhila/fuzzer_packages/fuzzedpackages/lpirfs/inst/testfiles/newey_west_pw/inputs/hhat_mat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "hhat_mat values: "<< hhat_mat << std::endl;
  NumericMatrix xpxi_mat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(xpxi_mat,"/home/akhila/fuzzer_packages/fuzzedpackages/lpirfs/inst/testfiles/newey_west_pw/inputs/xpxi_mat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "xpxi_mat values: "<< xpxi_mat << std::endl;
  NumericMatrix D_mat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(D_mat,"/home/akhila/fuzzer_packages/fuzzedpackages/lpirfs/inst/testfiles/newey_west_pw/inputs/D_mat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "D_mat values: "<< D_mat << std::endl;
  IntegerVector h(1);
  h[0]  = RcppDeepState_int();
  qs::c_qsave(h,"/home/akhila/fuzzer_packages/fuzzedpackages/lpirfs/inst/testfiles/newey_west_pw/inputs/h.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "h values: "<< h << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    newey_west_pw(hhat_mat,xpxi_mat,D_mat,h[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}