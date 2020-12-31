#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix lr2phs(NumericMatrix lr);

TEST(propr_deepstate_test,lr2phs_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix lr  = RcppDeepState_NumericMatrix();
  qs::c_qsave(lr,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/lr2phs/inputs/lr.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lr values: "<< lr << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    lr2phs(lr);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
