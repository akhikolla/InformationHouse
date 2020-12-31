#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List weightedSlopeOneRM(NumericMatrix x);

TEST(rrecsys_deepstate_test,weightedSlopeOneRM_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/weightedSlopeOneRM/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    weightedSlopeOneRM(x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
