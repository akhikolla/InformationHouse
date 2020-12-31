#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix half2mat(NumericVector X);

TEST(propr_deepstate_test,half2mat_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector X  = RcppDeepState_NumericVector();
  qs::c_qsave(X,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/half2mat/inputs/X.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "X values: "<< X << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    half2mat(X);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
