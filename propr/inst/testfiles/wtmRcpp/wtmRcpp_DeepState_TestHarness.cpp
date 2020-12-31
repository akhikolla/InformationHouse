#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double wtmRcpp(NumericVector x, NumericVector w);

TEST(propr_deepstate_test,wtmRcpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/wtmRcpp/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector w  = RcppDeepState_NumericVector();
  qs::c_qsave(w,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/wtmRcpp/inputs/w.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "w values: "<< w << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    wtmRcpp(x,w);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
