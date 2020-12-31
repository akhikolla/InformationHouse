#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector CapplyB(NumericVector us, NumericMatrix X, NumericVector mu);

TEST(robFitConGraph_deepstate_test,CapplyB_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector us  = RcppDeepState_NumericVector();
  qs::c_qsave(us,"/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/CapplyB/inputs/us.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "us values: "<< us << std::endl;
  NumericMatrix X  = RcppDeepState_NumericMatrix();
  qs::c_qsave(X,"/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/CapplyB/inputs/X.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "X values: "<< X << std::endl;
  NumericVector mu  = RcppDeepState_NumericVector();
  qs::c_qsave(mu,"/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/CapplyB/inputs/mu.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "mu values: "<< mu << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    CapplyB(us,X,mu);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
