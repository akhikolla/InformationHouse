#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix rhoRcpp(NumericMatrix X, NumericMatrix lr, const int ivar);

TEST(propr_deepstate_test,rhoRcpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix X  = RcppDeepState_NumericMatrix();
  qs::c_qsave(X,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/rhoRcpp/inputs/X.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "X values: "<< X << std::endl;
  NumericMatrix lr  = RcppDeepState_NumericMatrix();
  qs::c_qsave(lr,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/rhoRcpp/inputs/lr.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lr values: "<< lr << std::endl;
  std::ofstream ivar_stream;
  int ivar  = RcppDeepState_int();
  ivar_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/rhoRcpp/inputs/ivar");
  ivar_stream << ivar;
  std::cout << "ivar values: "<< ivar << std::endl;
  ivar_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    rhoRcpp(X,lr,ivar);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
