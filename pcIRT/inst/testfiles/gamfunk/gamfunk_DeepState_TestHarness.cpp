#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector gamfunk(NumericVector epsmat, NumericMatrix gammat);

TEST(pcIRT_deepstate_test,gamfunk_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector epsmat  = RcppDeepState_NumericVector();
  qs::c_qsave(epsmat,"/home/akhila/fuzzer_packages/fuzzedpackages/pcIRT/inst/testfiles/gamfunk/inputs/epsmat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "epsmat values: "<< epsmat << std::endl;
  NumericMatrix gammat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(gammat,"/home/akhila/fuzzer_packages/fuzzedpackages/pcIRT/inst/testfiles/gamfunk/inputs/gammat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "gammat values: "<< gammat << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    gamfunk(epsmat,gammat);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
