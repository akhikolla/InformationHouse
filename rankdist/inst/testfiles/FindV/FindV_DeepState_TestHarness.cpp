#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix FindV(NumericMatrix obs, NumericVector pi0);

TEST(rankdist_deepstate_test,FindV_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix obs  = RcppDeepState_NumericMatrix();
  qs::c_qsave(obs,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/FindV/inputs/obs.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "obs values: "<< obs << std::endl;
  NumericVector pi0  = RcppDeepState_NumericVector();
  qs::c_qsave(pi0,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/FindV/inputs/pi0.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "pi0 values: "<< pi0 << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    FindV(obs,pi0);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
