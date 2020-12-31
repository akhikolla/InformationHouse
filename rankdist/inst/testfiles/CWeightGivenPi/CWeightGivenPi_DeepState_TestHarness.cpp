#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector CWeightGivenPi(NumericVector r1, NumericVector r2);

TEST(rankdist_deepstate_test,CWeightGivenPi_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector r1  = RcppDeepState_NumericVector();
  qs::c_qsave(r1,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/CWeightGivenPi/inputs/r1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "r1 values: "<< r1 << std::endl;
  NumericVector r2  = RcppDeepState_NumericVector();
  qs::c_qsave(r2,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/CWeightGivenPi/inputs/r2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "r2 values: "<< r2 << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    CWeightGivenPi(r1,r2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
